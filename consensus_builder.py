#!/usr/bin/env python3
"""

Usage example:
    ./consensus_builder.py --fasta ref.fa --vcf vars.vcf --fragment_length 1000 \
        --num_fragments 5 --out outdir --af_threshold 0.01

    To start at a given position:
    ./consensus_builder.py --fasta ref.fa --vcf vars.vcf --fragment_length 1000 \
        --num_fragments 5 --start_chr chr1 --start_pos 752721
"""

from __future__ import annotations
import argparse
import os
import sys
import gzip
from typing import List, Dict, Tuple, Optional


# Chromosome sort order helper

# Basic ordering
def chrom_sort_key(chrom: str):
    # Normalize common forms of chromosomes annotation
    c = chrom
    if c.startswith("chr"):
        c = c[3:]

    try:
        val = int(c)
        return (0, val)
    except:
        pass
    if c == "X":
        return (0, 23)
    if c == "Y":
        return (0, 24)
    if c in ("M", "MT"):
        return (0, 25)

    return (1, c)


# handeling I / O

def open_text(path: str):
    """Open text file; handle gzipped files"""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8")

def load_fasta(path: str) -> Dict[str, str]:
    """Load FASTA into dict {chrom: sequence} (uppercase)."""
    genome: Dict[str, List[str]] = {}
    cur = None
    chunks: List[str] = []
    with open_text(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if cur is not None:
                    genome[cur] = "".join(chunks)
                cur = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip().upper())
        if cur is not None:
            genome[cur] = "".join(chunks)
    return genome

def wrap_fasta(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

# VCF parsing 

def parse_info_af(info: str) -> List[float]:
    """
    Extract AF= values from INFO field (as floats).
    Returns [] if AF not present or malformed.
    """
    if not info or info == ".":
        return []
    for part in info.split(";"):
        if part.startswith("AF="):
            raw = part[3:]
            if raw == "":
                return []
            try:
                return [float(x) for x in raw.split(",")]
            except ValueError:
                return []
    return []

def compute_af_from_genotypes(genotypes: List[str], max_alt_guess: Optional[int] = None) -> List[float]:
    """
    Compute allele frequency per ALT allele (could be multi alleles) from a list of genotype strings for row.
    genotypes: list of sample genotype fields like "0/1:..." or "./.:..."
    Returns list of AFs for alleles 1..k (k inferred).
    """
    counts: Dict[int, int] = {}
    total = 0
    for g in genotypes:
        if not g or g == ".":
            continue
        gt = g.split(":")[0]
        if gt in (".", "./.", ".|."):
            continue
        # normalize phasing to '/'
        gt = gt.replace("|", "/")
        for a in gt.split("/"):
            if a == ".":
                continue
            try:
                ai = int(a)
            except:
                continue
            counts[ai] = counts.get(ai, 0) + 1
            total += 1
    if total == 0:
        return []
    # handeling multi alleles -_-
    max_idx = max([i for i in counts.keys() if i > 0] + ([max_alt_guess] if max_alt_guess else [0]) + [0])
    afs = []
    for alt_i in range(1, max_idx + 1):
        afs.append(counts.get(alt_i, 0) / total)
    return afs

def load_vcf(path: str) -> Tuple[List[str], Dict[str, List[Tuple[int,str,List[str],str,List[str]]]]]:
    """
    Parse VCF into:
      - samples: list of sample names
      - variants_by_chrom: dict chrom -> list of variant tuples
        variant tuple: (pos:int, ref:str, alts:List[str], info:str, genotypes:List[str])
    Notes:
      - stores variants in the order they appear; we will sort later.
      - supports gzipped VCF when filename ends with .gz (even though in task it said teh input would be unzipped, i just added this option just in case)
      - this would have been easier with pandas :(
    """
    samples: List[str] = []
    variants_by_chrom: Dict[str, List[Tuple[int,str,List[str],str,List[str]]]] = {}
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            # check needed column names from header
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").lstrip("#").split("\t")
                samples = cols[9:]
                continue
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom = cols[0]
            pos = int(cols[1])
            ref = cols[3]
            alt_field = cols[4]
            alts = [] if alt_field == "." else alt_field.split(",")
            info = cols[7]
            genotypes = cols[9:] if len(cols) > 9 else []
            variants_by_chrom.setdefault(chrom, []).append((pos, ref, alts, info, genotypes))
    return samples, variants_by_chrom


# Consensus part

def select_alt_by_af(gt_allele: int, afs: List[float], af_threshold: float) -> Optional[int]:
    """
    Decide whether to use given genotype allele (0 = ref, >=1 = ALT index).
    Returns:
        0  => use reference
        k>=1 => use ALT index k (1-based)
        None => skip (AF below threshold or missing)
    """
    if gt_allele == 0:
        return 0
    if not afs:
        return None
    alt_idx = gt_allele - 1
    if alt_idx < 0 or alt_idx >= len(afs):
        return None
    if afs[alt_idx] < af_threshold:
        return None
    return gt_allele

def apply_variants_to_reference(ref_seq: str,
                                chrom_variants: List[Tuple[int,str,List[str],str,List[str]]],
                                sample_index: int,
                                window_start: int,
                                window_len: int,
                                af_threshold: float) -> str:
    """
    Build consensus for [window_start .. window_start+window_len-1] on a single chromosome.
    chrom_variants: list of variants for this chromosome (pos, ref, alts, info, genotypes)
    sample_index: index in genotypes list for the sample
    """
    start0 = window_start - 1
    end0 = start0 + window_len
    seq_list = list(ref_seq[start0:end0])

    # iterate variants in this chrom and apply those in window 
    for pos, ref, alts, info, genotypes in chrom_variants:
        if pos < window_start or pos >= window_start + window_len:
            continue

        # Get GT for this sample
        if sample_index >= len(genotypes):
            continue
        gt_field = genotypes[sample_index].split(":")[0]
        if gt_field in (".", "./.", ".|."):
            continue

        # parse alleles (e.g., "0/1" -> [0,1])
        alleles = []
        for a in gt_field.replace("|", "/").split("/"):
            if a == ".":
                alleles.append(None)
            else:
                try:
                    alleles.append(int(a))
                except:
                    alleles.append(None)

        # AF: prefer INFO AF, else compute from genotypes row
        afs = parse_info_af(info)
        if not afs:
            afs = compute_af_from_genotypes(genotypes, max_alt_guess=len(alts))

        # apply alleles (choose first ALT allele found in genotype that passes AF threshold)
        # ALT wins rule:
        # If any allele in the genotype is ALT (>=1), pick the first ALT allele, if it passes AF threshold.

        # alt_choice = None
        # for allele in alleles:
        #     if allele is None:
        #         continue
        #     if allele >= 1:  # ALT detected
        #         alt_choice = allele
        #         break

        pos0 = pos - window_start
        if not (0 <= pos0 < len(seq_list)):
            continue

        # determine genotype ALT allele (first ALT >0 found)
        alt_allele = None
        for allele in alleles:
            if allele is not None and allele >= 1:
                alt_allele = allele
                break

        # case 1: no ALT allele in genotype → use REF
        if alt_allele is None:
            seq_list[pos0:pos0+len(ref)] = list(ref)
            continue

        # Parse AF list into numeric floats
        try:
            af_values = [float(x) for x in afs]
        except:
            af_values = []

        # AF index maps to ALT allele index: allele 1 → af_values[0]
        allele_index = alt_allele - 1

        # Check AF threshold
        if allele_index >= len(af_values) or af_values[allele_index] < af_threshold:
            # AF fails → use REF
            seq_list[pos0:pos0+len(ref)] = list(ref)
            continue

        # ALT passes AF threshold → apply ALT
        alt_seq = alts[allele_index]
        seq_list[pos0:pos0+len(ref)] = list(alt_seq)

    return "".join(seq_list)


# Fragment selection and sorting


def build_global_variant_list(variants_by_chrom: Dict[str, List[Tuple[int,str,List[str],str,List[str]]]]):
    """
    Build and return a sorted list of (chrom, pos) from variants_by_chrom.
    Sorting uses chrom_sort_key then position.
    Also returns a mapping chrom -> sorted variant list (by pos).
    """
    global_list: List[Tuple[str,int]] = []
    sorted_variants_by_chrom: Dict[str, List[Tuple[int,str,List[str],str,List[str]]]] = {}
    for chrom, vlist in variants_by_chrom.items():
        # sort variants in this chrom by position
        v_sorted = sorted(vlist, key=lambda x: x[0])
        sorted_variants_by_chrom[chrom] = v_sorted
        for pos, ref, alts, info, genotypes in v_sorted:
            global_list.append((chrom, pos))
    # sort list by chrom order then position
    global_list.sort(key=lambda x: (chrom_sort_key(x[0]), x[1]))
    return global_list, sorted_variants_by_chrom

def select_fragments_from_global(global_variants: List[Tuple[str,int]],
                                 sorted_variants_by_chrom: Dict[str, List[Tuple[int,str,List[str],str,List[str]]]],
                                 fragment_length: int,
                                 n_fragments: Optional[int],  # None means "all"
                                 start_chr: Optional[str],
                                 start_pos: Optional[int]) -> List[Tuple[str,int,int]]:
    """
    VCF-driven global fragment selection.
    - If start_chr & start_pos provided: start from it.
    - Else start from the first variant in global_variants.
    - For each fragment: take variant at current index i -> (chrom,pos).
        fragment = (chrom, pos, pos + fragment_length - 1)
      then advance i to first variant whose (chrom,pos) is outside this fragment window.
    - Continue until n_fragments reached or variants exhausted.
    - If n_fragments is None, treat as "all" (cover all variants).
    """
    frags: List[Tuple[str,int,int]] = []
    if not global_variants:
        return frags

    # find starting index in global_variants
    start_idx = 0
    if start_chr and start_pos:
        # find first global variant with chrom == start_chr and pos >= start_pos
        found = False
        for idx, (c, p) in enumerate(global_variants):
            if c == start_chr and p >= start_pos:
                start_idx = idx
                found = True
                break
        if not found:
            # no variant at or after start_pos on that chrom -> return empty
            return frags
    else:
        start_idx = 0

    i = start_idx
    total = len(global_variants)
    needed = None if n_fragments is None else int(n_fragments)

    while i < total and (needed is None or len(frags) < needed):
        chrom, pos = global_variants[i]
        frag_start = pos
        frag_end = pos + fragment_length - 1
        frags.append((chrom, frag_start, frag_end))

        # advance i to first variant strictly > frag_end or different chrom
        j = i + 1
        while j < total:
            c2, p2 = global_variants[j]
            # if different chromosome, it's outside this window (since window is chrom-specific)
            if c2 != chrom:
                break
            if p2 > frag_end:
                break
            j += 1
        i = j

    return frags

# Main CLI logic


def main():
    p = argparse.ArgumentParser(description="Pure-Python consensus generator (VCF-driven fragments).")
    p.add_argument("--fasta", required=True, help="Reference FASTA (plain or .gz)")
    p.add_argument("--vcf", required=True, help="VCF file (plain or .gz)")
    p.add_argument("--out", default=".", help="Output directory (default: current dir)")
    p.add_argument("--fragment_length", type=int, required=True, help="Fragment length (bp)")
    p.add_argument("--num_fragments", required=True, help="Number of fragments or 'all'")
    p.add_argument("--af_threshold", type=float, default=0.0, help="AF threshold (default 0.0)")
    p.add_argument("--sample", default=None, help="Sample name to run (default: all samples)")
    p.add_argument("--start_chr", default=None, help="Optional start chromosome (e.g. chr1)")
    p.add_argument("--start_pos", type=int, default=None, help="Optional start position (1-based)")
    args = p.parse_args()

    # Prepare output dir
    os.makedirs(args.out, exist_ok=True)

    # Load fasta
    print("Loading FASTA...", file=sys.stderr)
    genome = load_fasta(args.fasta)

    # Load VCF
    print("Loading VCF...", file=sys.stderr)
    samples, variants_by_chrom = load_vcf(args.vcf)

    if not samples:
        print("No samples found in VCF header (#CHROM line). Exiting.", file=sys.stderr)
        sys.exit(1)

    # If a specific sample in input, validate
    if args.sample:
        if args.sample not in samples:
            print(f"Sample {args.sample} not found in VCF samples.", file=sys.stderr)
            sys.exit(1)
        sample_list = [args.sample]
    else:
        sample_list = samples

    # build global sorted variant list and per-chrom sorted lists
    global_variants, sorted_variants_by_chrom = build_global_variant_list(variants_by_chrom)

    # Determine requested number of fragments
    if args.num_fragments == "all":
        n_fragments = None
    else:
        try:
            n_fragments = int(args.num_fragments)
            if n_fragments <= 0:
                raise ValueError()
        except ValueError:
            print("num_fragments must be a positive integer or 'all'.", file=sys.stderr)
            sys.exit(1)

    # Select fragments 
    fragments = select_fragments_from_global(
        global_variants=global_variants,
        sorted_variants_by_chrom=sorted_variants_by_chrom,
        fragment_length=args.fragment_length,
        n_fragments=n_fragments,
        start_chr=args.start_chr,
        start_pos=args.start_pos
    )

    if not fragments:
        print("No fragments selected (no variants or start position beyond variants). Exiting.", file=sys.stderr)
        sys.exit(0)

    # For each sample, create output fasta and write selected fragments (apply variants)
    for sample_idx, sample_name in enumerate(samples):
        if args.sample and sample_name != args.sample:
            continue

        # output filename: include start chr/pos if user provided them
        if args.start_chr and args.start_pos:
            out_name = f"{sample_name}.consensus.{args.start_chr}_{args.start_pos}.fasta"
        else:
            out_name = f"{sample_name}.consensus.fasta"
        out_path = os.path.join(args.out, out_name)

        with open(out_path, "w") as outfh:
            for chrom, start, end in fragments:
                if chrom not in genome:
                    print(f"Warning: chromosome {chrom} not found in FASTA; skipping fragment {chrom}:{start}", file=sys.stderr)
                    continue
                ref_seq = genome[chrom]
                chrom_variants = sorted_variants_by_chrom.get(chrom, [])
                # apply variants for this sample (sample_idx is index in VCF samples)
                consensus = apply_variants_to_reference(
                    ref_seq=ref_seq,
                    chrom_variants=chrom_variants,
                    sample_index=sample_idx,
                    window_start=start,
                    window_len=args.fragment_length,
                    af_threshold=args.af_threshold
                )
                header = f">{chrom}_{start}_{end}"
                outfh.write(header + "\n")
                outfh.write(wrap_fasta(consensus) + "\n")

        print(f"Wrote {out_path}", file=sys.stderr)

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()
