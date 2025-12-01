# vcf_consensus
*Test task* 

## consensus_builder.py
Pure-Python consensus sequence generator from multi-sample VCF + reference FASTA.

## Input test vcf: 
https://drive.google.com/file/d/1p0PLCRzhxGRlZ3fwh9KZVMKFABirDPsQ/view?usp=sharing
### Note: This vcf has 10 samples from 1kg dataset in hg19 build, and only a subset of positions of 6 chromosomes. VCF was manipulated to add consequent variats to test consensus.

## Conditions:
- No external libraries (only Python standard library).
- VCF is read and sorted by chromosome and position.
- Fragments are chosen based on VCF positions:
    * Start at user-provided chr:pos (if given) or the first VCF position.
    * Make a fragment of requested length around that start (start..start+len-1).
    * Skip forward to the next VCF position outside the current fragment and repeat.
    * Continue until requested number of fragments is produced (or "all" is chosen).
- Produces one FASTA per sample, filenames:
    {sample}.consensus.fasta
    or if start_chr/start_pos provided:
    {sample}.consensus.{chr}_{pos}.fasta
- AF read from INFO (AF=...) when present; otherwise AF computed from genotype counts.
- ALT wins over REF (in passes AF threshold) if at least one ALT allele is present in genotype.
