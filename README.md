# Whole genome sequencing of Crimean Congo hemorrhagic fever virus circulating in Pakistan during 2022 - reproducible segmented analysis

Reproducible pipeline that mirrors the study design and methods. It builds contigs, uses NCBI BLAST to identify the likely clade for each segment, then performs targeted reference based assembly per segment, and finally infers phylogeny for S, M, and L. It also flags potential reassortment if a sample segments into different clades.

## Program summary
One Snakemake pipeline. Discovery first, then per segment targeted analysis.

1) Inputs
   - Paired end FASTQ from metagenomic RNA libraries sequenced on Illumina (MiSeq 2x150 or similar).

2) QC and trimming
   - FastQC for initial QC.
   - Trimmomatic with typical metagenomics settings: sliding window 4:15, minimum length 50, leading and trailing quality 3.

3) Pathogen discovery
   - SPAdes de novo assembly of trimmed reads.
   - Contig QC table with N50, GC percent, and percent N for contigs. Length cutoff default 300 nt.

4) BLAST based clade identification
   - BLASTN contigs against nt in remote mode or a local nt mirror.
   - Classify contigs as S, M, or L using a combination of contig length ranges and hit descriptions.
   - For each segment, select a best nucleotide reference accession from top hits and record a heuristic clade label if a known representative is matched.
   - Write a per segment summary with the chosen reference and provisional clade call.

   Expected segment sizes used by the classifier
   - S: about 1670 nt
   - M: about 5360 nt
   - L: about 12150 nt

5) Targeted reference based assembly per segment
   - BWA MEM mapping to the chosen reference for S, M, and L.
   - Sort, index, call variants, mask sites with depth below a threshold (default 10), and write masked consensus per segment with bcftools.
   - Combine all sample segments in `results/consensus/` for downstream analyses.

6) Per segment phylogeny and clade confirmation
   - Build a context set automatically by fetching the top N GenBank accessions from the segment specific BLAST results.
   - MAFFT alignment of sample consensus plus context.
   - IQ TREE ML tree with 1000 ultrafast bootstraps and optional ModelFinder.
   - Write a compact report per segment that lists the nearest labeled sequences in the tree to support a clade interpretation.

7) Reassortment flag
   - If the provisional clade labels differ across S, M, and L, write a `reassortment_flag.txt` for that sample.

## Why this design
CCHFV segments are known to fall into lineages that can differ across S, M, and L, and reassortment is documented. This pipeline follows the paper design for Pakistan 2022 and general best practice for segmented bunyaviruses.

## Requirements
- Python 3.11 or newer
- Option A: conda or mamba
- System tools installed via conda: fastqc, trimmomatic, spades, blast, bwa, samtools, bcftools, mafft, iqtree, entrez direct, seqkit, snakemake, biopython, pyyaml

### NCBI usage note
Set a contact email once per shell for E utilities. Optional API key improves rate limits.
```bash
export NCBI_EMAIL="you@example.com"
export NCBI_API_KEY="xxxxxxxxxxxxxxxxxxxxxxxxxxxx"   # optional
```

## One command run
```bash
conda env create -f env/environment.yml
conda activate cchfv-seg-env
export NCBI_EMAIL="you@example.com"
snakemake -s workflow/Snakefile -c 4 --printshellcmds
```

## Configuration
Edit `config/config.yaml`. Minimal example:
```yaml
pairs:
  - sample: CCHF_2022_ISB_01
    r1: data-private/CCHF_2022_ISB_01_R1.fastq.gz
    r2: data-private/CCHF_2022_ISB_01_R2.fastq.gz

discovery:
  min_len_contig: 300
  blast_remote: true
  blast_db: nt
  top_hits: 25

phylogeny:
  bootstrap: 1000
  use_model_finder: true
  min_depth_consensus: 10

context:
  per_segment_hits: 25
```

## Outputs
- `results/spades/<sample>/contigs.fasta`
- `results/spades/<sample>/contigs.qc.tsv` and `.summary.txt`
- `results/blast/<sample>.contig_top_hits.tsv`
- `results/segments/<sample>_segment_calls.tsv` with provisional clade labels
- `results/refs/<sample>.S.fasta`, `.M.fasta`, `.L.fasta` chosen references
- `results/consensus/<sample>.S.fa`, `.M.fa`, `.L.fa` masked consensus
- `results/aln/S_alignment.fasta`, `M_alignment.fasta`, `L_alignment.fasta`
- `results/iqtree/S.treefile`, `M.treefile`, `L.treefile`
- `results/reports/<sample>_segment_report.txt`
- `results/reports/<sample>_reassortment_flag.txt` present only if mixed clades

## Quick start
1) Drop your FASTQs under `data-private/` and set them in `config/config.yaml`.
2) Ensure you set `NCBI_EMAIL` in your shell.
3) Run `make run` or the Snakemake command above.
4) Inspect `results/reports/` and the per segment trees in `results/iqtree/`.

## References
- Umair M, Rehman Z, Haider SA, et al. Whole genome sequencing of Crimean Congo hemorrhagic fever virus circulating in Pakistan during 2022. Journal of Medical Virology. 2023. doi:10.1002/jmv.28604
- Deyde VM, Khristova ML, Rollin PE, Ksiazek TG, Nichol ST. CCHFV genomics and global diversity. J Virol. 2006;80(17):8834-8842.
- Bente DA, Forrester NL, Watts DM, et al. Crimean Congo hemorrhagic fever: history, epidemiology, pathogenesis, clinical syndrome and genetic diversity. Antiviral Res. 2013;100(1):159-189.
- Guo R, Shen S, Zhang Y, et al. A new strain of CCHFV isolated from Xinjiang, China. Virologica Sinica. 2017;32(1):80-88.
- Duh D, Saksida A, Petrovec M, et al. Complete genome sequence of a CCHFV strain. Virol J. 2008;5:11.
- Tools: FastQC, Trimmomatic, SPAdes, BLAST+, BWA, SAMtools, BCFtools, MAFFT, IQ TREE, Entrez Direct, Snakemake.
