# Du_2023

## Condensate cooperativity underlies transgenerational gene silencing

This repository holds all scripts necessary to perform small RNA alignment and plot construction described in Du et. al. 2023.

## Summary of scripts

### Preparation and requirements

To run the pipeline script, bowtie v1.2.1.1, bedtools v2.30.0, and R v4.0.3 must be installed.
The reference directory should contain subdirectories for the WS230 reference genome, W230 reference genome with oma-1::gfp and cdk-1::gfp transgenes added, and WS230 reference genome with mex-5p::gfp added. All subdirectories should contain the WS230 genome in FASTA format, available from Wormbase (and supplemented with pseudochromosomes corresponding to transgene sequences).
A bowtie index is necessary for alignment to genome, splice junctions, structural RNAs and miRNA hairpins. To build all necessary bowtie index files, run
```
HOME=`pwd`
cd $HOME/reference/WS230/
path_to_genome_fasta=$HOME/reference/WS230/c_elegans.WS230.genomic.fa
bowtie-build ${path_to_genome_fasta} c_elegans.WS230.genomic
bowtie-build ce_WS230.coding_transcript.juncs.fa ce_WS230.coding_transcript.juncs
bowtie-build ce_WS230.rna.knownRNA.fa ce_WS230.rna.knownRNA
bowtie-build ce_hairpin.dna.fa ce_hairpin.dna

cd $HOME/reference/WS230.reporter/
path_to_genome_fasta=$HOME/reference/WS230.reporter/c_elegans.WS230.genomic.fa
bowtie-build ${path_to_genome_fasta} c_elegans.WS230.genomic
bowtie-build ce_WS230.coding_transcript.juncs.fa ce_WS230.coding_transcript.juncs
bowtie-build ce_WS230.rna.knownRNA.fa ce_WS230.rna.knownRNA
bowtie-build ce_hairpin.dna.fa ce_hairpin.dna

cd $HOME/reference/WS230.mex5.gfp/
path_to_genome_fasta=$HOME/reference/WS230.mex5.gfp/c_elegans.WS230.genomic.fa
bowtie-build ${path_to_genome_fasta} c_elegans.WS230.genomic
bowtie-build ce_WS230.coding_transcript.juncs.fa ce_WS230.coding_transcript.juncs
bowtie-build ce_WS230.rna.knownRNA.fa ce_WS230.rna.knownRNA
bowtie-build ce_hairpin.dna.fa ce_hairpin.dna
```

Compressed fastq.gz files must be located in the fastq directory.

### Alignment script

To map reads contained in the fastq directory and calculate normalized reads per million against genes, run bowtie_alignment.sh
```
sh ./bowtie_alignment.sh
```
The directory config should contain a TXT file with 3 comma separated fields fields: path to the fastq.gz file, name of the library, and the genome to use for alignment.

### Coverage script

To generate browser images and metagene plots, coverage must first be calculated using compute_coverage.sh
```
sh ./compute_coverage.sh
```
The alignment script must be run prior to this step.

### Browser image and metagene construction

Run browser_meta_plots.R to use coverage data to generate browser images against GFP transgenes (piRNA reporter and RNAi reporter), and to generate metagene plots using IP data.
This script depends on R libraries:
ggplot2
scales
reshape2
eulerr
gridExtra
dplyr

```
Rscript ./browser_meta_plots.R
```
