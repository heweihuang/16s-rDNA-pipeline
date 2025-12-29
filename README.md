# 16s-rDNA-pipeline
Perform automated analysis of microbial (bacterial) 16s rDNA data (V3-V4), including: 1. Quality control: removal of chimeras, short sequences, etc. 2. OTU annotation after splicing. 3. Species abundance analysis, such as alpha and beta diversity. 4. Phylogenetic tree plotting. 5. Analysis of species differences between groups, etc.


Usage:

1. Install QIIME 1 and Nextflow framework on the server.

2. Download the pipeline code.

3. Run the command:
srun nextflow run -with-trace -with-timeline -with-report -with-dag -resume -offline -dsl1 16S.nf
