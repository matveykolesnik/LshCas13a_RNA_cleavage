# LshCas13a\_RNA\_cleavage
Analysis of RNA-Seq data of _E. coli_ cells expressing targeting/nontargeting Type VI CRISPR-Cas system

Each directory contains data and scripts for the particular experiment:

* LshCas13a\_C3000 - RNA-Seq of total RNA extracted from _E. coli_ C3000 cells carrying activated/nonactivated LshCas13a enzyme;
* LshCas13a\_d10LVM - RNA-Seq of total RNA extracted from _E. coli_ $\Delta$ 10LVM cells carrying activated/nonactivated LshCas13a enzyme;
* LshCas13a\_in\_vitro\_total\_RNA - RNA-Seq of total RNA extracted from _E. coli_ C3000 cells after _in vivo_ incubation with activated/nonactivated LshCas13a enzyme;
* LshCas13a\_in\_vitro\_tRNAs - RNA-Seq of total tRNA sample after _in vivo_ incubation with activated/nonactivated LshCas13a enzyme;

Each directory contains the following subdirectories:

* Data - directory containing the raw reads data;
* Annotations - directory containing GFF tables with genomic features;
* Alignments - directory containing alignments produced with read_mapping.sh script;
* Reference_sequences - directory containing FASTA files of sequences used for reads mapping;
* Scripts - directory containing scripts for data processing;
* Results - directory containing the results of data processing.

The "Results" directory contains the following subdirectories:

* Tables
    * Ends\_counts - contains files with coordinates of 5' ends of fragments;
    * Fragment\_coords - contains files with coordinates of fragments (SeqID - Fragment\_start - Fragment\_end - Strand)
    * Merged\_ends\_counts - contains tables with 5' ends counts derived from samples designeted for comparison
    * Read\_pairs\_TABs - contains tables with coordinates of read pairs.
* WIG_files - contains wig-files with 5' ends coverage.

The "Scripts" directory contains a set of scripts for the data processing. There is a "basic" set of scripts which is common for all experiments:

* raw\_data\_processing.sh - performs reads quality assessment, removes adapters and discards low-quality reads.
    * Requirements:
        * fastqc, trimmomatic
* read\_mapping.sh - maps paired-end reads to the reference sequences. Since the SAM alignments file are quite large, the output data is compressed using gzip.
    * Requirements:
        * bowtie2
* return\_fragment_coords\_table.py - receives alignment files (in gzipped SAM format) and generates all tables deposed in "Result-Tables" directory (except "Merged\_ends\_counts") and produces WIG files with 5' ends coverage;
    * Requirements:
        * python3 with gzip and pandas modules
* merge\_ends\_count\_tables.py - combines 5' ends counts tables from different tables into one table
    * Requirements:
        * python3 with pandas, gzip, re and functools modules 
* TCS_calling.R - performs statistical test producing table with the position, logFC and p-value values.
    * Requirements:
        * R with dplyr, data.table, tidyr and edgeR modules 
