# cfRNA-SEEK

> Website hosted by github page: https://lulab.github.io/cfRNA-SEEK

## Scripts for cfRNA Sequencing Data Analysis

This repo contains scripts used in cfRNA-SEEK, a pipeline for identification of cancer relevant features from cfRNA sequencing (cfRNA-seq) data generated by SMART-Total protocol.  Scripts for small cfRNA-seq data processing are also available.
We also provide [snakemake files](https://github.com/lulab/cfRNA-SEEK/blob/master/snakefiles/) to tidy up the scripts.

- **SMART-total cfRNA-seq Reads Processing**

To quantify different RNA variations for human genome, as well as the abundance of different microbial taxos, see [analysis-step-smart-total-libraries.md](analysis-step-smart-total-libraries.md).

- **Feature Selection**

To identify cancer relevant features and evaluate the classification performance, see [classification.md](classification.md).

- **Small cfRNA-seq Reads Processing**

To analyze small cfRNA-seq data, see [analysis-step-small-RNA-libraries.md](analysis-step-small-RNA-libraries.md).

---

> [Dependencies]:
>
> Data analysis scripts are based on bash, python (version 3.6.8) and R (version 3.5.1). 
> Customized scripts are available in [bin/](https://github.com/lulab/cfRNA-SEEK/blob/master/bin).
>
> Adapter trimming was performed using cutadapt (version 2.3). Reads alignment was performed using STAR (version 2.5.3a_modified) and bowtie2 (version 2.3.5). Duplications in bam files are removed using MarkDuplicates in picard toolkit (version 2.20.0). bedtools (version 2.28.0) was used to assign mapped reads to different genomic regions. Count matrix was generated using featureCounts (version 1.6.2) . rMATs (version 4.0.2),  DaPar (version 0.91) and RNAeditor (version 1.0) was used for RNA alternative splicing, APA and RNA editing analysis, respectively. edgeR (version 3.24.3), RUVSeq (version 1.6.1) were used for expression data preprocessing. Differential analysis was performed using edgeR (version 3.24.3) and python package scipy.stats (version 1.4.1). Python packages scipy.stats (version 1.4.1), sklearn (version 0.22.2), skrebate (version 0.6) and imblearn (version 0.6.2) were used for feature selection, classification and performance evaluation. 
>
> The raw RNA-seq reads can be acquired using illumina's bcl2fastq software. Fragment sizes in RNA-seq libraries can be acquired using Agilent Bioanalyzer 2100 TapeStation analysis software.



> Illustration of cfRNA-SEEK work flow
![cfRNA-SEEK Pipeline](cfRNA-SEEK.png)



