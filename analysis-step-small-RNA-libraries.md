# Small cfRNA-seq Data Analysis
## 1. Trim adaptor
```bash
## For NEB libraries:
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 16 --trim-n -q 30 \
-o >(pigz -c -p {threads} > {trimmed}) {input} > {log} 2>&1

## For smart libraries:
cutadapt -a AAAAAAAAAAA -g GGGGG -m 16 --trim-n -q 30 \
-o >(pigz -c -p {threads} > {trimmed}) {input} > {log} 2>&1

```
## 2. Sequence alignment
- Priority of sequential alignment
- spikein, univec, rRNA, lncRNA, miRNA, mRNA, piRNA, snoRNA, snRNA, srpRNA, tRNA, tucpRNA, Y RNA, genome, circRNA
```bash
## See snakefiles/sequential_mapping.snakemake for details in sequential alignment
pigz -d -c {inputFastq} \
        | bowtie2 -p 1 --norc --sensitive --no-unal \
            --un-gz {unmappedFastq} -x {sequenceIndex} - -S - \
        | samtools view -b -o {bam}
```


## 3. Quantification
```bash

## Count reads aligned to miRNA
bin/count_reads.py count_mature_mirna -i {miRNABam} -a {annotation}  -o {output}

## Count reads aligned to other transcripts
bin/count_reads.py count_transcript -i {transcriptBam} -s {strandness}  -o {output}

```

