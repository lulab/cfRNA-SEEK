# SMART-total cfRNA-seq Data Analysis
## 1.1 Trim adaptor 
```bash
cutadapt --pair-filter any  -q 30,30 \
            -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
            --trim-n -m 16  -o >(gzip -c > {output.fastq1}) -p >(gzip -c > {output.fastq2}) \
            {input.fastq1} {input.fastq2} > {log} 2>&1
```
## 1.2 Trim GC oligo introduced in template switching
- For reverse stranded (first strand) libraries, set strandness to reverse, for forward stranded (second strand) libraries, set strandness to forward
- Read pairs shorter than 30 nt are discarded
- Input should be compressed by gzip, named as `{sample_id}_{1,2}.fastq.gz`
```bash
python bin/trimGC.py -s {strandness} -o {output_preffix} -i {input_preffix} > {log} 2>&1
```

## 2.1 Reads alignment
- Sequentially align reads to spikeIn,UniVec,rRNA,hg38,and circRNA using the following command
  - For spikeIn,UniVec,rRNA,and circRNA alignment, set seedPerWindowNmax to 20
  - For hg38 alignment, set seedPerWindowNmax to 50
  - Note unmapped reads of STAR 2.5.4a (also called STAR 2.5.3a\_modified) may out of order in multi-thread alignment (https://github.com/alexdobin/STAR/issues/222), the latest version of STAR is recommended. Alternatively, unmapped fastq files should be repaired (using repair.sh in bbmap suite for example https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/repair-guide/) before downstream analysis, or the genome alignment rate could be extremely low.
```bash
STAR --genomeDir {sequenceIndex} \
            --readFilesIn {input.reads1} {input.reads2} \
            --runThreadN 4 \
            --outFileNamePrefix {output_prefix} \
            --outSAMtype BAM Unsorted \
            --outReadsUnmapped Fastx \
            --readFilesCommand gzip -d -c \
            --outSAMmultNmax 1 \
            --seedPerWindowNmax {seedPerWindowNmax}
```

## 2.2 Remove duplication in circRNA.bam and genome.bam with picard tools
```bash
java -jar {picardDir}/picard.jar \
            MarkDuplicates REMOVE_DUPLICATES=true \
            ASSUME_SORT_ORDER=queryname \
            I={input.bam} \
            O={output.bam} \
            M={output.metrics} \
            READ_NAME_REGEX=null
```

## 2.3 Get intron spanning reads from genome aligned reads
- If any mate in a read pair contains N in its CIGAR string, this fragment was considered as intron-spanning
```bash
bash bin/getIntron-spanning.sh {inbam} {outbam} > {log} 2>&1
```

## 3.1 Assign reads in genome bam file to certain genomic regions
- The priority of the assignment is defined in `config/priority.txt`
- bed files corresponding to regions in `config/priority.txt` should be present in {beddir}, named as {region}.bed, in bed6 format
```bash
# For forward stranded libraries, '-S' in bedtools coverage -counts -S -sorted -a - -b ${beddir}/${region}.bed should be replaced by '-s'
bash bin/sequential.assign.long.sh {bam} {outdir} {beddir}
```

## 3.2 Quantify gene and circRNA expression
```bash
# Quantify gene expression
# For forward stranded libraries,strandness=1, for reverse stranded libraries, strandness=2 
featureCounts -O -t exon -g gene_id -M -s {strandness} -p -a {gtf} -o {counts} {bam} > {log}
# Quantify circRNA expression
# Strandness: forward, reverse
bin/count_reads.py count_circrna -s {strandness} --paired-end -i {inbam} -o {output}
```

## 4.1 Calculate PSI scores
```bash
## Running rMATs (rMATS.4.0.2/rMATS-turbo-Linux-UCS4)
python2 {rmats_path} --b1 {pos_bam_paths} --b2 {neg_bam_paths} --gtf {gtfFile} --od {outdir} -t paired --readLength 150
## Summarize rMATs results for {splicing_type} (one of MXE A3SS A5SS SE RI)
python3 bin/summarize-splicing.py --input {splicing_type}.MATS.JC.txt  --outdir outdir --type {splicing_type} --method JC  --pos pos_ids.txt  --neg  neg_ids.txt
```

## 4.2 Calculate PDUI scores
```bash
## Sort bam file by coordinate
samtools sort  -o {bamUnsorted} {bamSorted}
samtools index  {bamSorted}

## Calculate coverage, strand in +,-
bedtools genomecov -strand {strand}  -split -ibam {bamSorted} \
            | LC_COLLATE=C sort -T {temp_dir} -k1,1 -k2,2n > {bedgraph}
bedGraphToBigWig {bedgraph} {chrom_sizes} {bigwig}

## Merge Coverage of + strand and - strand
bigWigMerge {bigwigPosStrand} {bigwigNegStrand} {bedgraphMerged}
sort -k1,1 -k2,2n {bedgraphMerged} > {bedgraphMergedSorted}
bedGraphToBigWig {bedgraphMergedSorted}  {chrom_sizes} {bigwigMerged}
bigWigToWig {bigwigMerged} {wigMerged}

## Run DaPar
{DaParsDir}/src/DaPars_main.py {config}

## Summarize results
python  bin/parseDapar.py -i {input} -c {config} -l {long} -s {short} -p {PDUI}

```

## 4.3 Calcualte RNA editing level by gene
```bash

## Identify recurrent RNA editing sites
## For each sample, run
{RNAEditorDir}/RNAEditor.py  -i  {fastq_1} {fastq_2} -c {config}

# Calculate detection recurrency for each editing sites:  bin/recurrent-editing.py

# Filter intergenic editing sites: bin/get-intragene-editing-sites.sh

## Calculate coverage for each sample at recurrently edited sites, include samples which no editing events were reported by RNAEditor
samtools mpileup -l  editing-sites-pos-list.txt -o ${output} --reference ${ref} ${bam}

# Parse pileup result (get reads number support editing/support not editing)
python bin/cal-coverage.py -i {pileup} -o {coverage}

# Summarize the coverage by genomic positions: bin/summarize-editing-coverage.py

# Calculate editing level for each gene and each editng site: bin/editing-level.py 
```

## 4.4 Classification of unmapped reads
```bash
## Classify unmapped reads with kraken2
kraken2 --db {kraken2db}  --unclassified-out {outprefix}  --report {report} --paired --use-names {unmapped_1} {unmapped_2}

## Summarize kraken2 result
## Genus level data was used for classification 
python bin/summarize-kraken.py -i {report} -l {taxoLevel} -o {output}

```
## 5.1 Filtering
- See `bin/filter.py`

- Gene expression data: retained genes with CPM > 2 for more than 50% of samples in at least one class

  ```bash
  bin/filter.py -i {input} -o {output} --stratify metadata/sample_classes.txt --proportion 0.5 --stratify_key label --threshold 2 --pass_gene_ids  {detected_gene} --method by_value
  ```

- PSI scores and PDUI scores: in each class, the values of more than 80% of samples should not be null

  ```bash
  python bin/filter.py --input {input} -o {output} --stratify metadata/sample_classes.txt --stratify_key label --collapse intersection --pass_gene_ids passed-ids.txt --proportion 0.8 --method by_na
  ```

- RNA editing level: genes with more than 50% of samples with fewer than 6
  reads supporting the editing level calculation were filtered out.
  
  ```bash
  bin/filter.py --input {coverage-by-gene} -o {coverage-by-gene-filtered} --stratify metadata/sample_classes.txt --stratify_key label --collapse intersection --pass_gene_ids {gene-passed-filter} --threshold 6 --proportion 0.5 --method by_value
  ```

## 5.2 Differential analysis:
- For counts data ( expression counts and kraken2 counts ), see default method (edger-glmlrt) in `bin/differential_expression.R`

  - --positive-ids/--negative-ids: path of file contain sample ids, one id per line

  ```bash
  Rscript bin/differential_expression.R \
  -i {input}   \
  --method edger_glmlrt \
  --positive-ids {pos-ids} \
  --negative-ids {neg-ids} \
  -o {output}
  ```

- For analysis of ratio based data (psi, PDUI and editing level), see `bin/ranksum.py`

  ```bash
  python bin/ranksum.py --input {matrix} --output {diff-table} --pos_ids {pos-ids} --neg_ids {neg-ids} 
  ```

  

