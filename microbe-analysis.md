## Analysis of unmapped reads

The unmapped reads were independently analyzed using an alignment free and an alignment based approach

### The alignment free approach

- download kraken2 standard database (https://benlangmead.github.io/aws-indexes/k2)
- note the database should contain human sequence
- process unmapped reads with kraken2

```{bash}
kraken2 --db {input.kraken2db}  --unclassified-out {params.outprefix}  --report {output.report} --paired --use-names  {input.unmapped_1} {input.unmapped_2}
```

### The alignment based approach

- As rRNA dominates bacterial reads, first extract rRNA reads with sortMeRNA

- extract rRNA reads, follow https://github.com/sortmerna/sortmerna
```{bash}
sortmerna --paired_in -ref rRNA.fasta --idx-dir sortmerna-index --reads unmapped/${sample_id}_1.fastq.gz --reads unmapped/${sample_id}_2.fastq.gz --fastx --workdir output/microbe-rRNA-split/${sample_id} --index 0 --out2 --other --threads 20
```

- map rRNA reads to SLIVA rRNA

```{bash}
bowtie2 --no-unal -p 4 -1 ${fastq_1} -2 ${fastq_2}  --un-conc-gz output/rRNA/unmapped/${sample_id}_%.fastq.gz --no-discordant --end-to-end -x rRNA/bowtie2-index/rRNA | samtools view -b > ${output}
```

- count rRNA reads
```{bash}
# -s means strand of the RNA-seq library. Take a look at metadata.txt
bin/count-rRNA.py -m 0 -b $bam -s $lib -c output/rRNA/count/control-MAPQ0/${sample_id}.txt --stats output/rRNA/count/control-MAPQ0/${sample_id}.stat
```
