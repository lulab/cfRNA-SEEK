#!/usr/bin/env python
import argparse
import numpy as np
from collections import OrderedDict, defaultdict
import HTSeq
from tqdm import tqdm

def count_rRNA(args):
    bam = HTSeq.BAM_Reader(args.bam)
    counts = OrderedDict()
    min_mapping_quality = args.min_mapping_quality
    strandness = {'no': 0, 'forward': 1, 'reverse': 2}.get(args.strandness, 0)
    stats = defaultdict(int)
    for read1,read2 in tqdm(HTSeq.pair_SAM_alignments(bam)):
        stats["total"] += 1
        # ignore singletons
        if (read1 is None) or (read2 is None):
            stats['singleton'] += 1
            continue
        # ignore unmapped reads
        if not (read1.aligned and read2.aligned):
            stats['unmapped'] += 1
            continue
        if read1.aQual < args.min_mapping_quality and read2.aQual < args.min_mapping_quality:
            stats['low-MAP'] += 1
            continue
        if read1.iv.chrom != read2.iv.chrom:
            stats['diff_chrom'] += 1
            continue
        if (strandness == 1) and (not ((read1.iv.strand == '+') and (read2.iv.strand == '-'))):
            stats['improper_strand'] += 1
            continue
        if (strandness == 2) and (not ((read1.iv.strand == '-') and (read2.iv.strand == '+'))):
            stats['improper_strand'] += 1
            continue
        if read1.iv.chrom not in counts:
            counts[read1.iv.chrom] = 0
        counts[read1.iv.chrom] += 1 
        stats["counted"] += 1

    with open(args.counts,"w") as f:
        f.write("seq_id\tcount\n")
        for seq_id, count in counts.items():
            f.write(f"{seq_id}\t{count}\n")

    if args.stats is not None:
        with open(args.stats,"w") as f:
            for stat in ["total","counted","unmapped","diff_chrom","improper_strand","low-MAP"]:
                count = stats[stat]
                f.write(f"{stat}\t{count}\n")
        

def main():
    parser = argparse.ArgumentParser(description='Count paired end reads mapped to rRNA database')
    parser.add_argument('--bam','-b',type=str,required=True,help="Input bam file, should in rRNA coordiante")
    parser.add_argument('--strandness','-s',type=str,default="no",help="Strandness",choices=["forward","reverse","no"])
    parser.add_argument('--counts','-c',type=str,help="Output count matrix",required=True)
    parser.add_argument('--stats',type=str,help="Counting statistics")
    parser.add_argument('--min-mapping-quality','-m',default=0,type=int,help="Min mapping quality")
    args = parser.parse_args()
    count_rRNA(args)    


if __name__ == "__main__":
    main()


