import pandas as pd
import numpy as np
import argparse
parser=argparse.ArgumentParser(description='Extract Data From rMATs result')
parser.add_argument("--input","-i",type=str,required=True,help="Input rMATs result table for parsing")
parser.add_argument("--outdir","-od",type=str,required=True,help="Output dir")
parser.add_argument("--type","-t",type=str,required=True,help="Required Splicing type")
parser.add_argument("--method","-m",type=str,required=True,help="Count method, either JC or JCEC")
parser.add_argument("--pos","-p",type=str,required=True,help="Sample 1 ids file")
parser.add_argument("--neg","-n",type=str,required=True,help="Sample 2 ids file")
args = parser.parse_args()
input_mat=args.input
outdir=args.outdir
splicing_type=args.type
method=args.method
load_ids = lambda x:open(x).read().strip().split("\n")
neg_ids = load_ids(args.neg)
pos_ids = load_ids(args.pos)
df = pd.read_csv(input_mat,sep="\t")
def se2mat(se,counts=True):
    lines=se.apply(lambda x:np.array(x.strip().split(",")))#.astype(datatype))
    return np.vstack(lines.values)

idInCommon = ['GeneID', 'geneSymbol', 'chr', 'strand']
upDownExEe=['upstreamES','upstreamEE', 'downstreamES', 'downstreamEE']
MXE = idInCommon+['1stExonStart_0base','1stExonEnd', '2ndExonStart_0base', '2ndExonEnd']+upDownExEe
SE = idInCommon+['exonStart_0base', 'exonEnd']+upDownExEe
RI = idInCommon+['riExonStart_0base','riExonEnd']+upDownExEe
A3SS = idInCommon+['longExonStart_0base','longExonEnd', 'shortES', 'shortEE', 'flankingES', 'flankingEE']
A5SS = A3SS
columns_dict={"MXE":MXE,"SE":SE,"RI":RI,"A3SS":A3SS,"A5SS":A5SS}
id_columns = columns_dict[splicing_type]
sample_ids = pos_ids + neg_ids
junction_ids = splicing_type + "|" + df.loc[:,id_columns].apply(lambda x:"|".join([str(each) for each in x]),axis=1).values
 
inc_1 = se2mat(df.loc[:,'IJC_SAMPLE_1'])
inc_2 = se2mat(df.loc[:,'IJC_SAMPLE_2'])
inc = np.hstack([inc_1,inc_2])
inc_df = pd.DataFrame(index=junction_ids,columns=sample_ids,data=inc)
print("{} positive samples".format(inc_1.shape[1]))
print("{} negative samples".format(inc_2.shape[1]))
print("{} junctions in total".format(inc.shape[0]))

skip_1 = se2mat(df.loc[:,'SJC_SAMPLE_1'])
skip_2 = se2mat(df.loc[:,'SJC_SAMPLE_2'])
skip = np.hstack([skip_1,skip_2])
skip_df = pd.DataFrame(index=junction_ids,columns=sample_ids,data=skip)

inc_level1 = se2mat(df.loc[:,'IncLevel1'],counts=False)
inc_level2 = se2mat(df.loc[:,'IncLevel2'],counts=False)
inc_level = np.hstack([inc_level1,inc_level2])
inc_level_df = pd.DataFrame(index=junction_ids,columns=sample_ids,data=inc_level)

stats_columns = ['IncFormLen','SkipFormLen','PValue', 'FDR','IncLevelDifference']
stats_df = df.loc[:,stats_columns]
stats_df.index = junction_ids

outpath = outdir + "/{}_{}_{}.txt"
inc_df.to_csv(outpath.format(splicing_type,method,"inc"),sep="\t")
skip_df.to_csv(outpath.format(splicing_type,method,"skip"),sep="\t")
inc_level_df.to_csv(outpath.format(splicing_type,method,"inc_level"),sep="\t")
stats_df.to_csv(outpath.format(splicing_type,method,"stats"),sep="\t")
