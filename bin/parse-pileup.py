import argparse
parser = argparse.ArgumentParser(description="Extract Editing Info From pileup results")
parser.add_argument("--input","-i",help="input mpileup",required=True)
parser.add_argument("--output","-o",help="output",required=True)
args = parser.parse_args()

editing = {"A":"G","T":"C"}

fout = open(args.output,"w")
print("chrom","pos","ref","n_ref","n_edit","depth",sep="\t",file=fout)
with open(args.input) as fin:
    for line in fin:
        chrom,pos,ref,depth,match,qual = line.strip().split("\t")
        ref = ref.upper()
        match = match.upper()
        n_ref = match.count(".")
        n_edit = match.count(editing[ref])
        print(chrom,pos,ref,n_ref,n_edit,depth,sep="\t",file=fout)

fout.close()

