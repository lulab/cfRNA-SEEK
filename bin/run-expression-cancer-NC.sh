#!/bin/bash
labels=metadata/covariates.txt
data=expression
input=variations/${data}.txt
for trend in up both both_balance;do
for cancer in CRC STAD LUAD HCC ESCA;do
outdir=output-shuffle/${data}/${cancer}-${trend}-NC
[ -d ${outdir} ] || mkdir -p ${outdir}
bsub -q Z-LU "python  -W ignore bin/classification.py --input ${input} --labels ${labels} --pos ${cancer} --neg NC --model  ${outdir}/model.pickle --performance ${outdir}/performance.txt --roc ${outdir}/roc.txt --features ${outdir}/features.txt --trend ${trend} > log/${trend}-${data}-${cancer}-NC.txt 2>&1 "
done
outdir=output-shuffle/${data}/pan-cancer-${trend}-NC
[ -d ${outdir} ] || mkdir -p ${outdir}
bsub -q Z-LU "python  -W ignore bin/classification.py --input ${input} --labels ${labels} --pos CRC,STAD,LUAD,HCC --neg NC --model  ${outdir}/model.pickle --performance ${outdir}/performance.txt --roc ${outdir}/roc.txt --features ${outdir}/features.txt --trend ${trend} > log/${trend}-${data}-pan-cancer-NC.txt 2>&1 "
done
