# Scripts for Classification
## Normalization, batch correction and data clean
```bash
## counts: input count matrix
## classes: metadata
## tmm: output tmm normalized data
## ruv: output ruvg normalized data
## anova: output anova statistics table
Rscript bin/normalization.R --counts ${input} --classes metadata/sample_classes.txt --tmm ${tmm} --ruv ${ruv} --anova ${anova}
```

## Evaluation of feature selection methods and classifiers
- Evaluate the stability of feature selection using KI index
```bash
## Method: ranksum RF SURF LR-L1  MI random ranksum-SURF
python bin/FS.py --input ${input} --pos metadata/CRC-dis.txt  --neg metadata/NC-dis.txt --resampling 100 --recurrency 1 --method ${method} --features stability/${method}.txt
```
- Evaluate the performance of feature selection
```bash
## Selector: ranksum RF SURF LR-L1  MI random ranksum-SURF 
## Classifier: RF-balanced LR KNN DT SVM
python bin/test-classification.py  --input ${input} --pos CRC --neg NC --selector ${selector} --classifier ${clf} --auroc performance/${selector}:${clf}.txt --labels ${labels} 
```

## Binary feature selection and performance evaluation

- `bin/classification.py` was used for identifying cancer relevant features. Feature selection was integrated in a cross validation procedure on discovery set. By default,  this each cross validation run, distinguishable features in two groups were identified by ranksum test filtering, followed by SURF feature selection. The most frequently selected features were finally reported. Additional constraints, such as the trend of alterations, or the initial searching space, can be specified.

- For cancer vs. HD classification, see `bin/run-variation-cancer-NC.sh` and `bin/run-expression-cancer-NC.sh`

- For one vs. rest classification for different cancers, see `bin/run-OvR5.sh`

- For performance evaluation of combining selected features (mixing features or logistic regression), see `bin/get-probability.sh` and `bin/probability-integration.py`

  ```bash
  ## For example, pan cancer binary classification call the script
  python bin/classification-selected.py --input {data-matrix} --features {binary-features} --pos CRC,STAD,LUAD,HCC --neg NC  --labels metadata/covariates.txt -proba pan-cancer-prob.txt --roc {roc} --auroc {auroc}
  ```

## Multi-class classification
- See `bin/multi-class.py` for performance evaluation of multi-class classification between different cancer types.  Output LOO-CV confusion matrix on discovery set and confusion matrix  on validation set

  ```bash
  python bin/multi-class.py  --input {data-matrix} --features {features} --performance {performance} --confusion_loo {confusion-loo} --confusion_test {confusion-test}  --labels {labels} 
  ```

  
