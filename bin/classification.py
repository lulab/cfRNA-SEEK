from BalancedBinaryClassifier import *
from fsutils import *
from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV,StratifiedShuffleSplit,cross_validate,StratifiedKFold
import argparse
import pickle
from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from sklearn.inspection import permutation_importance
from sklearn.utils import resample
from copy import copy

parser = argparse.ArgumentParser(description="Cross Validation")
parser.add_argument("--input","-i",help="Input matrix",required=True)
parser.add_argument("--scaling","-s",help="sacling",action="store_true",default=False)
parser.add_argument("--labels","-l",help="label of the samples",required=True)
parser.add_argument("--pos","-p",help="positive class",required=True)
parser.add_argument("--neg","-n",help="negative class",required=True)
parser.add_argument("--label_field","-lf",help="which field in the covariate table to use as label",default="labels")
parser.add_argument("--model","-m",help="where to save model",required=True)
parser.add_argument("--trend","-t",help="trend of selected markers",choices=["up","both","both_balance"],default="both")
parser.add_argument("--restrain",help="restrain the searching space to a gene set",default=None)
parser.add_argument("--performance",help="performance",required=True)
parser.add_argument("--features","-f",help="information of features finaly used",required=True)
parser.add_argument("--roc",help="path of testing roc curve",required=True)
parser.add_argument("--selector",help="which selector to use",default="ranksum-SURF",choices=["ranksum","relief","ranksum-SURF"])
parser.add_argument("--cv_features","-cf",help="cv features",default=None)
parser.add_argument("--n_jobs",type=int,help="number of workers for cross validation",default=1)
args = parser.parse_args()


fout = open(args.performance,"w")
i = 0

if args.cv_features is not None:
    out_cv_features = open(args.cv_features,"w")

def autoTune(estimator,X,y,params_grid):
    cv = GridSearchCV(estimator, params_grid,scoring='roc_auc',refit=True,cv=StratifiedKFold(shuffle=True,random_state=777))
    cv.fit(X,y)
    print(cv.best_params_)
    return cv.best_estimator_

def _my_scorer(clf, X_val, y_true_val):  
    global fout
    global i
    # do all the work and return some of the metrics
    y_pred_val = clf.predict_proba(X_val)[:,1]
    accuracy = metrics.roc_auc_score(y_true_val,y_pred_val)
    #print("\t".join(params)+"\troc_auc:{}".format(accuracy))
    i += 1
    print("cross validation: {}".format(str(i).zfill(3)))
    print("{}\tcross_validation\t{}".format(str(i).zfill(3),accuracy),file=fout,sep="\t")
    return accuracy

print("Load data ...")
labelPath = args.labels
matPath = args.input
mat = pd.read_csv(matPath,sep="\t",index_col=0).fillna(0)
impute = np.repeat(mat.median(axis=1).values.reshape((-1,1)),mat.shape[1],axis=1)
impute = pd.DataFrame(index=mat.index,columns=mat.columns,data=impute)
mat = mat.fillna(impute)

labels = pd.read_csv(labelPath,sep="\t",index_col=0)
label_field = args.label_field
pos_classes = args.pos.strip().split(",")
neg_classes = args.neg.strip().split(",")

sample_ids = set(labels.index).intersection(mat.columns)
labels = labels.loc[sample_ids,:]

pos = labels.apply(lambda x:x.loc[label_field] in pos_classes,axis=1)
neg = labels.apply(lambda x:x.loc[label_field] in neg_classes,axis=1)
discovery = labels["dataset"]=="discovery"
validation = labels["dataset"]=="validation"


pos_ids_train = sorted(list(labels[(pos&discovery)].index))
pos_ids_test = sorted(list(labels[(pos&validation)].index))
neg_ids_train = sorted(list(labels[(neg&discovery)].index))
neg_ids_test = sorted(list(labels[(neg&validation)].index))

cvOnly = False
if len(pos_ids_test) == 0:
    cvOnly = True
    print("The positive class is not present in test set, only perform cross validation")

sample_ids = pos_ids_train + pos_ids_test + neg_ids_train + neg_ids_test
X = mat.T.loc[sample_ids,:].values
print("Done .")

if args.scaling:
    X = StandardScaler().fit_transform(X)

space = set()
mask=None
if args.restrain is not None:
    with open(args.restrain) as f:
        for line in f:
            space.add(line.strip().split("|")[0].split(".")[0])
    mask = np.array(mat.index.map(lambda x:x.split("|")[0].split(".")[0] not in space))
    print("{} features in {} were masked".format(sum(mask),len(mask)))


X = pd.DataFrame(X,index=sample_ids,columns=mat.index)

print(pos_ids_test[:5])
print(neg_ids_test[:5])
X_train = X.loc[pos_ids_train+neg_ids_train,:].values
y_train = np.array([1]*len(pos_ids_train)+[0]*len(neg_ids_train))

if not cvOnly:
    X_test = X.loc[pos_ids_test+neg_ids_test,:].values
    y_test = np.array([1]*len(pos_ids_test)+[0]*len(neg_ids_test))

selector = selectors[args.selector]


clf = BalancedBinaryClassifier(trend=args.trend,space_mask=mask,selector=selector)
#splitter = StratifiedShuffleSplit(n_splits=100,test_size=0.2,random_state=666)
splitter = StratifiedShuffleSplit(n_splits=100,test_size=0.2,random_state=777)
print("Cross validation ...")
cv = cross_validate(clf,X_train,y_train,return_estimator=True,scoring=_my_scorer,cv=splitter,verbose=0,n_jobs=args.n_jobs)
print("Done .")
feature_counts = np.zeros(X_train.shape[1])
for estimator in cv['estimator']:
    features = estimator.features
    feature_counts[features] += 1
    features_ = [str(x) for x in features]
    if args.cv_features is not None:
        out_cv_features.write("\t".join(features_)+"\n")
if args.cv_features is not None:
    out_cv_features.close()
if args.trend == "both_balance":
    posMean = X.loc[pos_ids_train+pos_ids_test,:].mean(axis=0).values
    negMean = X.loc[neg_ids_train+neg_ids_test,:].mean(axis=0).values
    upFreq = feature_counts.copy()
    upFreq[np.where((posMean-negMean)<0)[0]] = 0 
    downFreq = feature_counts.copy()
    downFreq[np.where((posMean-negMean)>0)[0]] = 0
    features = np.array(list(upFreq.argsort()[::-1][:5]) + list(downFreq.argsort()[::-1][:5]))
    
else:
    features = feature_counts.argsort()[::-1][:10]
if not cvOnly:
    clf2 = BalancedRandomForestClassifier(n_estimators=100,random_state=666)
    grid = {"max_depth":[2,4,8,16,None]}
    print(X_test[:,features][:4,:])
    clf2 = autoTune(clf2,X_train[:,features],y_train,grid) 
    y_pred = clf2.predict_proba(X_test[:,features])[:,1]
    print("ROC-AUC on test set")
    print("full-training-set\ttest_set\t{}".format(metrics.roc_auc_score(y_test,y_pred)),sep="\t",file=fout)
    fpr, tpr, thresholds = metrics.roc_curve(y_test,y_pred)
    roc_data = {}
    roc_data["fpr"] = fpr
    roc_data["tpr"] = tpr
    roc_data["thresholds"] = thresholds
    pd.DataFrame(roc_data).to_csv(args.roc,index=False,sep="\t")
fout.close()
feature_info = pd.DataFrame(index=mat.index,columns=["frequency","importance"])
feature_info["frequency"] = feature_counts
if not cvOnly:
    feature_info.loc[mat.index[features],"importance"] = clf2.feature_importances_
    with open(args.model,"wb") as f:
        pickle.dump(cv,f)
feature_info.to_csv(args.features,sep="\t")

