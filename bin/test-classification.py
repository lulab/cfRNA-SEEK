from fsutils import *
import argparse
from tqdm import tqdm
import numpy as np
from sklearn.model_selection import cross_val_predict,LeaveOneOut,StratifiedShuffleSplit,GridSearchCV
from sklearn.ensemble import RandomForestClassifier
from imblearn.ensemble import BalancedBaggingClassifier
from sklearn.utils.validation import check_X_y,check_is_fitted
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import roc_auc_score,confusion_matrix,auc
from imblearn.ensemble import BalancedRandomForestClassifier
import pandas as pd
from scipy.stats import ranksums 
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.neighbors import KNeighborsClassifier

def FS(X,y,selector,recurrency=1,n_features=10):
    if recurrency <= 1:
        recurrency == 1
        importance = selector(X,y).reshape(-1)
    else:
        importance = np.zeros(X.shape[1])
        sss = StratifiedShuffleSplit(n_splits=recurrency,test_size=0.3)
        for train_index, test_index in sss.split(X, y):
            importance_ = selector(X[train_index,:],y[train_index]).reshape(-1)
            idx = importance_.argsort()[::-1][:n_features]
            importance[idx] += 1
    return importance.argsort()[::-1][:n_features]

def autoTune(estimator,X,y,params_grid):
    cv = GridSearchCV(estimator, params_grid,scoring='roc_auc',refit=True)
    cv.fit(X,y)
    print(cv.best_params_)
    return cv.best_estimator_

classifiers = {
    "DT":DecisionTreeClassifier(),
    "RF":RandomForestClassifier(n_estimators=400),
    #"RF-balanced":BalancedRandomForestClassifier(),
    "RF-balanced":BalancedRandomForestClassifier(n_estimators=400),
    "LR":LogisticRegressionCV(),
    "SVM":SVC(probability=True),
    "KNN":KNeighborsClassifier()
}

param_grids = {
    "RF-balanced":{"max_depth":[2,4,8,16,None]},
    "LR":{},
    "KNN":{"n_neighbors":[2,4,8]},
    #"SVM":{'gamma':np.logspace(-8, 3, 12)},
    "SVM":{},
    "DT":{"max_depth":[2,4,8,16,None]}
}

parser = argparse.ArgumentParser(description="Cross Validation")
parser.add_argument("--input","-i",help="Input matrix",required=True)
parser.add_argument("--labels","-l",help="sample labels",required=True)
parser.add_argument("--pos","-p",help="positive class",required=True)
parser.add_argument("--neg","-n",help="negative class",required=True)
parser.add_argument("--recurrency","-rc",type=int,help="Times of recurrency",default=1)
parser.add_argument("--selector","-s",help="Method for feature selection",required=True,choices=selectors.keys())
parser.add_argument("--classifier","-c",help="Classifier",required=True,choices=classifiers.keys())
parser.add_argument("--auroc","-roc",help="Performance output",required=True)
args = parser.parse_args()


labelPath = args.labels
matPath = args.input
mat = pd.read_csv(matPath,sep="\t",index_col=0)
labels = pd.read_csv(labelPath,sep="\t",index_col=0)
label_field = "labels"
pos_classes = args.pos.strip().split(",")
neg_classes = args.neg.strip().split(",")
sample_ids = set(labels.index).intersection(mat.columns)
labels = labels.loc[sample_ids,:]
pos = labels.apply(lambda x:x.loc[label_field] in pos_classes,axis=1)
neg = labels.apply(lambda x:x.loc[label_field] in neg_classes,axis=1)
discovery = labels["dataset"]=="discovery"
validation = labels["dataset"]=="validation"
pos_ids_train = list(labels[(pos&discovery)].index)
pos_ids_test = list(labels[(pos&validation)].index)
neg_ids_train = list(labels[(neg&discovery)].index)
neg_ids_test = list(labels[(neg&validation)].index)
sample_ids = pos_ids_train + pos_ids_test + neg_ids_train + neg_ids_test
X = mat.T.loc[sample_ids,:].values

#if args.scaling:
X = StandardScaler().fit_transform(X)


X = pd.DataFrame(X,index=sample_ids,columns=mat.index)
X_train = X.loc[pos_ids_train+neg_ids_train,:].values
y_train = np.array([1]*len(pos_ids_train)+[0]*len(neg_ids_train))
train_index = np.arange(0,X_train.shape[0])
np.random.shuffle(train_index)
X_train = X_train[train_index,:]
y_train = y_train[train_index]
X_test = X.loc[pos_ids_test+neg_ids_test,:].values
y_test = np.array([1]*len(pos_ids_test)+[0]*len(neg_ids_test))
test_index = np.arange(0,X_test.shape[0])
np.random.shuffle(test_index)
X_test = X_test[test_index,:]
y_test = y_test[test_index]

selector = selectors[args.selector]
idx = FS(X_train,y_train,selector,recurrency=args.recurrency)
clf = autoTune(classifiers[args.classifier],X_train[:,idx],y_train,param_grids[args.classifier])
#clf = classifiers[args.classifier]
#clf.fit(X_train[:,idx],y_train)
y_pred = clf.predict_proba(X_test[:,idx])[:,1]

f = open(args.auroc,"w")
print("{}-{}-auroc:{}".format(args.selector,args.classifier,roc_auc_score(y_test,y_pred)),file=f)
f.close()



