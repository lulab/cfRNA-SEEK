import numpy as np
import pandas as pd
from scipy.stats import ranksums 
from skrebate import SURF
from sklearn.linear_model import LogisticRegressionCV
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import mutual_info_classif
from sklearn.preprocessing import StandardScaler

def ranksum(X,y):
    pvalue = []
    X0,X1 = X[y==0],X[y==1]
    for i in range(X.shape[1]):
        _,p = ranksums(X0[:,i],X1[:,i])
        pvalue.append(p)
    return -np.log(np.array(pvalue))

def RandomForestFS(X,y):
    return RandomForestClassifier(max_features="auto").fit(X,y).feature_importances_


def rankSURF(X,y,topN=100):
    if X.shape[1] < topN:
        return relief(X,y)
    else:
        importances_ = ranksum(X,y)
        top = importances_.argsort()[::-1][:topN]
        importances = np.zeros(X.shape[1])
        topImportance =  relief(X[:,top],y)
        for i,importance in zip(top,topImportance):
            importances[i] = importance
        return importances

def LogisticRegressionL1FS(X,y):
    X = StandardScaler().fit_transform(X)
    return np.abs(LogisticRegressionCV(penalty="l1",solver='liblinear').fit(X,y).coef_)

def relief(X,y):
    return SURF().fit(X,y).feature_importances_

def randomFS(X,y):
    return np.random.rand(X.shape[1])

selectors = {
    "random":randomFS,
    "ranksum":ranksum,
    "RF":RandomForestFS,
    "LR-L1":LogisticRegressionL1FS,
    "SURF":relief,
    "ranksum-SURF":rankSURF,
    "MI":mutual_info_classif
}


