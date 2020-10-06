import argparse
from sklearn.utils import resample
from sklearn.model_selection import StratifiedShuffleSplit
from tqdm import tqdm
from fsutils import *

def overlap(x,y,n):
    xy = set(x).intersection(set(y))
    k = len(x)
    r = len(xy)
    return (r*n-k**2)/(k*n-k**2)

def KI(L,n):
    K = len(L)
    score = 0
    for i in range(K-1):
        for j in range(i+1,K):
            score += overlap(L[i],L[j],n)
    score /= K*(K-1)/2
    return score

def FS(X,y,selector,recurrency=1,n_features=10):
    if recurrency <= 1:
        recurrency == 1
        importance = selector(X,y).reshape(-1)
    else:
        importance = np.zeros(X.shape[1])
        sss = StratifiedShuffleSplit(n_splits=recurrency,test_size=0.2)
        for train_index, test_index in sss.split(X, y):
            importance_ = selector(X[train_index,:],y[train_index]).reshape(-1)
            idx = importance_.argsort()[::-1][:n_features]
            importance[idx] += 1
    return importance.argsort()[::-1][:n_features]


parser = argparse.ArgumentParser(description="Cross Validation")
parser.add_argument("--input","-i",help="Input matrix",required=True)
parser.add_argument("--pos","-p",help="positive class",required=True)
parser.add_argument("--neg","-n",help="negative class",required=True)
parser.add_argument("--resampling","-rs",type=int,help="Times of resampling",default=100)
parser.add_argument("--recurrency","-rc",type=int,help="Recurrency cutoff",default=1)
parser.add_argument("--method","-m",help="Method for feature selection",required=True,choices=selectors.keys())
parser.add_argument("--features","-f",help="Selected features",required=True)
args = parser.parse_args()


mat = pd.read_csv(args.input,sep="\t",index_col=0)
pos_ids = open(args.pos).read().split("\n")
neg_ids = open(args.neg).read().split("\n")
pos_ids = list(set(pos_ids).intersection(set(mat.columns)))
neg_ids = list(set(neg_ids).intersection(set(mat.columns)))
sample_ids = pos_ids + neg_ids
X = mat.fillna(0).T.loc[sample_ids,:].values
y = np.array(len(pos_ids)*[1] + len(neg_ids)*[0])


selector = selectors[args.method]
records = []
#for i in tqdm(range(args.resampling)):
splitter = StratifiedShuffleSplit(n_splits=args.resampling,test_size=0.2)
for train_index, test_index in splitter.split(X,y): 
    #X_,y_ = resample(X,y)
    X_,y_ = X[train_index,:],y[train_index]
    idx = FS(X_,y_,selector,recurrency=args.recurrency)
    records.append(idx)
    print(idx)
kuncheva_index = KI(records,X.shape[1])
f = open(args.features,"w")
print("#Kuncheva-Index\t{}".format(str(kuncheva_index)),file=f)
for i in range(args.resampling):
    print("\t".join(list(records[i].astype(str))),file=f)

f.close()
