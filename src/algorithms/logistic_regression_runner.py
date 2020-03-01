import tqdm
from tqdm import tqdm, trange
#import scikit
from rpy2.robjects import *
import subprocess
from rpy2 import robjects as ro
import numpy as np
from scipy import sparse
import src.algorithms.fastsinksource_runner as fastsinksource
import sklearn
from sklearn import metrics
from sklearn.svm import LinearSVC
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from scipy import sparse
import scipy
import seaborn as sns
import src.algorithms.alg_utils as alg_utils
import matplotlib.pyplot as plt
from sklearn.calibration import CalibratedClassifierCV
import time
from sklearn.linear_model import LogisticRegression
#from scikit.learn.svm.sparse import SVC

def setupInputs(run_obj):
    run_obj.ann_matrix = run_obj.ann_obj.ann_matrix
    run_obj.goids = run_obj.ann_obj.goids
    run_obj.prots = run_obj.ann_obj.prots
    run_obj.hpoidx = run_obj.ann_obj.goid2idx
    run_obj.protidx = run_obj.ann_obj.node2idx
    #run_obj.net_mat = run_obj.net_obj.W

    # if swsn weighting is to be used, then obtain the current fold W
    if run_obj.net_obj.weight_swsn:
        W, process_time = run_obj.net_obj.weight_SWSN(run_obj.ann_matrix)
        run_obj.params_results['%s_weight_time'%(run_obj.name)] += process_time
    else:
        W = run_obj.net_obj.W
    
    # if influence matrix is to be used, then obtain the influence matrix of W
    if run_obj.net_obj.influence_mat:
        W = alg_utils.influenceMatrix(W, ss_lambda=run_obj.params.get('lambda', None))
    
    # finally, normalize the edge weights of the W matrix.
    run_obj.P = alg_utils.normalizeGraphEdgeWeights(W, ss_lambda=run_obj.params.get('lambda', None))
 
    #run_obj.P = alg_utils.normalizeGraphEdgeWeights(run_obj.net_obj.W, ss_lambda=run_obj.params.get('lambda', None))

    # unnormalized network
    #run_obj.net_mat = run_obj.net_obj.W

    # normalized network

    run_obj.net_mat = run_obj.P

    return

def setup_params_str(weight_str, params, name):
    return "-logReg-2000ter"

def setupOutputs(run_obj):
    return

def run(run_obj):
    
    
    print("Running Logistic Regression for 2000 iterations")

    params_results = run_obj.params_results
    P, alg, params = run_obj.P, run_obj.name, run_obj.params

    # get the labels matrix and transpose it to have label names as columns
    ann_mat = run_obj.ann_matrix
    labels = ann_mat.transpose()    # genes x hpo
    #print(labels.shape)

    # get the feature vectors of the samples and check the size
    feats = run_obj.net_mat
    #print(feats.shape)
    
    
    train_mat = run_obj.train_mat
    test_mat = run_obj.test_mat

 

    # first just try multiclass classification with the first label
    scores = sparse.lil_matrix(labels.transpose().shape, dtype=np.float)        #   dim: hpo x genes
    
    combined_scores = sparse.lil_matrix(labels.transpose().shape, dtype=np.float) # dim: hpo x genes terms
    for l in tqdm(range(labels.shape[1])):
        #print("******************working on label........: {}".format(l))
        
        # compute the test gene indices of the annotations for the given label
        
        train_pos, train_neg = alg_utils.get_goid_pos_neg(train_mat,l)
        if len(train_pos)==0:
            print("Skipping term, 0 positive examples")
            continue
        #train_set = set(train_pos) | set(train_neg)
        #train_set = sorted(list(train_set))
        
        test_pos, test_neg = alg_utils.get_goid_pos_neg(test_mat, l)
        test_set = set(test_pos) | set(test_neg)
        test_set = sorted(list(test_set))
        
        #print(run_obj.protidx.values())
        #train_set = sorted(list(set(run_obj.protidx.values()) - set(test_set)))
         

        train_set = sorted(list(set(train_pos)|set(train_neg)))
        

        X_train = feats[train_set, :]
        X_test = feats[test_set]
        #print(X_test.shape)
        # construct the training label data without these indices

        #y_train = np.delete(train_mat.toarray().transpose(), test_set, 0)
        y_train = train_mat.transpose()[train_set, :]
        y_train = sparse.lil_matrix(y_train)
        temp = []
        
        # this (SVC) uses internal 5-FCV when probability option is set to true
        #clf = SVC(probability=True, verbose=True, gamma='auto', kernel='linear')
    
        clf =  LogisticRegression(max_iter=2000)

        # get the column of training data for the given label 
        lab = y_train[:,l].toarray().flatten()
        #print(lab.shape)
        # now train the model on the constructed training data and the column of labels
        
        #print("Training")

        process_time = time.process_time()
        clf.fit(X_train, lab)
        #print(clf.classes_)
        
        #process_time = time.process_time()
        # make predictions on the constructed training set
        #print("Predicting")
        predict = clf.predict_proba(X_test)[:,1]

        process_time = time.process_time() - process_time
        predict = predict.tolist()
        
        assert(len(predict) == X_test.shape[0]), "Predicted labels set not the same as test set length"
        assert any(v > 0 for v in predict), "All values in predicted label are 0"

        # get the current scores for the given label l
        curr_score = scores[l].toarray().flatten()
        # for the test indices of the current label, set the scores
        for i in range(len(test_set)):
            curr_score[test_set[i]] = predict[i]

        # add the scores produced by predicting on the current label of test set to a combined score matrix
        scores[int(l)] = curr_score

        # keep track of running time per term, and add it to the total running time
        alg_name = "%s%s" % (alg, run_obj.params_str)
        params_results["%s_process_time"%alg_name] += process_time
        #break

    run_obj.goid_scores = scores
    run_obj.params_results = params_results

