import tqdm
from tqdm import tqdm, trange
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

def setupInputs(run_obj):
    run_obj.ann_matrix = run_obj.ann_obj.ann_matrix
    run_obj.goids = run_obj.ann_obj.goids
    run_obj.prots = run_obj.ann_obj.prots
    run_obj.hpoidx = run_obj.ann_obj.goid2idx
    run_obj.protidx = run_obj.ann_obj.node2idx
    #run_obj.net_mat = run_obj.net_obj.W
    run_obj.P = alg_utils.normalizeGraphEdgeWeights(run_obj.net_obj.W, ss_lambda=run_obj.params.get('lambda', None))

    run_obj.net_mat = run_obj.P

    return

def setup_params_str(weight_str, params, name):
    return "linear_svm-2000ter"

def setupOutputs(run_obj):
    return

def run(run_obj):
    

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

    '''
    test_set = set()
    
    for i in range(labels.shape[1]):
        test_pos, test_neg = alg_utils.get_goid_pos_neg(test_mat, i)
        for i in test_pos:
            test_set.add(i)
        for i in test_neg:
            test_set.add(i)
        break

    #print(len(test_set))

    
    # creating the training and test feature vectors
    print("Feature vector dimensions")
    X_train = np.delete(feats.toarray(), list(test_set), 0)
    X_train = sparse.lil_matrix(X_train)
    print(X_train.shape)
    
    X_test = feats.toarray()[list(test_set)]
    X_test = sparse.lil_matrix(X_test)
    print(X_test.shape)
    




    # creating the training and test label matrices
    print("Labels dimensions")
    y_train = np.delete(train_mat.toarray().transpose(), list(test_set), 0)
    y_train = sparse.lil_matrix(y_train)
    print(y_train.shape)

    y_test = test_mat.toarray().transpose()[list(test_set)]
    y_test = sparse.lil_matrix(y_test)
    print(y_test.shape)
    
    '''

    # first just try multiclass classification with the first label
    scores = sparse.lil_matrix(labels.transpose().shape, dtype=np.float)        #   dim: hpo x genes
    
    combined_scores = sparse.lil_matrix(labels.transpose().shape, dtype=np.float) # dim: hpo x genes terms
    #print("Shape of combined scores: {}".format(combined_scores.shape))
    #test_set = list(test_set)
    for l in tqdm(range(labels.shape[1])):
        #print("******************working on label........: {}".format(l))
        
        # compute the test gene indices of the annotations for the given label
        
        train_pos, train_neg = alg_utils.get_goid_pos_neg(train_mat,l)
        if len(train_pos)==0:
            print("Skipping term, 0 positive examples")
            continue
        
        
        
        
        #test_set = set()
        test_pos, test_neg = alg_utils.get_goid_pos_neg(test_mat, l)
        
        test_set = set(test_pos) | set(test_neg)
        
        #for i in test_pos:
        #    test_set.add(i)
        #for i in test_neg:
        #    test_set.add(i)
       
        test_set = sorted(list(test_set))
        
        # remove these indices from the training data
        X_train = np.delete(feats.toarray(), test_set, 0)
        X_train = sparse.lil_matrix(X_train)
        #print(X_train.shape)

        # construct the training set with only these indices
        X_test = feats.toarray()[test_set]
        X_test = sparse.lil_matrix(X_test)
        #print(X_test.shape)

        # construct the training label data without these indices

        y_train = np.delete(train_mat.toarray().transpose(), test_set, 0)
        y_train = sparse.lil_matrix(y_train)
        #print(y_train.shape)

        #test_set = list(test_set)
        #print(test_set)
        #print(len(test_set))
        #test_set = sorted(test_set)
        #print(len(test_set))
        # construct the classifier
        temp = []
        
        # this (SVC) uses internal 5-FCV when probability option is set to true
        #clf = SVC(probability=True, verbose=True, gamma='auto', kernel='linear')
    
        svm = LinearSVC(max_iter=2000)
        clf = CalibratedClassifierCV(svm) 

        # get the column of training data for the given label 
        lab = y_train[:,l].toarray().flatten()

        # now train the model on the constructed training data and the column of labels
        #print("Training")

        process_time = time.process_time()

        clf.fit(X_train.toarray(), lab)
        #print(clf.classes_)
        
        # make predictions on the constructed training set
            
        #print("Predicting")
        predict = clf.predict_proba(X_test.toarray())[:,2]
        
        process_time = time.process_time() - process_time

        predict = predict.tolist()
        
        #print(predict)
        
        #print("Length of the predicted array : {}".format(len(predict)))
        
        
        # get the current scores for the given label l
        curr_score = scores[l].toarray().flatten()
         
        # for the test indices of the current label, set the scores
        for i in range(len(test_set)):
            curr_score[test_set[i]] = predict[i]

        #combined_scores[l] = curr_score

        # add the scores produced by predicting on the current label of test set to a combined score matrix
        scores[int(l)] = curr_score

        # keep track of running time per term, and add it to the total running time
        alg_name = "%s%s" % (alg, run_obj.params_str)
        params_results["%s_process_time"%alg_name] += process_time
        #break

    #print(scores.shape)
    run_obj.goid_scores = scores
    run_obj.params_results = params_results
    #print("Number of non-zero values in the scores matrix: %s" %run_obj.goid_scores.count_nonzero())
    #print("Shape of scores matrix after this fold: {}".format(run_obj.goid_scores.shape))

    #print(run_obj.goid_scores)
    #print(scores)
