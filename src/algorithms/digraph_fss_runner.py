
import time
import src.algorithms.fastsinksource as fastsinksource
import src.algorithms.alg_utils as alg_utils
from tqdm import tqdm, trange
import matplotlib.pyplot as plt
import numpy as np
import csv
from collections import defaultdict
from scipy import sparse
import networkx as nx
from scipy.sparse import csr_matrix

def setupInputs(run_obj):
    # may need to make sure the inputs match
    ## if there are more annotations than nodes in the network, then trim the extra pos/neg nodes
    #num_nodes = self.P.shape[0] if self.weight_gm2008 is False else self.normalized_nets[0].shape[0]
    #if len(self.prots) > num_nodes: 
    #    positives = positives[np.where(positives < num_nodes)]
    #    negatives = negatives[np.where(negatives < num_nodes)]

    # extract the variables out of the annotation object
    run_obj.ann_matrix = run_obj.ann_obj.ann_matrix
    run_obj.goids = run_obj.ann_obj.goids

    if run_obj.net_obj.weight_swsn:
        W, process_time = run_obj.net_obj.weight_SWSN(run_obj.ann_matrix)
        run_obj.P = alg_utils.normalizeGraphEdgeWeights(W, ss_lambda=None)
        #run_obj.Pchild = alg_utils.normalizeGraphEdgeWeights(W, ss_lambda=run_obj.params.get('lambda', None))
        run_obj.params_results['%s_weight_time'%(run_obj.name)] += process_time
    elif run_obj.net_obj.weight_gm2008:
        # this will be handled on a GO term by GO term basis
        run_obj.P = None
    else:
        run_obj.P = alg_utils.normalizeGraphEdgeWeights(run_obj.net_obj.W, ss_lambda=None)
        #run_obj.Pchild = alg_utils.normalizeGraphEdgeWeights(run_obj.net_obj.W, ss_lambda=run_obj.params.get('lambda', None))

    run_obj.hpo2idx = run_obj.ann_obj.goid2idx
    run_obj.prots = run_obj.net_obj.nodes

    run_obj.neg_erased = 0
    run_obj.count = 0
    return

def compute_degree(W, factor, axis = 1):
    deg = np.asarray(W.sum(axis=axis)).flatten()
    deg = np.divide(1., factor + deg)
    if axis == 1:
        deg = csr_matrix(deg).T
    else:
        deg = csr_matrix(deg)
    return deg


# setup the params_str used in the output file
def setup_params_str(weight_str, params, name="fastsinksource"):
    # ss_lambda affects the network that all these methods use
    ss_lambda = params.get('lambda', 0)
    params_str = "%s-l%s" % (weight_str, ss_lambda)
    if name.lower() not in ["local", "localplus"]:
        a, eps, maxi = params['alpha'], params['eps'], params['max_iters']
        if params['itss_scores'] == True:
            params_str += "-itss-a%s-eps%s-maxi%s" % ( 
                str_(a), str_(eps), str_(maxi))
        elif params['frac_ann_scores'] == True:
            params_str += "-ann-a%s-eps%s-maxi%s" % (
                str_(a), str_(eps), str_(maxi))
        else:
            w = params['weight']
            params_str += "-w%s-a%s-eps%s-maxi%s" % (
                  str_(w), str_(a), str_(eps), str_(maxi))



    return params_str


def setupOutputFile(run_obj):
    return


# nothing to do here
def setupOutputs(run_obj):
    return




def update_degree(run_obj, child, parent):

        orig_ann = run_obj.orig_ann
        child_ann = orig_ann[child,:]
        child_pos = (child_ann > 0).nonzero()[1]

        par_ann = orig_ann[parent,:]
        parent_pos = (par_ann > 0).nonzero()[1]

        degree = len(child_pos)/len(parent_pos)

        return degree

def get_new_adj_matrix(run_obj, ss_lambda):
    if run_obj.net_obj.weight_swsn:
        run_obj.child_W, process_time = run_obj.net_obj.weight_SWSN(run_obj.ann_matrix)
        run_obj.Pchild = alg_utils.normalizeGraphEdgeWeights(run_obj.child_W, ss_lambda=ss_lambda)

    else:
        run_obj.child_W = run_obj.net_obj.W
        run_obj.Pchild = alg_utils.normalizeGraphEdgeWeights(run_obj.child_W, ss_lambda=ss_lambda)

    return



def setUpPair(run_obj, ss_lambda, calculate_weight=False, child=None, parent=None):
    
    weight = ss_lambda
    if ss_lambda is not None:
        # this happens when the option is either itss or a prespecified score
        degree = weight
        # set up new adjacency matrix
        get_new_adj_matrix(run_obj, ss_lambda=weight)           # sets up run_obj.Pchild and run_obj.child_W

        # get the degree of the the nodes of the genes for the child term, to be used for multiplying by score
        deg_nodes = compute_degree(run_obj.child_W, degree).toarray().flatten()

        return weight, deg_nodes

    # calculate weight based on the number of annotations of the child and parent pair
    elif calculate_weight==True:

        # calculate the weight of score carrying edge depending on the number of annotations in the parent and child terms
        weight = update_degree(run_obj, child, parent)
        degree=weight
        get_new_adj_matrix(run_obj, ss_lambda=weight)
        deg_nodes = compute_degree(run_obj.child_W, degree).toarray().flatten()
        return weight, deg_nodes

        



def run(run_obj):
    
    
    
    params_results = run_obj.params_results 
    goid_scores = run_obj.goid_scores 
    P = run_obj.P
    hpoidx = run_obj.hpo2idx
    
    alg = run_obj.name
    params = run_obj.params
    a, eps, max_iters = params['alpha'], float(params['eps']), params['max_iters']
    
    

    orig_ann = run_obj.orig_ann
    
    assert not(np.array_equal(orig_ann.toarray(), run_obj.ann_matrix.toarray())), "The complete annotation matrix and fold matrix are the same!"
 
    # Read the child-parent pairs from the csv file and construct a dictionary where the keys are the child nodes
    
    pairs = {}
    
    par = []
    child = []

    degree_pair = {}

    # TODO: add the file to the config file and use that
    # this file contains the child-parent pairs along with a weight containing the itss score
    with open("./data/itss-child_parent.csv", "r") as f:
        #with open("./data/bpln-child_parent.csv", "r") as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            pairs[hpoidx[row[0]]] = hpoidx[row[1]]
            child.append(hpoidx[row[0]])
            par.append(hpoidx[row[1]])
            degree_pair[hpoidx[row[0]]] = float(row[2])

    print("Number of children terms constructed: {}" .format(len(pairs)))
    
    print("Number of nodes that are only parents: {}" .format(len(set(par) - set(child))))

    assert len(pairs) == run_obj.ann_matrix.shape[0]- (len(set(par) - set(child))), "Number of modified children scores not the same as number of children in pairs"
    
    # used later to add the scores of the parent terms to the matrix
    # this is done since some terms are only parents, and thus are not handled during the parent-child run
    only_parents = list(set(par)-set(child))
    
    # if a prespecified weight is to be used for edges for all hpo term pairs
    if params['itss_scores'] == False and params['frac_ann_scores'] == False:
        degree, deg = setUpPair(run_obj, ss_lambda=run_obj.params.get('weight', None))


    # Run the parent-child model

    # get each child in the child-parent pairs dictionary
    for idx, c in enumerate(pairs.keys()):
        print(idx)
        
        # store the parent of the given child term
        parent = pairs[c]

        # store the goids of the parent and child terms
        goid_par = run_obj.goids[parent]
        goid_child = run_obj.goids[c]

        print(goid_child, goid_par)
        
        
        # if itss scores are to be used then degree_pair[c] is to be used as the weight of an edge

        if params['itss_scores'] == True:
            degree, deg = setUpPair(run_obj, ss_lambda=degree_pair[c])

        # calculate weight of edge using the fraction of proteins annotated to the child and the number of proteins annotated to parent
        elif params['frac_ann_scores'] == True:
            degree, deg = setUpPair(run_obj, ss_lambda = None, calculate_weight = True, child=c, parent=parent)
        #else:
        #    degree, deg = setUpPair(run_obj, ss_lambda=run_obj.params.get('weight', None))

        
        # get positive and negative indices for the parent terms using the pos/neg dictionary constructed above
        # NOTE: the parent terms annotations for cross validation come from the full annotation matrix 
        p = orig_ann[parent,:]
        positives_par = (p > 0).nonzero()[1]
        negatives_par = (p < 0).nonzero()[1]

        copy_par_pos = (p > 0).nonzero()[1]
        copy_par_neg = (p < 0).nonzero()[1]
        
        # get positive and negative indices for the child term
        # the child gets its annotations from the current fold
        y = run_obj.ann_matrix[c,:]
        positives = (y > 0).nonzero()[1]
        negatives = (y < 0).nonzero()[1]
        
         
        # For every negative annotation that has become unknown in the current fold annotations of the child term
            # erase the annotation of the corresponding gene from the parent term as well.
            # This is done since 

        # get the child terms original annotations
        child_orig = orig_ann[c,:]

        # obtain the annotations in the full annotation list of the child term
        child_orig_neg = (child_orig < 0).nonzero()[1]
        child_orig_pos = (child_orig > 0).nonzero()[1]

        # obtain the negative annotations that have been erased in the current fold of cross validation for the child term
        child_neg_erased = list(set(child_orig_neg) - set(negatives))
        child_pos_erased = list(set(child_orig_pos) - set(positives))

        # erase the annotations of these negative gene terms of the child not present in the current fold, from the full annotations of the parent term
        # TODO: convert this to a set operation
        
        # for each erased negative annotation
         
        for neg in child_neg_erased:
            # if it exists in the parents set of negative annotations
            if neg in negatives_par:
                # remove the neg annotations from the parent
                negatives_par = np.delete(negatives_par, np.where(negatives_par == neg))
                run_obj.neg_erased += 1

        for pos in child_pos_erased:
            if pos in positives_par:
                positives_par = np.delete(positives_par, np.where(positives_par == pos))
                run_obj.count += 1
        
             
        assert len(negatives_par) <= len(copy_par_neg), "Number of negative annotations after removing have increased"
        assert len(positives_par) <= len(copy_par_pos), "Number of positive annotations after removing have increased"
        
        #run_obj.neg_erased.append(len(set(copy_par_neg) - set(negatives_par)))
         
            
        # call fastsinksource on the parent term using the positive negative indices from the dictionary/above updated indices

        parent_scores, process_time, wall_time, iters = fastsinksource.runFastSinkSource(
                P, positives_par, negatives=negatives_par, max_iters=max_iters,
                eps=eps, a=a, verbose=run_obj.kwargs.get('verbose', False))
        


        tqdm.write("\t%s Parent Term converged after %d iterations " % (alg, iters) +
                "(%0.4f sec) for %s" % (process_time, goid_par))
        

        # update the P (adjacency matrix) to alter the weights of the edges due to the effect of adding a weight of "w" from each node to itself
        # NOTE: this is now obtained in the setupinput params function, assuming that the weight of every edges carrying a parent score is 1
        # this will need to be made into a function when considering the weight of a score carrying edge to be the probability.

        new_P = run_obj.Pchild
        
        # update the parent scores to be normalized by the degree obtained by considering the notion of adding an edge of weight 1 that represent the parent score

        # get the degree from the function defined above: compute_degree
        # this will have to be removed/changed to enable code reuse
        #print(degree)
        #deg = compute_degree(run_obj.child_W, degree).toarray().flatten()
    
        # now multiply the parent scores by this weighted degree
        # the scores will then be added to the fixed score vector "f" computed in alg_utils

        for i in range(len(parent_scores)):
            #print(parent_scores[i], deg[i], degree)
            parent_scores[i] = parent_scores[i]*deg[i]*degree

        

        # call fastsinksource for the child term and send the scores of the parent term to update the fixed "f" scores
        # the child term FastSinkSource is called with the new normalized adjacency matrix 
        scores, process_time, wall_time, iters = fastsinksource.runFastSinkSource(
                new_P, positives, negatives=negatives, max_iters=max_iters,
                eps=eps, a=a, verbose=run_obj.kwargs.get('verbose', False),  scores = parent_scores)
        

        assert scores.shape[0] == run_obj.ann_matrix.shape[1], "Number of genes in scores doesn't match total number of genes"

        tqdm.write("\t%s Child term converged after %d iterations " % (alg, iters) +
                "(%0.4f sec) for %s" % (process_time, goid_child))
        
        # store the scores of the child term in the matrix
        goid_scores[c] = scores
   
      
   
    # get the scores of all the nodes that are only parent terms
    # for this we use the fold annotation directly.

    for term in only_parents:
        p = run_obj.ann_matrix[term,:]
        positives_par = (p > 0).nonzero()[1]
        negatives_par = (p < 0).nonzero()[1]

        goid = run_obj.goids[term]
        parent_scores, process_time, wall_time, iters = fastsinksource.runFastSinkSource(
                P, positives_par, negatives=negatives_par, max_iters=max_iters,
                eps=eps, a=a, verbose=run_obj.kwargs.get('verbose', False))


        tqdm.write("\t%s Parent Term converged after %d iterations " % (alg, iters) +
                "(%0.4f sec) for %s" % (process_time, goid))

        goid_scores[term] = parent_scores
    
    #print(goid_scores)
    print(goid_scores.shape)

    # store the scores in the runner object of the algorithm
    run_obj.goid_scores = goid_scores
    run_obj.params_results = params_results
    #print(run_obj.neg_erased)
    
    '''
    with open('num_neg_erased.txt','a') as f:
        f.write('\n'.join(str(v) for v in run_obj.neg_erased))
    '''
    print("Total erased: {}".format(run_obj.neg_erased))
    print(run_obj.count)

    return
    


def str_(s):
    return str(s).replace('.','_')
