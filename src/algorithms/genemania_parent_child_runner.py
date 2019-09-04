
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
import src.algorithms.genemania as genemania

def setupInputs(run_obj):
    # get the annotation matrix

    # if running in prediction mode
    run_obj.orig_ann = run_obj.ann_obj.ann_matrix

    # stores the annotation matrix of the current fold if cross validation phase, and full ann if prediction phase
    run_obj.ann_matrix = run_obj.ann_obj.ann_matrix

    # get the hpo ids
    run_obj.hpoids = run_obj.ann_obj.goids

    # get the normalized adjacency matrix
    run_obj.P = alg_utils.normalizeGraphEdgeWeights(run_obj.net_obj.W, ss_lambda=None)

    # get the mappings from hpo terms to the idx being used
    run_obj.hpo2idx = run_obj.ann_obj.goid2idx

    # get the proteins being used in the network
    run_obj.prots = run_obj.net_obj.nodes

    if run_obj.params['matching_hpo']:
        run_obj.parent_full_ann = run_obj.orig_ann
        run_obj.child_full_ann = run_obj.orig_ann
    
    run_obj.neg_erased = 0
    run_obj.pos_erased = 0
    return



# setup the params_str used in the output file
def setup_params_str(weight_str, params, name="fastsinksource"):
    params_str = "%s-tol%s" % (
          weight_str, str(params['tol']).replace('.','_'))
    w = params['weight']

    params_str += "--w%s" % (str_(w))
    return params_str

def setupOutputFile(run_obj):
    return


# nothing to do here
def setupOutputs(run_obj):
    return

def callGenemania(term, net_matrix, term_ann, run_obj):

    # set the graph laplacian for genemania

    L = genemania.setup_laplacian(net_matrix)

    params_results = run_obj.params_results
    goid_scores = run_obj.goid_scores
    alg = run_obj.name

    # run GeneMANIA on each GO term individually
    #for i in tqdm(range(pc_ann_matrix.shape[0])):
    goid = run_obj.hpoids[term]
    # get the row corresponding to the current goids annotations
    #y = ann_matrix[term,:]
    y = term_ann
    # now actually run the algorithm
    scores, process_time, wall_time, iters = genemania.runGeneMANIA(L, y, tol=float(run_obj.params['tol']), verbose=run_obj.kwargs.get('verbose', False))
    tqdm.write("\t%s converged after %d iterations " % (alg, iters) +
            "(%0.3f sec, %0.3f wall_time) for %s" % (process_time, wall_time, goid))
        

        
        #print(goid_scores[i].shape, scores.shape)
    goid_scores = scores[:goid_scores.shape[1]]

    # also keep track of the time it takes for each of the parameter sets
    params_results["%s_wall_time"%alg] += wall_time
    params_results["%s_process_time"%alg] += process_time

    return goid_scores

def check_incoming_weight(net_matrix, params):
    #out_file = params['out_file']
    #with open(out_file, 'a') as :
    #    continue
    return


def run(run_obj):
    ann_matrix = run_obj.ann_matrix
    net_matrix = run_obj.net_obj.W
    goid_scores = run_obj.goid_scores 
    params = run_obj.params
    hpoidx = run_obj.hpo2idx
    alg = run_obj.name
    print("Running correct script")

    pairs = {}

    par = []
    child = []
   
    params_results = run_obj.params_results
    a, eps, max_iters = params['alpha'], float(params['eps']), params['max_iters']

    with open(params['matchings_file'], "r") as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            #print(len(hpoidx))
            #print(hpoidx)
            if(params['matching_hpo']) :
                if row[0] in hpoidx.keys() and row[1] in hpoidx.keys():
                    pairs[hpoidx[row[0]]] = hpoidx[row[1]]
                    par.append(hpoidx[row[1]])
                    child.append(row[0])
            #else:
                # set the value of corresponding key to be index mapping for GO term   
            if params['itss_scores'] == True:
                degree_pair[hpoidx[row[0]]] = float(row[2])

    print(len(pairs)) 
    if params['matching_hpo']:
        only_parents = list(set(hpoidx.keys()) - set(child))
    print(only_parents)
    print(len(only_parents))
    for idx, c in enumerate(pairs.keys()):
        print(idx)
        
        # store the parent of the given child term
        parent = pairs[c]

        # store the hpoids of the parent and child terms
        if params['matching_hpo']:
            id_par = run_obj.hpoids[parent]
            id_child = run_obj.hpoids[c]
        '''
        else:
            id_par = parent_terms_goid[parent]
            id_child = run_obj.goids[c]
        '''

        print(id_child, id_par)
        
        #TODO: revert back to using the setUpPair function
        # if itss scores are to be used then degree_pair[c] is to be used as the weight of an edge

        #if params['itss_scores']:
        #    weight, deg_nodes = setUpPair(run_obj, ss_lambda=degree_pair[c])

        # calculate weight of edge using the fraction of proteins annotated to the child and the number of proteins annotated to parent
        #elif params['frac_ann_scores']:
        #    weight, deg_nodes = setUpPair(run_obj, ss_lambda = None, calculate_weight = True, child=c, parent=parent)
        
        
        # get positive and negative indices for the parent terms using the pos/neg dictionary constructed above
        # NOTE: the parent terms annotations for cross validation come from the full annotation matrix
        
        print('Getting parent and child annotations')

        parent_full_ann = run_obj.parent_full_ann
        child_full_ann = run_obj.child_full_ann

        # obtain the parent terms full annotations
        p = parent_full_ann[parent,:]
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


        # get the child terms original annotations
        child_orig = child_full_ann[c,:]

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
                #run_obj.neg_erased += 1

        for pos in child_pos_erased:
            if pos in positives_par:
                positives_par = np.delete(positives_par, np.where(positives_par == pos))
                #run_obj.count += 1
        
 

            


        """ 
    	# creating the new adjacency matrix
        I = np.eye(net_matrix.shape[0])*params['weight']
        #print(I.shape)
        zero = np.zeros(net_matrix.shape)
    	#print(zero.shape)
    	#print(net_matrix.shape)
    	#csr_matrix(np.hstack([net_matrix, I]))
        print('Creating new adjacency matrix')
    	#csr_matrix(np.hstack([net_matrix.todense(), I]))
        p1 = np.hstack([net_matrix.toarray(), I])
        p2 = np.hstack([zero, net_matrix.toarray()])
        W_combined = csr_matrix(np.vstack([p1, p2]))
        #W_combined = csr_matrix(np.vstack([np.hstack([net_matrix.toarray(), I]), np.hstack([zero, net_matrix.toarray()])]))
    
        print('Normalizing network')        
        P = alg_utils.normalizeGraphEdgeWeights(W_combined, ss_lambda=None)
    
        """
        new_size = (net_matrix.shape[0]*2, 2*net_matrix.shape[1])
        W_combined = sparse.lil_matrix(new_size, dtype=np.float)

        print("Setting first quarter of new adjacency matrix")
        for i in tqdm(range(net_matrix.shape[0])):
            W_combined[i, : net_matrix.shape[1]] = net_matrix[i, : net_matrix.shape[1]]

        print("Setting second quarter of the new adjacency matrix")
        for i in tqdm(range(net_matrix.shape[0])):
                W_combined[i, i+net_matrix.shape[1]] = params['weight']


        #NOTE: the 3rd quarter of the new adjacency matrix is set to 0, since we do not consider any edges to the parent network from the child network
        print("Setting the fourth quarter of the new adjacency matrix")
        for i in tqdm(range(net_matrix.shape[0], W_combined.shape[0])):
                # set the fourth quarter which is equal to the original adjacency matrix
                W_combined[i, net_matrix.shape[1] : ] = net_matrix[i-net_matrix.shape[0], : net_matrix.shape[1]]

        #print('Normalizing Matrix')
        #P = alg_utils.normalizeGraphEdgeWeights(W_combined, ss_lambda=None)
        
        print('Creating new positive and negative annotations index lists')
    	# create the new annotations matrix
	# new annotations matrix : childs annotations | parents updated annotations 
        # the indices in the parents positive and negative annotation sets must be incremented by the number of genes

        negatives_par = [i + ann_matrix.shape[1] for i in negatives_par]
        positives_par = [i + ann_matrix.shape[1] for i in positives_par]

        
        # now create the row for the given "child" terms full network

        full_pos = list(set(positives_par) | set(positives))
        full_neg = list(set(negatives_par) | set(negatives))
    
        term_ann = np.zeros(W_combined.shape[1])
        term_ann[full_pos] = 1
        term_ann[full_neg] = -1

        # call genemania for the current term
        # send the new adjacency matrix to the function
        scores =  callGenemania(c, W_combined, term_ann, run_obj)

        goid_scores[c] = scores[:goid_scores.shape[1]]


        #break

    if params['matching_hpo']:
        for term in only_parents:
            t = hpoidx[term]
            p = run_obj.ann_matrix[t,:]
            positives_par = (p > 0).nonzero()[1]
            negatives_par = (p < 0).nonzero()[1]

            goid = run_obj.hpoids[t]

            term_ann = np.zeros(net_matrix.shape[1])
            term_ann[positives_par] = 1
            term_ann[negatives_par] = -1

            # call genemania for the current term

            scores =  callGenemania(t, net_matrix, term_ann, run_obj)

            goid_scores[t] = scores
    
    print(goid_scores.shape)

    # store the scores in the runner object of the algorithm
    run_obj.goid_scores = goid_scores
    run_obj.params_results = params_results
    
    return

def str_(s):
    return str(s).replace('.','_')
