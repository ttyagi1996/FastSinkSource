
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
    return

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


def run(run_obj):
    ann_matrix = run_obj.ann_matrix
    net_matrix = run_obj.net_obj.W
    goid_scores = run_obj.goid_scores 
    params = run_obj.params
    hpoidx = run_obj.hpo2idx
    print("Running correct script")

    pairs = {}

    par = []
    child = []


    with open(params['matchings_file'], "r") as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            if(params['matching_hpo']) :
                if row[0] in hpoidx and row[1] in hpoidx:
                    pairs[hpoidx[row[0]]] = hpoidx[row[1]]
                    par.append(hpoidx[row[1]])
                    child.append(hpoidx[row[0]])
            #else:
                # set the value of corresponding key to be index mapping for GO term   
            if params['itss_scores'] == True:
                degree_pair[hpoidx[row[0]]] = float(row[2])
    
    if params['matching_hpo']:
        only_parents = list(set(par)-set(child))

    # if a prespecified weight is to be used for edges for all hpo term pairs
    #if params['itss_scores'] == False and params['frac_ann_scores'] == False:
    #weight, deg_nodes = setUpPair(run_obj, ss_lambda=run_obj.params.get('weight'))
    
    weight = params['weight']
    

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
        
        parent_full_ann = run_obj.parent_full_ann
        child_full_ann = run_obj.child_full_ann
        #parent_full_ann = run_obj.parent_full_ann
        #child_full_ann = run_obj.child_full_ann

         # obtain the parent terms full annotations
        parent_ann = parent_full_ann[parent,:].toarray().flatten()
        #print(parent_ann[:20])
        
        # get positive and negative indices for the child term
        # the child gets its annotations from the current fold
        y = run_obj.ann_matrix[c,:]
        positives = (y > 0).nonzero()[1]
        negatives = (y < 0).nonzero()[1]



         # get the child terms original annotations
        child_orig = child_full_ann[c,:]

        # obtain the annotations in the full annotation list of the child term
        child_orig_neg = (child_orig < 0).nonzero()[1]
        child_orig_pos = (child_orig > 0).nonzero()[1]

        # obtain the negative annotations that have been erased in the current fold of cross validation for the child term
        child_neg_erased = list(set(child_orig_neg) - set(negatives))
        child_pos_erased = list(set(child_orig_pos) - set(positives))

        # for each erased negative annotation
        for neg in child_neg_erased:
            # if it exists in the parents set of negative annotations
            if parent_ann[neg] == -1 :
                # remove the neg annotations from the parent
                parent_ann[neg] = 0
                run_obj.neg_erased += 1

        for pos in child_pos_erased:
            if parent_ann[pos] == 1:
                parent_ann[pos] = 0
                run_obj.pos_erased += 1



        print(run_obj.pos_erased, run_obj.neg_erased)

        # call genemania for the parent term

        parent_scores = callGenemania(parent, net_matrix, parent_ann, run_obj)

        for i in range(len(parent_scores)):
            parent_scores[i] = parent_scores[i]*weight


        # setup the new graph which consists of parent child networks connected

        #print("Size of annotations matrix: {} x {}" .format(run_obj.ann_matrix.shape[0], run_obj.ann_matrix.shape[1]))
        new_size = (ann_matrix.shape[0], 2*ann_matrix.shape[1])
        pc_ann_matrix = sparse.lil_matrix(new_size, dtype=np.float)
        #print("Size of new annotations matrix: {} x {}" .format(pc_ann_matrix.shape[0], pc_ann_matrix.shape[1]))
    
        #print("Setting the new annotation matrix")
        child_ann = np.ones(pc_ann_matrix.shape[1])
        child_ann[:ann_matrix.shape[1]] = ann_matrix[c,:].toarray().flatten()
        #print(child_ann.shape)
        
        
        #for i in tqdm(range(ann_matrix.shape[0])):
        #    pc_ann_matrix[i, : ann_matrix.shape[1]] = ann_matrix[i, : ann_matrix.shape[1]]
   
        #for i in range(ann_matrix.shape[0]):
        #    pc_ann_matrix[i, ann_matrix.shape[1]: ] = 1        # set all the nodes of the graph corresponding to the parent term to be a "positive" example
    

        #assert(pc_ann_matrix.count_nonzero()-(3320*17839) , ann_matrix.count_nonzero()), "New annotations matrix not set correctly"
    
        #print("Size of the network: {} x {}" .format(net_matrix.shape[0], net_matrix.shape[1]))
        new_size = (net_matrix.shape[0]*2, 2*net_matrix.shape[1])
        pc_graph = sparse.lil_matrix(new_size, dtype=np.float)
        #print("Size of new network matrix: {} x {}" .format(pc_graph.shape[0], pc_graph.shape[1]))
    
    
        #parent_scores = np.zeros(17839)
        
        # if the parent scores are all zero then there is no influence from the network of the parent term, and we only consider the network corresponding to the child term
        if len(np.where(parent_scores==0)[0]) != len(parent_scores):

            
            print("Setting first quarter of new adjacency matrix")
            for i in tqdm(range(net_matrix.shape[0])):
                pc_graph[i, : net_matrix.shape[1]] = net_matrix[i, : net_matrix.shape[1]]
                #pc_graph[i, net_matrix.shape[1]: ] = parent_scores[i]
    


            print("Setting second quarter of the new adjacency matrix")
            for i in tqdm(range(net_matrix.shape[0])):
                pc_graph[i, i+net_matrix.shape[1]] = parent_scores[i]
                # this needs to be set to weight of score carrying edge*parent score
                # pc_graph[i, i] = parent_scores[i-net_matrix.shape[1]]
        
            print("Setting the fourth quarter of the new adjacency matrix")
            for i in tqdm(range(net_matrix.shape[0], pc_graph.shape[0])):
                # set the fourth quarter which is equal to the original adjacency matrix
                pc_graph[i, net_matrix.shape[1] : ] = net_matrix[i-net_matrix.shape[0], : net_matrix.shape[1]]
    
        else:
            pc_graph = run_obj.net_obj.W
            child_ann = y.toarray().flatten()
        
        goid_scores[c] = callGenemania(c, pc_graph, child_ann, run_obj)

    
    # RUN GENEMANIA FOR THOSE TERMS THAT ARE ONLY PARENTS
    if params['matching_hpo']:
        for term in only_parents:
            term_ann = run_obj.ann_matrix[term, :].toarray().flatten()
            goid_scores[term] = callGenemania(term, net_matrix, term_ann, run_obj)





    run_obj.goid_scores = goid_scores
    #run_obj.params_results = params_results





    return

def str_(s):
    return str(s).replace('.','_')
