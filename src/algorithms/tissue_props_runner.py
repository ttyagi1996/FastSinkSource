
import time
import src.algorithms.fastsinksource as fastsinksource
import src.algorithms.alg_utils as alg_utils
from tqdm import tqdm, trange
from scipy import sparse as sp
import numpy as np
from numpy import load
import pandas as pd
import time


pd.set_option('display.float_format', lambda x: '0.20f' %x)

def setupInputs(run_obj):
    # may need to make sure the inputs match
    ## if there are more annotations than nodes in the network, then trim the extra pos/neg nodes
    #num_nodes = self.P.shape[0] if self.weight_gmw is False else self.normalized_nets[0].shape[0]
    #if len(self.prots) > num_nodes: 
    #    positives = positives[np.where(positives < num_nodes)]
    #    negatives = negatives[np.where(negatives < num_nodes)]

    # extract the variables out of the annotation object
    run_obj.ann_matrix = run_obj.ann_obj.ann_matrix
    run_obj.goids = run_obj.ann_obj.goids

    run_obj.edge_density = []
    run_obj.segregation = []
    run_obj.num_genes = []
    if run_obj.net_obj.weight_swsn:
        # TODO if the net obj already has the W_SWSN object, then use that instead
        W, process_time = run_obj.net_obj.weight_SWSN(run_obj.ann_matrix)
        run_obj.P = alg_utils.normalizeGraphEdgeWeights(W, ss_lambda=run_obj.params.get('lambda', None))
        run_obj.params_results['%s_weight_time'%(run_obj.name)] += process_time
        run_obj.W = W
    elif run_obj.net_obj.weight_gmw:
        # this will be handled on a GO term by GO term basis
        run_obj.P = None
    else:
        run_obj.P = alg_utils.normalizeGraphEdgeWeights(run_obj.net_obj.W, ss_lambda=run_obj.params.get('lambda', None))
        run_obj.W = run_obj.net_obj.W

    return


# setup the params_str used in the output file
def setup_params_str(weight_str, params, name="fastsinksource"):
    # ss_lambda affects the network that all these methods use
    """
    ss_lambda = params.get('lambda', 0)
    params_str = "%s-l%s" % (weight_str, ss_lambda)
    if name.lower() not in ["local", "localplus"]:
        a, eps, maxi = params['alpha'], params['eps'], params['max_iters']
        params_str += "-a%s-eps%s-maxi%s" % ( 
            str_(a), str_(eps), str_(maxi))
    """
    params_str=""
    return params_str


def get_alg_type():
    return "term-based"


def setupOutputFile(run_obj):
    return


# nothing to do here
def setupOutputs(run_obj, **kwargs):
    return


def run(run_obj):
    """
    Function to run FastSinkSource, FastSinkSourcePlus, Local and LocalPlus
    *goids_to_run*: goids for which to run the method. 
        Must be a subset of the goids present in the ann_obj
    """

    start = time.process_time()

    hpoids = run_obj.goids
    ann_matrix = run_obj.ann_matrix
    net_matrix = run_obj.P
    
    df = pd.DataFrame()
    for i in trange(ann_matrix.shape[0]):
        
        # get the annotation row corresponding to the term
        y = ann_matrix[i,:]

        # indices of the genes that are positively annotated to the current hpo term
        positives_list = list((y>0).nonzero()[1])
        
        # indices of genes that are left out (not pos)
        left_out_list = list((y<1).nonzero()[1])


        #print(len(positives_list), len(left_out_list), net_matrix.shape[0])

        # number of genes annotated to this term
        num_genes = len(positives_list)
        
        #print(num_genes)
        # obtain the grid for a submatrix of the genes that are positively annotated to this term
        ixrange = np.ix_(positives_list, positives_list)

        # get the submatrix
        sub_net_matrix = net_matrix[ixrange]
        
        # get the sum of the edges in the submatrix
        edge_sum = sub_net_matrix.sum()

        # get the total number of possible edges given the number of genes
        total_edge_sum = num_genes*((num_genes-1)/2)
        
        # compute the edge density
        edge_density = edge_sum / total_edge_sum
        
        run_obj.num_genes.append(num_genes)
        run_obj.edge_density.append(edge_density)
         
        #print(len(positives_list), len(left_out_list))
        seg_range = np.ix_(positives_list, left_out_list)

        # get the submatrix
        
        seg_sub_net_matrix = net_matrix[seg_range]
        
        #print(seg_sub_net_matrix.shape)
        # get the sum of the edges in the left out submatrix

        left_out_sum = seg_sub_net_matrix.sum()
         
        #print(edge_sum, left_out_sum)
        
        if left_out_sum != 0:
            segregation = edge_sum/left_out_sum
        else:
            segregation = 0
        
        #segregation = 0 if left_out_sum == 0 else edge_sum/left_out_sum
    
        

        #print(edge_sum, total_edge_sum)
        #print(edge_sum, left_out_sum)
        segregation = edge_sum/left_out_sum
        run_obj.segregation.append(segregation)
        
    
    print("Time taken to compute properties: {}".format(time.process_time() - start))

    hpo_name = pd.read_csv('outputs/geneset_properties/pos-neg-10-summary-stats.tsv', sep='\t')


    df['HPO Term Name'] = hpo_name['HPO term name']
    df['HPO Id'] = run_obj.goids
    df['Number of genes'] = run_obj.num_genes
    df['Edge Density'] = run_obj.edge_density
    df['Segregation'] = run_obj.segregation



    df.to_csv('outputs/geneset_properties/liver_props-feb11.tsv', index=False, sep='\t')
    #df.to_csv('outputs/peripheral_nervous_system_props.tsv', index=False, sep='\t')
    #print(run_obj.edge_density)
    #run_obj.goid_scores = goid_scores
    #run_obj.params_results = params_results
    
    return


def str_(s):
    return str(s).replace('.','_')
