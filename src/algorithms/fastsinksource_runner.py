

import time
import src.algorithms.fastsinksource as fastsinksource
import src.algorithms.alg_utils as alg_utils
from tqdm import tqdm, trange
from scipy import sparse as sp
from scipy.sparse.linalg import LinearOperator, spilu
import numpy as np


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

    if run_obj.net_obj.weight_swsn:
        # TODO if the net obj already has the W_SWSN object, then use that instead
        W, process_time = run_obj.net_obj.weight_SWSN(run_obj.ann_matrix)
        run_obj.P = alg_utils.normalizeGraphEdgeWeights(W, ss_lambda=run_obj.params.get('lambda', None))
        run_obj.params_results['%s_weight_time'%(run_obj.name)] += process_time
    elif run_obj.net_obj.weight_gmw:
        # this will be handled on a GO term by GO term basis
        run_obj.P = None
    else:
        run_obj.P = alg_utils.normalizeGraphEdgeWeights(run_obj.net_obj.W, ss_lambda=run_obj.params.get('lambda', None))
        #run_obj.P = alg_utils._net_normalize(run_obj.net_obj.W)

    return


# setup the params_str used in the output file
def setup_params_str(weight_str, params, name="fastsinksource"):
    # ss_lambda affects the network that all these methods use
    ss_lambda = params.get('lambda', 0)
    params_str = "%s-l%s" % (weight_str, ss_lambda)
    if name.lower() not in ["local", "localplus"]:
        a, eps, maxi = params['alpha'], params['eps'], params['max_iters']
        tol = "-tol%s" % (str_(params['tol'])) if 'tol' in params else ''
        solver = "-%s" % (params['solver']) if 'solver' in params else ''
        params_str += "-a%s-eps%s-maxi%s%s%s" % ( 
            str_(a), str_(eps), str_(maxi), tol, solver)

    return params_str
