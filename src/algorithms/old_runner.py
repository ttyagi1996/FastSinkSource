
import sys
from collections import defaultdict
import src.setup_sparse_networks as setup
import src.algorithms.alg_utils
import src.algorithms.fastsinksource_runner as fastsinksource
import src.algorithms.genemania_runner as genemania
import src.algorithms.apt_birg_rank_runner as birgrank
import src.algorithms.sinksource_bounds
from src.algorithms.aptrank_birgrank.birgrank import birgRank
import src.algorithms.aptrank_birgrank.run_birgrank as run_birgrank
import src.algorithms.ranks_runner as ranks
import src.algorithms.hemdag_runner as hemdag
#import src.algorithms.sanity_hemdag as sanity
import src.algorithms.svm_runner as svm
import src.algorithms.digraph_fss_runner as digraph_fss
#import src.algorithms.semantic_similarity_runner as itss
import src.algorithms.fss_parent_child_runner as parent_child
import src.algorithms.genemania_parent_child_runner as genemania_pc
import numpy as np
from scipy import sparse


LibMapper = {
    'fastsinksource': fastsinksource,
    'genemania': genemania,
    'aptrank': birgrank,
    'birgrank': birgrank,
    'ranks' : ranks,
    'hemdag' : hemdag,
    #'sanity' : sanity,
    'svm' : svm,
    'digraph_fss' : digraph_fss,
    'parent_child' : parent_child,
    'genemania_pc' : genemania_pc,
    #'itss' : itss,
}


#AlgorithmMapper = {
#    'fastsinksource': fastsinksource.run,
#    'genemania': genemania.run,
#}


#OutputParser = {
#    'fastsinksource': fastsinksource.setupOutputs,
#    'genemania': genemania.setupOutputs,
#}


class Runner(object):
    '''
    A runnable analysis to be incorporated into the pipeline
    '''
    def __init__(self, name, net_obj, ann_obj, out_dir, params, **kwargs):
        self.name = name
        self.net_obj = net_obj
        self.ann_obj = ann_obj
        self.out_dir = "%s/%s/" % (out_dir, name)
        params.pop('should_run', None)  # remove the should_run parameter
        self.params = params
        self.kwargs = kwargs
        self.verbose = kwargs.get('verbose', False) 
        self.forced = kwargs.get('forcealg', False) 
        # for term-based algorithms, can limit the goids for which they will be run
        self.goids_to_run = kwargs.get('goids_to_run', ann_obj.goids)

        # track measures about each run (e.g., running time)
        self.params_results = defaultdict(int) 
        # store the node scores for each GO term in a sparse matrix
        self.goid_scores = sparse.csr_matrix(ann_obj.ann_matrix.shape, dtype=np.float)

        # keep track of the weighting method for writing to the output file later
<<<<<<< HEAD
        self.weight_str = '%s%s%s' % (
            'unw-' if net_obj.unweighted else '', 
            'gm2008-' if net_obj.weight_gm2008 else '',
            'swsn-' if net_obj.weight_swsn else '')
        #self.setupParamsStr()
        self.setupParamsStr(self.weight_str, params, name)
=======
        self.setupParamsStr(net_obj.weight_str, params, name)


>>>>>>> master
    # if the method is not in Python and needs to be called elsewhere, use this
    def setupInputs(self):
        LibMapper[self.name].setupInputs(self)

    # run the method
    def run(self):
        return LibMapper[self.name].run(self)

    # if the method is not in Python and was called elsewhere (e.g., R), 
    # then parse the outputs of the method
    def setupOutputs(self):
        LibMapper[self.name].setupOutputs(self)

    # setup the params_str used in the output file
    def setupParamsStr(self, weight_str, params, name):
        self.params_str = LibMapper[self.name].setup_params_str(weight_str, params, name)

def get_runner_params_str(name, dataset, params):
    """
    Get the params string for a runner without actually creating the runner object
    """
    # get the weight str used when writing output files
    unweighted = dataset['net_settings'].get('unweighted', False) if 'net_settings' in dataset else False
    weight_method = ""
    if dataset.get('multi_net',False) is True: 
        weight_method = "-"+dataset['net_settings']['weight_method']
    weight_str = '%s%s' % (
        '-unw' if unweighted else '', weight_method)

    return LibMapper[name].setup_params_str(weight_str, params, name=name)

