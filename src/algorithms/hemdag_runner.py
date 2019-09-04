import time
import subprocess
from rpy2 import robjects as ro
from tqdm import tqdm, trange
import pyreadr
import numpy as np
import pandas as pd
from scipy import sparse
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.conversion import localconverter
import os
pandas2ri.activate()
import src.evaluate.eval_utils as eval_utils
import os.path
from pathlib import Path
import src.algorithms.fastsinksource_runner as fastsinksource


def setup_params_str(weight_str,param_str,run_obj):
    return


def setupInputs(run_obj):
    run_obj.ann_matrix = run_obj.ann_obj.ann_matrix
    run_obj.hpoids = run_obj.ann_obj.goids
    run_obj.goid2idx = run_obj.ann_obj.goid2idx
    run_obj.net_matrix = run_obj.net_obj.W
    run_obj.prots = run_obj.net_obj.nodes
    run_obj.prot2idx = run_obj.net_obj.node2idx


def setupOutputs(run_obj):
    return


def setOutput(run_obj):
    ro.r['load']('./data/HEMDAG/htd.rda')
    scores = ro.r['S.htd']
    scores = np.array(scores)


    scores = scores.transpose()


    sparse_scores = sparse.lil_matrix(scores)
    print(sparse_scores[:1,:])
    
    run_obj.goid_scores = sparse_scores
    run_obj.params_str = "cv_HTD_scores"


    # delete the scores matrix once it has been processed and sent to the main trunk of the pipeline
    #os.remove("./data/HEMDAG/htd.rda")
    #os.remove("./data/HEMDAG/fss_scores.txt")
    #os.remove("/home/ttyagi007/projects/hpo/hpo-prediction/data/ann.txt")




def write_to_files(mat, filename):

    print("writing to file")
    with open(filename, 'w') as f:
        for pair in mat:
                f.write(' '.join(str(i) for i in pair))
                f.write('\n')



def net_file_to_ids(run_obj, filename, outfile, ont=0):

    if(ont):
        hpo2idx = run_obj.goid2idx
        net_file = filename
        mat = []
        with open(net_file, 'r') as f:
            for line in f:
                if line[0] == "#":
                    continue
                #   u,v,w = line.rstrip().split('\t')[:3]
                line = line.rstrip().split(' ')
                u,v = line[:2]
                if u != 'HP:0000118':
                    l = [hpo2idx[u], hpo2idx[v]]
                else:
                    l = ['HP:0000118', hpo2idx[v]]
                mat.append(l)
                #break
    else:
        hpoidx = run_obj.goid2idx
        protidx = run_obj.prot2idx
        net_file = filename
        mat = []
        with open(net_file, 'r') as f:
            for line in f:
                if line[0] == "#":
                    continue
                #   u,v,w = line.rstrip().split('\t')[:3]
                line = line.rstrip().split('\t')
                u, v, w = line[:3]
                l = [hpoidx[u], protidx[v], w]
                mat.append(l)

        




    write_to_files(mat, outfile)

def call_set_inputs():
    print("setting inputs for RANKS")
    subprocess.call("/home/ttyagi007/projects/github/FastSinkSource/src/algorithms/hemdag_set_inputs.R", shell=True)
    print("inputs for HEMDAG have been set")


def call_hemdag():
    print("Calling HEMDAG")
    subprocess.call("/home/ttyagi007/projects/github/FastSinkSource/src/algorithms/hemdag.R", shell=True)
    print("HEMDAG finished")


def write_scores_to_file(hpo_scores, hpoids, prots, out_file):
    with open(out_file, 'w') as out:
        for i in range(hpo_scores.shape[0]):
            scores = hpo_scores[i].toarray().flatten()
            scores = {j:s for j, s in enumerate(scores)}
            for n in sorted(scores, key=scores.get, reverse=True)[:-1]:
                out.write("%s %s %0.6e\n" % (i, n, scores[n]))


def conv_ann_to_adjList(run_obj):

    ann_matrix = run_obj.ann_matrix
    hpoid = run_obj.hpoids
    prots = run_obj.prots
    print("Making the adjacency list")

    adjList = []

    for i in trange(ann_matrix.shape[0]):
        hpoid_ann = ann_matrix[i,:]
        positives = (hpoid_ann > 0).nonzero()[1]
        for x in positives:
            adjList.append([hpoid[i], x, 1])

    return adjList



def write_to_files(mat, filename):

    print("writing to file")
    with open(filename, 'w') as f:
        for pair in mat:
                f.write(' '.join(str(i) for i in pair))
                f.write('\n')



def run(run_obj):


    #TODO: how to integrate with pipeline to get the required scores?
    
    # ontology edges conversion: Convert the HPO names to the indices
    # update: generate the ontology separately and use that
    #net_file_to_ids(run_obj, "./data/HEMDAG/hpo10.txt","./data/HEMDAG/ontology_ids.txt", 1)

    
    # call sinksource
    # TODO: add a runner object for sinksource within config for hemdag
    
    run_obj.name = "fastsinksource"
    fastsinksource.setupInputs(run_obj)
    fastsinksource.run(run_obj)
    run_obj.name = "hemdag"
    
    # check if scores are returned
    #print(run_obj.goid_scores[1,:])


    # write the flat scores to a file

    out_file = "./data/HEMDAG/fss_scores.txt"
    write_scores_to_file(run_obj.goid_scores, run_obj.hpoids, run_obj.prots, out_file)


    
    # write the annotation matrix in form of edges in a text file
    ann_file = Path("./data/HEMDAG/ann.txt")

    #if(not ann_file.exists()):
    print("construct triples")
    adjList = conv_ann_to_adjList(run_obj)

    print("write the triples to file")
    write_to_files(adjList, ann_file)
 
    # convert the FSS scores and annotation matrix to rda files
    

    call_set_inputs()
    
    # call hemdag
    call_hemdag()

    # make the scores of the run object to be equal to what hemdag returns
    setOutput(run_obj)
    
    return

