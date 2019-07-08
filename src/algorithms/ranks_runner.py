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
def setup_params_str(run_obj):
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


# from stackoverflow: https://stackoverflow.com/questions/15081858/can-i-use-rpy2-to-save-a-pandas-dataframe-to-an-rdata-file

def save_rdata_file(mat, filename):

    nr,nc = mat.shape
    Br = ro.r.matrix(mat, nrow=nr, ncol=nc)

    ro.r.assign("B", Br)
    r("save(B, file='{}')".format(filename))
    '''
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_df = ro.conversion.py2rpy(df)
    #r_df = r("data.matrix('{}')".format(r_df))
    r.assign("r_df", r_df)
    r("save(r_df, file='{}')".format(filename))
    '''



def conv_ann_to_adjList(run_obj):

    ann_matrix = run_obj.ann_matrix

    print("Making the adjacency list")

    adjList = []

    for i in trange(ann_matrix.shape[0]):
        hpoid_ann = ann_matrix[i,:]
        positives = (hpoid_ann > 0).nonzero()[1]
        for x in positives:
            adjList.append([x, i, 1])

    return adjList

    



    '''
    adjList = []
    ann = run_obj.ann_matrix
    print(ann)
    hpoid2idx = run_obj.goid2idx
    prot2idx = run_obj.prot2idx

    
    for hpo in run_obj.hpoids:
        row_ann = ann[hpoid2idx[hpo]].toarray().flatten()
        for prot in run_obj.prots:
            adjList.append([hpoid2idx[hpo], prot2idx[prot], row_ann[prot2idx[prot]]])

    return adjList
    '''

def write_to_files(mat, filename):

    print("writing to file")
    with open(filename, 'w') as f:
        for pair in mat:
                f.write(' '.join(str(i) for i in pair))
                f.write('\n')
    
def call_ranks():
    
    print("call RANKS")

    # TODO:
    # set an output filename for RANKS here
    # send the string to RANKS via the subprocess call
    print("Calling RANKS")
    subprocess.call ("/home/ttyagi007/projects/hpo/FastSinkSource/src/algorithms/ranks.R", shell=True)
    print("RANKS has completed")
    #setOutput()
    

def call_set_inputs():
    print("setting inputs for RANKS")
    subprocess.call("/home/ttyagi007/projects/hpo/FastSinkSource/src/algorithms/set_inputs.R", shell=True)
    print("inputs for RANKS have been set")



def setOutput(run_obj):
    ro.r['load']('/home/ttyagi007/projects/hpo/FastSinkSource/outputs/Scores.eav.score.p1.a2.hpo.all.rda')
    scores = ro.r['S']
    scores = np.array(scores)
    

    scores = scores.transpose()


    sparse_scores = sparse.lil_matrix(scores)
    #print(sparse_scores[:5,:])
    


    run_obj.goid_scores = sparse_scores
    run_obj.params_results = "cv_ranks_scores"

    # delete the scores matrix once it has been processed and sent to the main trunk of the pipeline
    #os.remove("/home/ttyagi007/projects/hpo/FastSinkSource/outputs/Scores.eav.score.p1.a2.hpo.all.rda")
    #os.remove("/home/ttyagi007/projects/hpo/hpo-prediction/data/ann.rda")
    #os.remove("/home/ttyagi007/projects/hpo/hpo-prediction/data/ann.txt")

def net_file_to_ids(run_obj, filename, outfile):
    prot2idx = run_obj.prot2idx
    net_file = filename
    mat = []
    with open(net_file, 'r') as f:
        for line in f:
            if line[0] == "#":
                continue
            #u,v,w = line.rstrip().split('\t')[:3]
            line = line.rstrip().split('\t')
            u,v,w = line[:3]
            l = [prot2idx[u], prot2idx[v], w]
            mat.append(l)
            #break
    #print(l)
    #write_to_files(mat, "/home/ttyagi007/projects/hpo/hpo-prediction/data/net.txt") 	
    write_to_files(mat, outfile)


def run(run_obj):

    # TODO: make the paths dynamic instead of being fixed
    
    # change the annotation matrix to only be composed of 0's and 1's

    '''
    run_obj.ann_matrix[run_obj.ann_matrix == -1] = 0
    ann = run_obj.ann_matrix.todense().transpose()
    ann = pd.DataFrame(ann)
    #print(run_obj.ann_matrix)
    print(ann)
    ann.to_csv("/home/ttyagi007/projects/hpo/FastSinkSource/data/ann.txt", header=True, index=True, sep='\t', mode='a')
    '''
    
    #wadj = run_obj.net_matrix.todense()
    #wadj = pd.DataFrame(wadj)

    print("writing adjacency matrix to file")
    
    net_filename = "/home/ttyagi007/projects/hpo/hpo-prediction/inputs/2017_10-string/2017_10-string-net.txt"
    net_outfile = Path("/home/ttyagi007/projects/hpo/hpo-prediction/data/net.txt")

    if(not net_outfile.exists()):
        net_file_to_ids(run_obj, net_filename, net_outfile)
	
    #save_rdata_file(wadj, "/home/ttyagi007/projects/hpo/FastSinkSource/data/net.rda")

    #wadj.to_csv("/home/ttyagi007/projects/hpo/FastSinkSource/data/wadj.txt", header=True, index=True, sep=' ', mode='a')
    '''
    New: 
    - Write the sparse matrices to text files
    - In the r script, implement a function to read the text files, and build the dense matrices
    - These matrices will then have to be written to .rda files since that's what RANKS needs
    Hopefully this makes the process faster
    '''

    # get the ann matrix as pairs of edges and a score

    ann_file = Path("/home/ttyagi007/projects/hpo/hpo-prediction/data/ann.txt")

    if(not ann_file.exists()):
        print("construct triples")
        adjList = conv_ann_to_adjList(run_obj)
    
        print("write the triples to file")
        write_to_files(adjList, ann_file)
    #write_to_files(run_obj.net_matrix,  "/home/ttyagi007/projects/hpo/FastSinkSource/data/net.txt")
    
    '''
    - Old method of converting to dense matrices -> dataframe -> r dataframe -> r matrix -> to file
    
    dense_ann = pd.DataFrame(run_obj.ann_matrix.todense().transpose())
    dense_net = pd.DataFrame(run_obj.net_matrix.todense())

    # make sure that the shapes are correct
    print(dense_ann.shape)
    print(dense_net.shape)


    # save the annotation matrix to .rda format to be used by RANKS
    print("saving the annotation matrix to file")
    save_rdata_file(dense_ann, "/home/ttyagi007/projects/hpo/FastSinkSource/data/ann.rda")

    # save the network matrix to .rda format to be used by RANKS
    print("saving the network matrix to file")
    #save_rdata_file(dense_net, "/home/ttyagi007/projects/hpo/FastSinkSource/data/net.rda")
    
    '''

    # call R script to set inputs for ranks
    call_set_inputs()
    # now we call RANKS
    call_ranks()
    setOutput(run_obj)
    out_file = "home/ttyagi007/projects/hpo/FastSinkSource/outputs/ranks_cv-folds.txt"

    return
    

    
