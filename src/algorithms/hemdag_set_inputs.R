#! /home/ttyagi007/R/bin/Rscript

#install.packages("tidyr", repos = "http://cran.us.r-project.org")
#setwd('../../../hpo-prediction');
library(tidyr)
write_scores <- function(filename)
{
    print("Reading the table of flat scores")
    m <- as.matrix(read.table(filename, colClasses="character", stringsAsFactors=FALSE));
    print(dim(m))
    print("making annotations matrix");
    #prot_nodes <- as.character(sort((unique(as.vector(m[,2])))));
    
    prot_nodes <- as.character(seq(0, 17838, 1))
    hpo_nodes <- as.character(sort((unique(as.vector(m[,1])))));

    p.nodes <- length(prot_nodes);
    h.nodes <- length(hpo_nodes);
    
    print(p.nodes)
    print(h.nodes)
    W <- matrix(0, nrow=p.nodes, ncol=h.nodes);

    dimnames(W)[1] <- list(prot_nodes);
    dimnames(W)[2] <- list(hpo_nodes);
    W[cbind(m[,2], m[,1])] <- as.numeric(m[,3]);
    print(dim(W))
    save(W, file="/home/ttyagi007/projects/hpo/hpo-prediction/data/HEMDAG/flat_scores.rda");

}


make_undir_adjMat <- function(filename)
{
        print("making annotations matrix");
        m <- as.matrix(read.table(filename));
        #print(head(m, n=1))
        prot_nodes <- as.character(sort(as.numeric(unique(as.vector(m[,2])))));
        #prot_nodes <- as.character(as.vector(seq(0, 17838, 1)))
        #print(prot_nodes)
        hpo_nodes <- as.character(sort(as.numeric(unique(as.vector(m[,1])))));

        p.nodes <- length(prot_nodes);
        h.nodes <- length(hpo_nodes);
        
        print(prot_nodes)
        print(hpo_nodes)        
        W <- matrix(0, ncol=p.nodes, nrow=h.nodes);

        dimnames(W)[1] <- list(hpo_nodes);
        dimnames(W)[2] <- list(prot_nodes);
        
        #print(dim(W))
        W[cbind(m[,1], m[,2])] <- as.numeric(m[,3]);

        W = t(W)
        print(dim(W))
        #W = gather_matrix(W)
        save(W, file="./data/HEMDAG/flat_scores.rda");


    }

mat <- function(filename)

{
    print("reading into a dataframe for matrix")
    m <- as.data.frame.matrix(read.table(filename));
    
    print("reading into matrix for dimnames")
    m1 <- as.matrix(read.table(filename));
    prot_nodes <- as.character(sort(as.numeric(unique(as.vector(m1[,2])))));
    hpo_nodes <- as.character(sort(as.numeric(unique(as.vector(m1[,1])))));
    
    print("creating the matrix")
    W <- matrix(0, nrow=length(prot_nodes), ncol=length(hpo_nodes));
    dimnames(W)[1] <- list(prot_nodes)
    dimnames(W)[2] <- list(hpo_nodes)

    apply(m, 1, function(row) 
    {
        hpo <- row["V1"]
        gene <- row["V2"]
        W[gene,hpo] <- as.numeric(row["V3"]);

    })
    save(W, file="./data/HEMDAG/ann.rda");

}


updated_mat <- function(filename, outfile)
{   

    print("reading into a matrix to get dim names")
    m1 <- as.matrix(read.table(filename));
    
    
    print("reading the table")
    

    m <- as.data.frame.matrix(read.table(filename));
    prot_nodes <- as.character(sort(as.numeric(unique(as.vector(m1[,2])))));
    hpo_nodes <- as.character(sort(as.numeric(unique(as.vector(m1[,1])))));

    print("constructing the matrix")
    W = spread(m, V1, V3)
    print("saving to file")
    W <- as.matrix(W);
    print(dim(W))
    W <- W[, -c(1)]
    print(length(prot_nodes))
    print(length(hpo_nodes))
    dimnames(W)[1] <- list(prot_nodes)
    dimnames(W)[2] <- list(hpo_nodes) 
    save(W, file=outfile);
    #save(W, file="./data/HEMDAG/ann.rda");


}

print("saving annotations file in rda form")
if(!file_test("-f", './data/HEMDAG/ann.rda'))
{
updated_mat('./data/HEMDAG/ann.txt', "./data/HEMDAG/ann.rda")
#mat('./s_ann_check.txt')
}


print("saving scores file in rda form")
updated_mat('./data/HEMDAG/fss_scores.txt', "./data/HEMDAG/flat_scores.rda")
#write_scores("/home/ttyagi007/projects/hpo/hpo-prediction/data/HEMDAG/fss_scores.txt")



