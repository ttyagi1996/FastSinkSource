#! /home/ttyagi007/R/bin/Rscript
library(HEMDAG);

#setwd("/home/ttyagi007/projects/hpo/hpo-prediction/")

#print(getwd())
#return

loadRDA <- function(filename)
{
        load(filename)
        get(ls()[ls() != "filename"])


}

# convert the list of HPO edges to graphNEL object

print("Reading the ontology edges into a graphNEL object")
g <- read.graph(file="./data/HEMDAG/ontology_ids.txt");


# build ancestors of HPO

print("Building ancestors of the ontology")
anc <- build.ancestors(g);



# load the scores produced by fast sink source

print("Loading the flat scores")
Flat <- loadRDA("./data/HEMDAG/flat_scores.rda")
#S <- t(S)
#print(dim(S))
colnamerange = dim(Flat)[2] - 1

dimnames(Flat)[2] <- list((seq(0, colnamerange,1)))

# normalize the scores

print("Normalizing the flat scores")
Flat <- normalize.max(Flat);

# find the root of HPO

print("Finding root of the ontology")
root <- root.node(g);
print(root)

print("Removing root from flat scores, if exists")

if(root %in% colnames(Flat))
    Flat <- Flat[,-which(colnames(Flat)==root)];

# finally, run hierarchical top down for DAG

print(head(Flat,n=1))

print("Running HTD")
Hem <- htd(Flat,g,root);

print(head(Hem, n=1))
#Hem   <- TPR.DAG(Flat, g, root, positive="descendants", bottomup="weighted.threshold.free", topdown="HTD", w=0.5);



# save the scores in rda format to the data directory

print("Writing the scores to file")
save(Hem, file="./data/HEMDAG/hemdag.rda")


# check hierarchy of the flat vs htd scores
print(check.hierarchy(Flat, g, root)$Status)
print(check.hierarchy(Hem, g, root)$Status)


print("1")
