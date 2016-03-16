#===========================================================================
# DETECTION OF POLYGENIC SELECTION IN GENESETS (POLYSEL)
#
# MODULE: polysel_shufflesets_cluster.R
#
# Use this module on a cluster to test (with pruning) randomized 
# gene sets and create a null distribution of p-values
#
# Created by: Josephine Daub, 
#             CMPG lab, University of Bern, Switzerland
#             Department of Ecology and Evolution, University of Lausanne, 
#             Switzerland
#             Current address: Institute for Evolutionary Biology, UPF-CSIC, 
#             Barcelona, Spain
# Created: March 16, 2016
#
#===========================================================================

# Get parameters from commandline
args <- commandArgs(TRUE)
runnr<-as.character(as.numeric(args[1]))
nrand<-as.numeric(args[2])
minsetsize<-as.numeric(args[3])
approx.null<-as.logical(args[4])
test<-as.character(args[5])
seq.rnd.sampling<-as.logical(args[6])
use.bins<-as.logical(args[7])
project.txt<-as.character(args[8])

data.path=as.character(args[9])
code.path=as.character(args[10])
emp.fdr.path<-as.character(args[11])

# Load all functions needed for the pipeline
source(file.path(code.path,"polysel.R"))
# Load the necessary objects set.info, set.obj and obj.stat
load(file.path(data.path,"polysel_objects.RData"))

# Test randomized sets
TestShuffledSets(set.info, set.obj, obj.stat, 
                 nrand=nrand, minsetsize=minsetsize, 
                 approx.null=approx.null, test=test,
                 seq.rnd.sampling=seq.rnd.sampling,
                 use.bins=use.bins, out.path=emp.fdr.path, 
                 runnr=runnr, project.txt=project.txt)

