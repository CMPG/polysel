#===========================================================================
# DETECTION OF POLYGENIC SELECTION IN GENESETS (POLYSEL)
#
# MODULE: polysel_cluster.R
#
# Use this module to run a gene set enrichment analysis on a cluster
# NOTE: first run polysel_shufflesets_cluster.R to create a null
# distribution of p-values
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
nrand<-as.numeric(args[1])
minsetsize<-as.numeric(args[2])
approx.null<-as.logical(args[3])
test<-as.character(args[4])
seq.rnd.sampling<-as.logical(args[5])
use.bins<-as.logical(args[6])
do.pruning<-as.logical(args[7])
project.name<-as.character(args[8])
do.emp.fdr<-as.logical(args[9])
emp.fdr.nruns<-as.numeric(args[10]) 
emp.fdr.est.m0<-as.logical(args[11])
qvalue.method<-as.character(args[12])

data.path=as.character(args[13])
code.path=as.character(args[14])
results.path=as.character(args[15])
emp.fdr.path<-as.character(args[16])

# Load all functions needed for the pipeline
source(file.path(code.path,"polysel.R"))
# Load the necessary objects set.info, set.obj and obj.stat
load(file.path(data.path,"polysel_objects.RData"))

# Test randomized sets
r<-EnrichmentAnalysis(set.info, set.obj, obj.stat,
                      nrand=nrand, approx.null=approx.null, 
                      seq.rnd.sampling=seq.rnd.sampling,
                      use.bins=use.bins, test=test,
                      do.pruning=do.pruning, minsetsize=minsetsize,
                      project.txt=project.name, do.emp.fdr=do.emp.fdr, 
                      emp.fdr.path=emp.fdr.path, 
                      emp.fdr.nruns=emp.fdr.nruns, 
                      emp.fdr.est.m0=emp.fdr.est.m0,
                      qvalue.method=qvalue.method)

set.scores.prepruning <- r$set.scores.prepruning
set.scores.postpruning <- r$set.scores.postpruning

save(set.scores.prepruning, 
     file = file.path(results.path, 
                      paste(project.txt=project.name,
                            "_setscores_prepruning_",
                            formatC(nrand,format="d"),".RData",
                            sep="")))

write.table(set.scores.prepruning, quote=FALSE, sep="\t", row.names=FALSE,            
            file = file.path(results.path,
                             paste(project.txt=project.name,
                                   "_setscores_prepruning_",
                                   formatC(nrand,format="d"),".txt",
                                   sep=""))) 

save(set.scores.postpruning, 
     file = file.path(results.path, 
                      paste(project.txt=project.name,
                            "_setscores_postpruning_",
                            formatC(nrand,format="d"),".RData",
                            sep="")))

write.table(set.scores.postpruning, quote=FALSE, sep="\t", row.names=FALSE,            
            file = file.path(results.path,
                             paste(project.txt=project.name,
                                   "_setscores_postpruning_",
                                   formatC(nrand,format="d"),".txt",
                                   sep="")))
