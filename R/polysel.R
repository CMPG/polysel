#===========================================================================
# DETECTION OF POLYGENIC SELECTION IN GENESETS (POLYSEL)
#
# MODULE: polysel.R
#
# This module contains functions to run a gene set enrichment analysis to 
# find signals of polygenic selection on pathways 
#
# Created by: Josephine Daub, 
#             CMPG lab, University of Bern, Switzerland
#             Department of Ecology and Evolution, University of Lausanne, 
#             Switzerland
#             Current address: Institute for Evolutionary Biology, UPF-CSIC, 
#             Barcelona, Spain
# Created: Wednesday 16, 2016
#
#===========================================================================

#===========================================================================
# Function: ReadSetObjTables(in.path, set.info.file,set.obj.file,
#                            obj.info.file)
# Read in all required gene (object) and gene set (set) tables
#
# - in.path      : path to directory with input files
# - set.info.file: tab seperated file with fields: 
#                  setID, setName, ...
# - set.obj.file : tab seperated file with fields: 
#                  setID, objID
# - obj.info.file: tab seperated file with fields: 
#                  objID, objName, objStat, (objBin), (objSNPcnt), ...
# - minsetsize   : exclude gene sets with size below minsetsize
# - maxsetsize   : exclude gene sets with size above maxsetsize
# - obj.in.set   : exclude genes that are not part of a set
#
# These files must contain headers, IDs can be strings
# Internal numeric IDs will be assigned to objects and sets to improve 
# further computations
#===========================================================================
ReadSetObjTables<-function(in.path, set.info.file, set.obj.file,
                           obj.info.file, minsetsize=10, maxsetsize=1000,
                           obj.in.set=FALSE){
  
  # Read in information on gene sets
  set.info<-read.table(file=file.path(in.path,set.info.file),
                       header=T, sep="\t", stringsAsFactors=F,
                       quote = "")  
  set.info$setID.orig<-as.character(set.info$setID)
  set.info$setID<-seq(nrow(set.info)) 
  
  # Read in information on genes
  obj.info<-read.table(file=file.path(in.path,obj.info.file),
                       header=T, sep="\t", stringsAsFactors=F,
                       quote = "")
  obj.info$objID.orig<-as.character(obj.info$objID)
  obj.info$objID<-seq(nrow(obj.info))
  if (!("objBin" %in% names(obj.info))) obj.info$objBin<-1
    
  # Read in information on which genes are in which set
  set.obj<-read.table(file=file.path(in.path,set.obj.file),
                      header=T, sep="\t", stringsAsFactors=F,
                      quote = "")
  set.obj$setID<-as.character(set.obj$setID)
  set.obj$objID<-as.character(set.obj$objID)
  m<-match(set.obj$setID,set.info$setID.orig)
  set.obj$setID<-set.info$setID[m]
  m<-match(set.obj$objID,obj.info$objID.orig)
  set.obj$objID<-obj.info$objID[m]
  
  # Cleaning data:
  
  # Remove sets from set.obj that are not in set.info
  ix<-set.obj$setID %in% set.info$setID
  cat("Found and removed ", sum(!ix), " gene sets in set.obj",
      " which where not in set.info\n", sep="")
  set.obj<-set.obj[ix,]  
  
  # Remove genes from set.obj that are not in obj.info
  ix<-set.obj$objID %in% obj.info$objID
  cat("Found and removed ", sum(!ix), " objects in set.obj",
      " which where not in obj.info\n", sep="")
  set.obj<-set.obj[set.obj$objID %in% obj.info$objID,]
  
  # Remove gene sets that have size outside [minsetsize,maxsetsize]
  t<-table(set.obj$setID)
  ix<-(t<minsetsize | t>maxsetsize)
  cat("Found and removed ", sum(ix), " gene sets in set.obj with size ",
      "outside [",minsetsize, ", ", maxsetsize,"]\n",
      sep="")
  set.obj<-set.obj[set.obj$setID %in% names(t)[!ix],]
  
  # Exclude genes that are not part of a set
  if (obj.in.set){
    ix<-obj.info$objID %in% set.obj$objID
    cat("Found and removed ", sum(!ix), " objects in obj.info",
        " which where not in set.obj\n", sep="")
    obj.info<-obj.info[ix,]
  }
  
  return(list(set.info=set.info,set.obj=set.obj,obj.info=obj.info))
}


#===========================================================================
# Function: CheckStatDistribution(obj.stat, n.rand, setsize, bin)
# Check whether objStat follows normal distribution
#
# - obj.stat  : dataframe with fields objID and objStat
# - n.rand    : number of randomizations
# - setsize   : number of genes in set
# - show.stats: show expected and real mean and variance?
# - bin       : only show distribution of genes in this bin
#===========================================================================
CheckStatDistribution<-function(obj.stat, n.rand=5000, setsize=50, 
                                xlab="objStat", show.stats=F,
                                bin){

  # if bin is given restrict to genes in bin
  if (!missing(bin))
    stats<-obj.stat$objStat[obj.stat$objBin==bin]
  else 
    stats<-obj.stat$objStat

  # calculate mean and sd of expected normal distribution
  n.obj<-length(stats)
  stat.mean<-setsize*mean(stats)
  stat.var<-setsize*var(stats)  
  stat.sd<-sqrt(stat.var)    

  # create null distribution
  scores<-rep(0,n.rand)
  for (i in 1:n.rand){
    ix<-sample(1:n.obj,setsize)
    scores[i]<-sum(stats[ix])
  }
  
  # plot null distribution
  main.txt<-paste("setsize = ",setsize,sep="")
  if (!missing(bin))
    main.txt<-paste(main.txt, ", bin = ", bin, sep="")
  hist(scores,freq=F,breaks=200, main=main.txt,
       xlab=xlab)
  
  # plot normal distribution 
  x<-seq(stat.mean-4*stat.sd,stat.mean+4*stat.sd,length=1000)
  y<-dnorm(x,mean=stat.mean,sd=sqrt(stat.var))
  lines(x,y,col="red")
  abline(v=mean(scores),col="blue")
  
  legend("topright", legend=c("estimated normal distr.",
                              "mean SUMSTAT scores"),
         col=c("red","blue"), lty=1, bty="n", cex=0.9)
  
  if (show.stats){
    cat(main.txt, "\n", sep="")
    cat("obs mean:",mean(scores), "\nexp mean:",stat.mean,
        "\n\nobs var:",var(scores),"\nexp var:",stat.var,"\n\n",sep="")
  }
}

#===========================================================================
# Function: PlotStatField(obj.info, fld, xlab, ylab, logaxis, show.bins)
# Plot the correlation between original score of a gene and the number 
# of SNPs per gene
#
# - obj.info : dataframe with fields objID, objStat, objBin and fld
# - fld      : gene scores will be plotted against this column of obj.info 
# - xlab     : x-axis label
# - ylab     : y-axis label
# - logaxis  : which of the axis should be shown in logvalues?
# - show.bins: plot dividers between bins?
#===========================================================================

PlotStatField<-function(obj.info, fld="objSNPcnt", 
                          xlab="stat", ylab="SNPs per gene", 
                          logaxis="y", show.bins=T) {
  
  # plot genes
  plot(obj.info$objStat,obj.info[[fld]],xlab=xlab,
       ylab=ylab, log=logaxis)
  
  if (show.bins) {
  # plot bins  
    bin.min<-sapply(split(seq(1:nrow(obj.info)),obj.info$objBin),
                    FUN=function(x){min(obj.info[[fld]][x])})
    bin.max<-sapply(split(seq(1:nrow(obj.info)),obj.info$objBin),
                    FUN=function(x){max(obj.info[[fld]][x])})  
    
    nbins<-length(bin.min)
    bin.borders<-bin.max[1:(nbins-1)]+0.5
    for (i in 1:length(bin.borders))
      abline(h=bin.borders[i],col="red")
  }
  
}


#===========================================================================
# Function: EnrichmentAnalysis(set.info, set.obj, obj.stat, nrand, 
#                              approx.null, seq.rnd.sampling, use.bins, test,
#                              do.pruning, minsetsize, project.txt, 
#                              do.emp.fdr, emp.fdr.path, emp.fdr.nruns, 
#                              emp.fdr.est.m0
#
# Complete pipeline to test gene sets for a significant shift 
# in the distribution of their gene scores
#
# - set.info        : dataframe with fields setID, ..
# - set.obj         : dataframe with fields setID, objID, ..
# - obj.stat        : dataframe with fields objID, objStat, ..
# - nrand           : number of randomizations
# - approx.null     : approximate null distribution with normal distribution?
# - seq.rnd.sampling: use sequential random sampling to improve speed?
# - use.bins        : use bins (sample sets with same bin distribution)?
# - test            : choose "highertail", "lowertail", or "twosided"
# - do.pruning      : apply 'pruning' to remove overlap between sets?
# - minsetsize      : during pruning only keep sets with at least this size
# - project.txt     : name of the project
# - do.emp.fdr      : calculate empirical fdr after pruning?
# - emp.fdr.path    : path to files with p-value null distributions
# - emp.fdr.nruns   : number of runs (randomizations) to calculate emp fdr
# - emp.fdr.est.m0  : estimate proportion of true nulls to calculate emp fdr?
# - qvalue.method   : to set parameter pi0.method of function qvalue
#===========================================================================

EnrichmentAnalysis<-function(set.info, set.obj, obj.stat,
                       nrand=10000, approx.null=FALSE, 
                       seq.rnd.sampling=FALSE,
                       use.bins=FALSE, test="highertail",
                       do.pruning=FALSE, minsetsize=10,
                       project.txt="DAGI", do.emp.fdr=FALSE, 
                       emp.fdr.path=".", 
                       emp.fdr.nruns=10, 
                       emp.fdr.est.m0=TRUE,
                       qvalue.method="bootstrap"){
  
  set.seed(nrand)
  
  
  set.scores.prepruning<-TestGeneSets(set.info, set.obj, obj.stat,
                 nrand=nrand, approx.null=approx.null, 
                 seq.rnd.sampling=seq.rnd.sampling,
                 use.bins=use.bins, test="highertail")
    
  pack.qvalue.missing <- !(require(qvalue))
  if (pack.qvalue.missing) {
    set.scores.prepruning$setQ <- -1
  } else {
    set.scores.prepruning$setQ <- qvalue(set.scores.prepruning$setP, 
                                         pi0.method=qvalue.method)$qvalues
  }
  
  if (do.pruning) {
    # check whether obj.stat only contains objects that are part of a set
    # if so, set keep.obj.stat=TRUE    
    keep.obj.stat<-ifelse(sum(!(obj.stat$objID %in% set.obj$objID))==0,
                          TRUE, FALSE)
    
    set.scores.postpruning <- PruneSets(set.info, set.obj, obj.stat,
                                        set.scores.prepruning, 
                                        nrand=nrand, minsetsize=minsetsize,
                                        approx.null=approx.null, test=test,
                                        keep.obj.stat=keep.obj.stat,
                                        seq.rnd.sampling=seq.rnd.sampling,
                                        use.bins=use.bins)
    
    # calculate emperical FDR    
    if (do.emp.fdr){
      set.scores.postpruning<-
        GetEmpericalFDR(set.scores=set.scores.postpruning,
                        in.path=emp.fdr.path, 
                        nruns=emp.fdr.nruns, 
                        est.m0=emp.fdr.est.m0,
                        project.txt=project.txt,
                        nrand=nrand)
    }
        
    return (list(set.scores.prepruning=set.scores.prepruning, 
                 set.scores.postpruning=set.scores.postpruning))
  } else {
    return (list(set.scores.prepruning=set.scores.prepruning))
  }
  
}


#===========================================================================
# Function: TestGeneSets(set.info, set.obj, obj.stat, nrand, approx.null,
#                        seq.rnd.sampling, use.bins, test)
# Test gene sets for a significant shift in the distribution of gene scores

# - set.info        : dataframe with fields setID, ..
# - set.obj         : dataframe with fields setID, objID, ..
# - obj.stat        : dataframe with fields objID, objStat, ..
# - nrand           : number of randomizations
# - approx.null     : approximate null distribution with normal distribution?
# - seq.rnd.sampling: use sequential random sampling to improve speed?
# - use.bins        : use bins (sample sets with same bin distribution)?
# - test            : choose "highertail", "lowertail", or "twosided"
#===========================================================================
TestGeneSets<-function(set.info, set.obj, obj.stat,
                       nrand=10000, approx.null=FALSE, 
                       seq.rnd.sampling=FALSE,
                       use.bins=FALSE, test="highertail"){
  
  M<-nrow(set.info)
  N<-nrow(obj.stat)
  
  # Create matrix setobjmat with nrows = # sets, ncols = # objs
  # each row is logical vector indicating which obj is in set
  set.obj.mtx<-CreateSetObjMtx(set.obj, set.info$setID,obj.stat$objID)
  # Calculate setsize for all sets
  set.sizes <- rowSums(set.obj.mtx)
  # Calculate sumstat for all sets
  set.stats <- as.vector(set.obj.mtx %*% obj.stat$objStat)  
    
  if (approx.null==FALSE) { 
    
    if (use.bins | seq.rnd.sampling) {      
            
      if (!("objBin" %in% names(obj.stat)))
        obj.stat$objBin<-1
      
      ix.order<-order(obj.stat$objBin)      
      obj.stat<-obj.stat[ix.order,]
      set.obj.mtx<-set.obj.mtx[,ix.order]
        
      # create null dist and count number of randomized sets with SUMSTAT 
      # as high (low) as scores
      # create matrix indicating which object is in which bin
      max.bin<-max(obj.stat$objBin)
      bin.mtx<-matrix(0,nrow=N,ncol=max.bin)    
      for (i in seq(max.bin)){
        ix<-which(obj.stat$objBin==i)                
        bin.mtx[min(ix):max(ix),i]<-1
      }
      # do matrix multipl.
      # create matrix indicating per set how many genes in each bin
      bin.dist.mtx<-set.obj.mtx %*% bin.mtx
      
      # create null distr in steps of maxt
      # only select candidates for the next run
      # which have less than min.hits hits
      
      # get some bin properties    
      #bin.ix<-rep(0,max.bin)
      bin.ix<-sapply(seq(max.bin),function(x){min(which(obj.stat$objBin==x))})
      bin.ix<-c(bin.ix,N+1)
      
      if (seq.rnd.sampling) {
        step.nrand<-min(nrand,10000)
      } else {
        step.nrand<-nrand
      }
      
      min.hits<-5000
      null.cnt<-rep(0,M)    
      null.size<-rep(0,M)
      tmp.nrand<-0
      ix.cand<-seq(M)
      
      while(tmp.nrand<nrand & length(ix.cand)>0){
                
        null.cnt[ix.cand]<-null.cnt[ix.cand] +        
          CreateNullDist_Binned(stats=obj.stat$objStat,bin.ix=bin.ix,
                                maxbin=max.bin,bindistmtx=bin.dist.mtx[ix.cand,],
                                nrand=step.nrand,scores=set.stats[ix.cand],
                                test=test)
                
        tmp.nrand<-tmp.nrand+step.nrand
        null.size[ix.cand]<-tmp.nrand
        
        ix.cand<-which(null.cnt<min.hits)
      }      
      
      # calculate p values
      setP <- (null.cnt+1)/(null.size+1)
      if (test=="twosided") setP <- 2* setP
      
    } else {
    
      # Create null distribution without using bins or 
      # sequential random sampling
      # We need one null distribution per setsize
      setsizes.unique<-sort(unique(set.sizes))
      nullmtx<-matrix(0,nrow=length(setsizes.unique),ncol=nrand)
      nullmtx<-replicate(nrand,
                          cumsum(obj.stat$objStat[sample(seq(N),
                          max(setsizes.unique))])[setsizes.unique])
      
      # Get p-values
      setsize.ix<-match(set.sizes, setsizes.unique)
      
      if (test %in% c("highertail","twosided")){
        cnt.right<-sapply(seq(M),
                    FUN = function(x) {
                      sum(nullmtx[setsize.ix[x],]>=set.stats[x])
                    })
        p.right<-(cnt.right+1)/(nrand+1)
      } 
      if (test %in% c("lowertail","twosided")){
        cnt.left<-sapply(seq(M),
                          FUN = function(x) {
                            sum(nullmtx[setsize.ix[x],]<=set.stats[x])
                          })
        p.left<-(cnt.left+1)/(nrand+1)
      }
      
      # get proper p-value
      if (test=="highertail") {
        setP<- p.right
      } else if (test=="lowertail") {
        setP<- p.left
      } else {
        setP <- pmin(2*p.right, 2*(1-p.right))
      }      
      
    } 
    
  } else {
    
    # or infer p-values from approximate null distribution
    # using mean and var objStat
    
    stat.mean<-mean(obj.stat$objStat)
    stat.var<-var(obj.stat$objStat)  
    stat.sd.corr <- sqrt(stat.var*(N-set.sizes)/(N-1))
        
    p.right<-sapply(seq(M), 
               FUN=function(x){
                 pnorm(set.stats[x],
                       mean=stat.mean*set.sizes[x],
                       sd=stat.sd.corr[x]*sqrt(set.sizes[x]),
                       lower.tail=FALSE)
               })
    p.left<-1-p.right
    
    # get proper p-value
    if (test=="highertail") {
      setP<- p.right
    } else if (test=="lowertail") {
      setP<- p.left
    } else {
      setP <- pmin(2*p.right, 2*(1-p.right))
    }
  
  } 

  
  set.scores <- data.frame(set.info$setID,rep(0,M),rep(0,M),
                           rep(0,M),rep(0,M), stringsAsFactors=F)
  names(set.scores) <- c("setID","setSize","setScore","setP", "setQ")
  set.scores$setSize<-set.sizes
  set.scores$setScore<-set.stats
  set.scores$setP<-setP
    
  set.scores$setName <- set.info$setName
  set.scores$setID.orig <- set.info$setID.orig
  set.scores<-set.scores[order(set.scores$setP),]
  
  return(set.scores)
  
}

#===========================================================================
# Function: PruneSets(obj.stat, set.obj, set.info, set.scores, 
#                     nrand, minsetsize,approx.null, 
#                     test, keep.obj.stat)
#
# Iteratively remove genes in top scoring set from other sets
# Inspired by topGO method.
# Reference:
# Alexa, A., Rahnenfuhrer, J., and Lengauer, T. (2006). 
# Improved scoring of functional groups from gene expression data by 
# decorrelating GO graph structure. Bioinformatics 22, 1600–1607.
#
# - set.info        : dataframe with fields setID, ..
# - set.obj         : dataframe with fields setID, objID, ..
# - obj.stat        : dataframe with fields objID, objStat, ..
# - set.scores      : dataframe with test scores
# - nrand           : number of randomizations
# - minsetsize      : during pruning only keep sets with at least this size
# - approx.null     : approximate null distribution with normal distribution?
# - seq.rnd.sampling: use sequential random sampling to improve speed?
# - use.bins        : use bins (sample sets with same bin distribution)?
# - test            : choose "highertail", "lowertail", or "twosided"
# - keep.obj.stat   : if TRUE, don't remove objects from obj.stat
#                     during pruning (recommended when obj.stat
#                     only contains objects that are part of 
#                     at least one geneset)
#===========================================================================
PruneSets<-function(set.info, set.obj, obj.stat, set.scores, 
                    nrand=10000, minsetsize=10,
                    approx.null=FALSE, seq.rnd.sampling=TRUE,
                    use.bins=FALSE, test="higher.tail",
                    keep.obj.stat=FALSE) {  

  # start with new dataframe
  set.scores.nw<-data.frame()
    
  while (nrow(set.scores)>1 & nrow(set.obj)>0) {
    # add first row of set.scores to new set.scores
    # and remove it from old set.scores
    set.scores.nw<-rbind(set.scores.nw,set.scores[1,]) 
    # get genes from first row    
    setID<-set.scores$setID[1]        
    set.scores<-set.scores[-1,]    
    
    # remove genes in setID from other sets
    genes<-set.obj$objID[set.obj$setID==setID]    
    set.obj<-set.obj[!(set.obj$objID %in% genes),]
    # remove setID from set.info
    set.info<-set.info[set.info$setID!=setID,]
    
    if (!keep.obj.stat){
      # remove genes in setID from obj.stat
      obj.stat<-obj.stat[!(obj.stat$objID %in% genes),]            
    }
    
    # remove sets size<minsetsize from set.obj and set.info
    t<-table(set.obj$setID)
    small.sets<-unique(names(t[t<minsetsize]))
    set.obj<-set.obj[!(set.obj$setID %in% small.sets),]
    set.info<-set.info[!(set.info$setID %in% small.sets),]
    # remove sets from set.info without obj's in set.obj
    set.info<-set.info[set.info$setID %in% set.obj$setID,]
    
    if (nrow(set.obj)>0) {
      # test the remaining sets
      set.scores<-TestGeneSets(set.info, set.obj, obj.stat,
                               nrand=nrand, approx.null=approx.null,
                               test=test, seq.rnd.sampling=seq.rnd.sampling,
                               use.bins=use.bins)       
    }
  }    
  
  set.scores.nw<-rbind(set.scores.nw,set.scores)
  
  # order on p-value
  set.scores.nw<-set.scores.nw[order(set.scores.nw$setP),] 
  set.scores.nw$setQ <- -1  
  
  return(set.scores.nw)
}

#===========================================================================
# Function: TestShuffledSets(set.info, set.obj, obj.stat, 
#                            nrand, minsetsize, approx.null, test, 
#                            seq.rnd.sampling, use.bins, keep.obj.stat,
#                            out.path, runnr, project.txt)
#
# Shuffle stat in obj.stat and redo enrichment test
# The network structure remains unchanged
#
# - set.info        : dataframe with fields setID, ..
# - set.obj         : dataframe with fields setID, objID, ..
# - obj.stat        : dataframe with fields objID, objStat, ..
# - nrand           : number of randomizations
# - minsetsize      : during pruning only keep sets with at least this size
# - approx.null     : approximate null distribution with normal distribution?
# - test            : choose "highertail", "lowertail", or "twosided"
# - seq.rnd.sampling: use sequential random sampling to improve speed?
# - use.bins        : use bins (sample sets with same bin distribution)?
# - out.path        : path to store output files
# - runnr           : run number
# - project.txt     : name of project
#===========================================================================
TestShuffledSets<-function(set.info, set.obj, obj.stat, 
                           nrand=10000, minsetsize=10, 
                           approx.null=FALSE, test="higher.tail",
                           seq.rnd.sampling=TRUE,
                           use.bins=FALSE,
                           out.path=".", runnr=1,
                           project.txt="DAGI"){
  
  # set seed to reproduce results
  set.seed(runnr)  
  
  # shuffle stats
  N<-nrow(obj.stat)
  obj.stat.shuffle<-obj.stat  
  ix.shuffle<-sample(1:N,N,replace=FALSE)
  obj.stat.shuffle$objStat<-obj.stat$objStat[ix.shuffle]
  
  # check whether obj.stat only contains objects that are part of a set
  # if so, set keep.obj.stat=TRUE   
  keep.obj.stat<-ifelse(sum(!(obj.stat$objID %in% set.obj$objID))==0,
                        TRUE, FALSE)
    
  # test sets
  set.scores.shuffle<-TestGeneSets(set.info, set.obj, 
                                  obj.stat.shuffle,
                                  nrand=nrand, 
                                  approx.null=approx.null,
                                  seq.rnd.sampling=seq.rnd.sampling,
                                  use.bins=use.bins,
                                  test=test)
  set.scores.shuffle.postpruning<-PruneSets(set.info, set.obj, 
                                  obj.stat.shuffle, set.scores.shuffle, 
                                  nrand=nrand, 
                                  minsetsize=minsetsize,
                                  approx.null=approx.null, 
                                  seq.rnd.sampling=seq.rnd.sampling,
                                  use.bins=use.bins,
                                  test=test,
                                  keep.obj.stat=keep.obj.stat)

  p.shuffle.postpruning<-set.scores.shuffle.postpruning$setP
  
  f<-paste(project.txt,"_shuf",runnr,"_",formatC(nrand,format="d"),
           ".RData",sep="")
  save(p.shuffle.postpruning,file=file.path(out.path,f))
  
}

#===========================================================================
# Function: GetEmpericalFDR(set.scores, in.path, nruns, est.m0, project.txt)
# Calculate the fdr emperically from a null distribution of p-values
#
# - set.scores : dataframe with set scores
# - in.path    : path to files with p-value distributions
# - nruns      : number of runs (files are numbered from 1 to nruns)
# - est.m0     : if TRUE the proportion of true H0's will be estimated using
#                a histogram based method
# - project.txt: name of project
#===========================================================================
GetEmpericalFDR<-function(set.scores, in.path, nruns=3, est.m0=T,
                          project.txt="DAGI", nrand) {
  
  p.shuffle.postpruning.ls<-list()
  
  # remove p-values lower than 0
  ix.p<-which(set.scores$setP>=0)
  p.lim<-set.scores$setP[ix.p]
  
  # Get p-values from tests with shuffled objStats 
  for (i in seq(nruns)){
    f<-paste(project.txt,"_shuf",i,"_",formatC(nrand,format="d"),
             ".RData",sep="")
    load(file=file.path(in.path,f))    
    p.shuffle.postpruning.ls[[i]]<-p.shuffle.postpruning      
  }     
  
  # remove p-values lower than 0 from shuffle list
  for (i in seq(nruns)){
    ls.tmp<-p.shuffle.postpruning.ls[[i]]
    p.shuffle.postpruning.ls[[i]]<-ls.tmp[ls.tmp>=0]
  }
    
  # estimate pi.0, the proportion of true H0's
  nsc<-nrow(set.scores)
  if (est.m0==T){
    p.perm.pooled<-numeric()
    for (p_ in seq(nruns)){
      p.perm.pooled<-c(p.perm.pooled,
                       p.shuffle.postpruning.ls[[p_]])
    }
    
    # take mean per run:
    M0<-M0_from_hist(bins=20,p=set.scores$setP[ix.p],e=0.0001,
                     p.rand=p.shuffle.postpruning.ls)
    pi.0<-min(1,M0/nsc)
    
  } else {
    pi.0<-1
  }
  
  # get ratio #p-values in real tests/#p-values shuffled tests
  fdr.est1<-numeric()
  for (p in 1:length(p.lim)){
    # calculate running mean of N0
    n <- 0
    mean.N0 <- 0
    for (run in 1:nruns){
      nsh<-length(p.shuffle.postpruning.ls[[run]])
      N0<-sum(p.shuffle.postpruning.ls[[run]]<=p.lim[p])/nsh
      n <- n + 1
      delta <- N0 - mean.N0
      mean.N0 <- mean.N0 + delta/n
    }
    N.t<-sum(set.scores$setP<=p.lim[p])
    
    fdr.est1[p]<-pi.0*nsc*mean.N0/N.t
  }
  
  # get the smallest ratio per p-value
  fdr.est2<-numeric()
  p.lim.n<-length(p.lim)
  for (i in 1:length(p.lim)){
    fdr.est2[i]<-min(fdr.est1[i:p.lim.n])
  }
  
  set.scores$setQ <- -1
  for (i in ix.p){
    set.scores$setQ[i]<-
      fdr.est2[min(which(p.lim>=set.scores$setP[i]))]
  }
  
  return(set.scores=set.scores)
  
}


#===========================================================================
# Function: CreateNullDist_Binned(stats, bin.ix, maxbin, bindistmtx, 
#                                 scores, nrand, test)
# Create null distribution of random gene sets that have the same bin 
# distribution as tested sets
#
# - stats     : vector of object scores (ordered on bin)
# - bin.ix    : vector of indices reflecting the start and end of each bin
#               example: bin.ix=c(1,4,6) means: bin 1 contains objects 1,2,3,
#               bin 2 contains objects 4,5
# - maxbin    : max bin number
# - bindistmtx: matrix indicating per set how many genes in each bin 
# - scores    : the sumstat scores of the gene sets
# - nrand     : number of randomizations
# - test      : choose "highertail", "lowertail", or "twosided"
#===========================================================================
CreateNullDist_Binned<-function(stats, bin.ix, maxbin, bindistmtx, 
                                scores, nrand, test="highertail"){
  
  # if no binning create null distribution differently
  if (length(bin.ix)==2){
    nms<-names(bindistmtx)
    dim(bindistmtx)<-c(length(scores),1)
    dimnames(bindistmtx)<-list(nms,"1")      
  }
  
  
  if (length(scores)==1){
    # if only one set is tested (nrows bindistmtx=1)
    # get null distribution differently    
    
    # requires package Matrix
    setobjmtx_null<-CreateSetObjMtx_Null(n_objs=length(stats), 
                                         setsize=sum(bindistmtx), 
                                         nrand=nrand, bin.ix=bin.ix,
                                         bins_ngenes=bindistmtx)
    testnull<-as.vector(setobjmtx_null %*% stats)
    
    return(sum(testnull>=scores))
    
  } else {  
    
    bindistmtx_ix<-bindistmtx  
    
    #create empty list, of length maxbin
    binlist<-vector("list",maxbin) 
    
    for (i in 1:maxbin){  
      bindist<-sort(unique(bindistmtx[,i]),decreasing=T)
      bindistmtx_ix[,i]<-factor(bindistmtx[,i],levels=bindist)
      bdl<-length(bindist)
      if (bindist[bdl]!=0) {
        bindist<-c(bindist,0) 
        bdl<-bdl+1
      }
      bindist<-c(bindist,-1)  
      bdl<- bdl+1    
      
      # create mtx to put the results in
      nullmtx<-matrix(0,nrow=bdl,ncol=nrand)
      
      cnt<-bin.ix[i]:(bin.ix[i+1]-1)
      max.ngenes<-bindist[1]
      bindist_tmp<-bindist[1:(bdl-1)]
      if (bdl>2){    
        nullmtx[1:(bdl-2),]<-replicate(nrand,
                              cumsum(stats[sample(cnt,max.ngenes)])[bindist_tmp])
      }        
      binlist[[i]]<-nullmtx       
    }
    
    
    setIDs<-unlist(dimnames(bindistmtx)[1])  
    N<-length(setIDs)  
    
    testnull<-matrix(0,nrow=N,ncol=nrand)
    testnullmtx<-matrix(0,nrow=maxbin,ncol=nrand) 
    
    testnull<-t(sapply(1:N, function(j){
      bindist.set<-bindistmtx_ix[j,]
      testnullmtx<-sapply(1:maxbin,                          
                          function(b){binlist[[b]][bindist.set[b],]})
      return(rowSums(testnullmtx))
    }))
    
    if (test=="highertail"){
      return(rowSums(testnull>=scores))
    } else if (test=="lowertail"){
      return(rowSums(testnull<=scores))
    } else {
      return (pmin(rowSums(testnull<=scores), rowSums(testnull>=scores)))
    }
  }
} 

#===========================================================================
# HELPER FUNCTIONS
#===========================================================================


#===========================================================================
# Function: CreateSetObjMtx(set.obj, objIDs, setIDs)
# create matrix setobjmtx with
# nrows = # sets, ncols = # objs
# each row is logical vector indicating which obj is in set
#
# -set.obj  : dataframe with fields setID, objID
# -setIDs   : (optional) dataframe with fields setID, ...
# -objIDs   : (optional) dataframe with fields objID, ...
#===========================================================================

CreateSetObjMtx<-function(set.obj, setIDs, objIDs){
    
  if (missing(setIDs))
    setIDs<-unique(set.obj$setID)
  k <- length(setIDs)
  
  if (missing(objIDs))
    objIDs<-unique(set.obj$objID)
  l<-length(objIDs)  

  
  my_func<-function(setID){
    obj_i <- set.obj$objID[set.obj$setID==setID]
    return(objIDs %in% obj_i)
  }
  setobjmtx<-matrix(unlist(lapply(setIDs,my_func)),
                    nrow=k,ncol=l,byrow=T,
                    dimnames=list(setIDs,objIDs))
  
  return(setobjmtx)
}

#===========================================================================
# Function: CreateSetObjMtx_Null(n_objs,setsize,nrand,bin.ix,bins_ngenes)
# Create null distribution setobjmtx matrix when testing only one gene set
# nrows = # randomizations, ncols =  # objs
# Each row is logical vector indicating which obj is in random set
#
# n_objs     : total number of objects
# setsize    : number of genes in set
# nrand      : number of randomizations
# bin.ix     : vector of indices reflecting the start and end of each bin
#              example: bin.ix=c(1,4,6) means: bin 1 contains objects 1,2,3,
#              bin 2 contains objects 4,5
# bins_ngenes: vector indicating for the set how many genes are in each bin 
#===========================================================================

CreateSetObjMtx_Null<-function(n_objs, setsize, nrand, bin.ix, bins_ngenes){  
  require(Matrix)  
  
  if (missing(bin.ix)){
    spmtx<-spMatrix(nrand, n_objs, i=rep(1:nrand,each=setsize), 
                    j=replicate(nrand,sample(1:n_objs,setsize)), 
                    x=rep(1,nrand*setsize))
  } else {   
    n_bins<-length(bin.ix)-1
    spmtx<-spMatrix(nrand, n_objs, 
                    i=rep(1:nrand,each=setsize), 
                    j=replicate(nrand,unlist(lapply(1:n_bins,
                        function(x){sample(bin.ix[x]:(bin.ix[x+1]-1),
                                           bins_ngenes[x])}))), 
                    x=rep(1,nrand*setsize)
    )
  }  
  
  return(spmtx)
}

#===========================================================================
# Function: AssignBins(obj.info, fld, nbins, min.bin.size, max.bin.size)
# Assign bins to objects based on fld 
#
# obj.info    : dataframe with columns: geneID and fld
# fld         : column of obj.info which is used to bin the genes 
# nbins       : number of bins; if -1, bins will be created with size 
#               between min.bin.size and max.bin.size, while trying to 
#               keep genes with the same fld value in the same bin 
# min.bin.size: the minimum number of genes in a bin
# max.bin.size: the maximum number of genes in a bin
#===========================================================================

AssignBins<- function (obj.info, fld="SNPcount", nbins=-1,
                        min.bin.size=1000, max.bin.size=1800) 
{
  # order on fld 
  obj.info <- obj.info[order(obj.info[[fld]]),]
  
  if (nbins>0){
    # method 1: assign genes to equal sized bins
    N<-nrow(obj.info)
    obj.info$objBin <- as.integer((1:N/(N/nbins)))+1 
    obj.info$objBin[obj.info$objBin==nbins+1] <- nbins
  } else {
    # method 2: for genes with low value of fld, keep 
    # same together
    # define under and upper limit and keep adding groups until limit
    # is reached
    # note that bins larger than max.bin.size can occurr, when there are
    # more than max.bin.size gene with the same fld value
    obj.info$objBin<-0
    bin <- 1
    tot.size <- 0
    t<-table(obj.info[[fld]])
    for (i in 1:length(t)){
      tot.size <- tot.size+t[i]
      # assign bin to fld
      obj.info$objBin[obj.info[[fld]]==names(t[i])]<-bin
      # check if bin is full
      if (tot.size >= min.bin.size){
        bin <- bin+1
        tot.size <- 0
      }
    }
    if (sum(obj.info$objBin %in% c(bin,bin-1))<max.bin.size) 
      obj.info$objBin[obj.info$objBin==bin] <- bin-1
  }
  
  return(obj.info)
}

#===========================================================================
# Function: RescaleBins(obj.info, fld)
# Rescale the values in field fld per bin 
#
# obj.info    : dataframe with columns: geneID and fld
# fld         : column of obj.info which contains values to be scaled 
#===========================================================================

RescaleBins <- function(obj.info, fld="objStat") 
{

  obj.info[[paste(fld,".orig",sep="")]]<-obj.info[[fld]]
  maxbin<-max(obj.info$objBin)
  for (b in 1:maxbin){
    ix <- obj.info$objBin==b
    # compute modified z
    obj.info[[fld]][ix]<-modz(obj.info[[fld]][ix])
  }
  
  return(obj.info)
}

#===========================================================================
# Function: modz(x)
# Returns modified z-score of vector x
#
# Reference:
# Iglewicz, B. and Hoaglin, D. (1993). 
# How to detect and handle outliers.
# Milwaukee (WI): ASQC Quality Press.
#===========================================================================

modz<- function (x) 
{
  y<-x-median(x)
  return(0.6745*y/median(abs(y)))
}


#===========================================================================
# Function: M0_from_hist(bins, p, e, p.rand)
# Estimate the number of true H0's from histogram
#
# Reference:
# Nettleton, D., Hwang, J., Caldo, R., and Wise, R. (2006). 
# Estimating the number of true null hypotheses from a histogram 
# of p values. Journal of Agricultural, Biological, and Environmental 
# Statistics 11, 337–356. 
#
# bins    : number of histogram bins
# p       : vector with observed p-values
# e       : convergence parameter
# p.rand  : vector with p-values from randomizations, used to infer expected 
#           frequencies per bin
# plothist: plot a histogram showing the frequency of observed p-values after
#           SUMSTAT testing with pruning, and a blue line depicting the 
#           expected frequencies of true null hypotheses 
#===========================================================================
M0_from_hist <- function(bins=20, p, e=0.0001, p.rand, plothist=F) {

  # calculate upper and lower bounds bins
  bins.tmp<-seq(from=0,to=1,length.out=bins+1)
  # calculate M, total number of observations
  M<-length(p)
  # set M0 to M
  M0<-M  
  
  # if p.rand is list, calculate means per bin
  if (is.list(p.rand)){    
    N<-length(p.rand)
    fmtx<-matrix(rep(0,N*bins),nrow=N)
    for (n in 1:N){
      h.rand<-hist(p.rand[[n]],breaks=bins.tmp,plot=F)
      fmtx[n,]<-h.rand$counts/length(p.rand[[n]])
    }
    f<-apply(fmtx,MARGIN=2,mean)
  } else {  
    h.rand<-hist(p.rand,breaks=bins.tmp,plot=F)
    f<-h.rand$counts/length(p.rand)    
  }
  
  # count number of p-values per bin
  h<-hist(p,breaks=bins.tmp,plot=F)
  p.bin.obs<-h$counts

  ploti=1  
  if (plothist){
    ploti=1
    graphics.off()
    windows(width=10,height=4)
    op<-par(mfrow=c(1,3))    
  }
  repeat{
    M0.old<-M0
    # calculate the expected number of p-values per bin
    p.bin.exp<-f*M0
    # find the first bin where exp > obs
    bin.lim<-min(which(p.bin.exp > p.bin.obs))-1
    M0<-M-(sum(p.bin.obs[1:bin.lim]-p.bin.exp[1:bin.lim]))
    if (abs(M0-M0.old) < e) break
    if (plothist & ploti<=3){
      plot(h,col=colors()[425],main=paste("iteration ",ploti,sep=""))  
      lines(bins.tmp,c(p.bin.exp,p.bin.exp[length(p.bin.exp)]),type="s",
            col="blue",lwd=2)      
      ploti=ploti+1
    }
    
  }  
  
  if (plothist) par(op)  
  
  return(M0)
}
  