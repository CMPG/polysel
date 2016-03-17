## Detection of polygenic selection in gene sets (polysel) ##

Welcome to the *polysel* GITHUB repository. Here you can find code, examples, scripts and data sets to run a gene set enrichment test to detect pathways enriched with signals of polygenic selection.

We use the following folder structure:

* **polysel**
    + **R** - *R files with all functions needed to run the pipeline*
    + **data** - *Input text files with gene and gene set information, but also the R dataobject files that are created during data preparation. Each project has its own subfolder*
    + **empfdr** - *Output files from the creation of the p-value null distribution (see below)*
    + **log** - *Location of error and output log files, mainly used on the cluster*
    + **script** - *Script files to run the pipeline on a cluster*
    + **results** - *Results of the gene set enrichment test, namely a list of gene sets with their scores and p and q values*

When you are running the pipeline with your own data sets or one of the example sets, we advice to use the same structure, both on your PC and on a computer cluster.

In the root folder of this project, you will find examples explaining step by step how to run the enrichment pipeline. We strongly encourage you to try one or more of them before running the pipeline with your own data.

####Examples

**Primates** project: *Detection of pathways affected by positive selection in primate lineages ancestral to humans*

*	[example_primates_homininae.html](https://github.com/CMPG/polysel/blob/master/example_primates_homininae.html "example_primates_homininae.html"); click [here](http://htmlpreview.github.com/?https://github.com/CMPG/polysel/blob/master/example_primates_homininae.html "example_primates_homininae") for a html preview

*	[example_primates_hominidae.html](https://github.com/CMPG/polysel/blob/master/example_primates_hominidae.html "example_primates_hominidae.html"); click [here](http://htmlpreview.github.com/?https://github.com/CMPG/polysel/blob/master/example_primates_hominidae.html "example_primates_hominidae") for a html preview


####References

  Daub, J.T. et al. (in preparation). **Detection of pathways affected by positive selection in primate lineages ancestral to humans**  

  Daub, J.T. et al. (2013). **Evidence for Polygenic Adaptation to Pathogens in the Human Genome.** Mol Biol Evol 30, 1544–1558.
