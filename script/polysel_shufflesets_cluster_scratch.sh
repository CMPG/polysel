#!/bin/bash

#-------------------------------------------------------------------------------------------
# project specific parameters
#-------------------------------------------------------------------------------------------

# POSITIVE SELECTION IN PRIMATES PROJECT
# Daub et al. (in preparation)
# Detection of pathways affected by positive selection in primate lineages ancestral to humans
# Total nr of runs: 300
if [ ${projectname:0:8} = "primates" ]; then
	branch=${projectname:9}
	projectdir="primates/$branch"
	nrand=1000000
	approxnull="FALSE"
	seqrndsampling="TRUE"
	usebins="FALSE"
fi

# LOCAL ADAPTATION IN HUMAN POPULATIONS PROJECT
# Daub et al. (2013)
# Evidence for Polygenic Adaptation to Pathogens in the Human Genome. 
# Mol Biol Evol 30, 1544–1558.
# Total nr of runs: 200
if [ $projectname = "humanpops" ]; then
	projectdir="humanpops"
	nrand=1
	approxnull="TRUE"
	seqrndsampling="FALSE"
	usebins="FALSE"
fi

# ADAPTATION TO HIGH ALTITUDE PROJECT
# Foll et al. (2014)
# Widespread Signals of Convergent Adaptation to High Altitude in Asia and America. 
# The American Journal of Human Genetics 95, 394–407.
# Total nr of runs: 200
if [ $projectname = "altitude" ]; then
	projectdir="altitude"
	nrand=500000
	approxnull="FALSE"
	seqrndsampling="TRUE"
	usebins="TRUE"
fi

# ADAPTATION TO FORESTS PROJECT
# Amorim et al. (2015)
# Detection of Convergent Genome-Wide Signals of Adaptation to Tropical Forests in Humans. 
# PLoS One 10.
# Total nr of runs: 300
if [ ${projectname:0:6} = "forest" ]; then
	popset=${projectname:7}
	projectdir="forest/$branch"
	nrand=500000
	approxnull="FALSE"
	seqrndsampling="TRUE"
	usebins="TRUE"
fi

#-------------------------------------------------------------------------------------------
# general test parameters
#-------------------------------------------------------------------------------------------
minsetsize=10
test="highertail"

#-------------------------------------------------------------------------------------------
# cluster parameters
#-------------------------------------------------------------------------------------------
datapath="./data"
codepath="./R"
empfdrpath="./empfdr"

#-------------------------------------------------------------------------------------------
# Create Scratch Folder
#-------------------------------------------------------------------------------------------
directory=`pwd`
scratchFolder=${TMPDIR}
cd $scratchFolder

#-------------------------------------------------------------------------------------------
# Copy Files into scratchFolder
#-------------------------------------------------------------------------------------------

cd $directory
cp ${codepath}/*.R $scratchFolder/.
cp ${datapath}/${projectdir}/*.RData $scratchFolder/.

#-------------------------------------------------------------------------------------------
# Lauch R
#-------------------------------------------------------------------------------------------

cd $scratchFolder

chmod +x $main_script
cat $main_script | R --vanilla --args $runnr $nrand $minsetsize $approxnull $test \
$seqrndsampling $usebins ${projectname} . . .

cd $directory

#-------------------------------------------------------------------------------------------
# Copy Results to Master
#-------------------------------------------------------------------------------------------

cp $scratchFolder/${projectname}_shuf*_${nrand}.RData $resultspath

