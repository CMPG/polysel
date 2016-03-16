## Examples of how to run the sumstat enrichment test on a cluster

Before using these commands, find out what job submission system is installed on your organization's cluster and what are the guidelines for using your cluster. Here we list examples when using a SUN Grid Engine (SGE).  

Run these commands from the root of your project folder (*polysel*).

###SUN Grid Engine job submission system###

**STEP 1: run polysel_shufflesets_cluster.sh script with:**

Change to `for i in {1..3}` if you just want to try out if everything is working correctly.

	  projectname="primates_homininae"
	  for i in {1..300}
	  do
	  qsub -l h_cpu=5:00:00 -l h_vmem=4G -o ./log/shuffle${i}_${projectname}.out \
	  -e ./log/shuffle${i}_${projectname}.err -cwd -N shuffle${i}_${projectname} \
	  ./script/polysel_shufflesets_cluster.sh $i $projectname
	  done


**or (USING SCRATCH): run polysel_shufflesets_cluster_scratch.sh script with:**

	  projectname="primates_homininae"
	  for i in {1..300}
	  do
	  qsub -l h_cpu=5:00:00 -l h_vmem=4G -o ./log/shuffle${i}_${projectname}.out \
	  -e ./log/shuffle${i}_${projectname}.err -cwd -N shuffle${i}_${projectname} \
	  ./script/polysel_shufflesets_cluster_scratch.sh $i $projectname
	  done


**STEP 2: run polysel_cluster.sh script with:**

	  projectname="primates_homininae"
	  qsub -l h_cpu=5:00:00 -l h_vmem=4G -o ./log/${projectname}.out \
	  -e ./log/${projectname}.err -cwd -N ${projectname} \
	  ./script/polysel_cluster.sh $projectname

**or: (USING SCRATCH): run polysel_cluster_scratch script with:**

	  projectname="primates_homininae"
	  qsub -l h_cpu=5:00:00 -l h_vmem=4G -o ./log/${projectname}.out \
	  -e ./log/${projectname}.err -cwd -N ${projectname} \
	  ./script/polysel_cluster_scratch.sh $projectname