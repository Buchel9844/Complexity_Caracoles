#!/bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Thu Apr 13 2023 16:11:40 GMT+1100 (Australian Eastern Daylight Time)

# Partition for the job:
#SBATCH --partition=physical

# Multithreaded (SMP) job: must run on one node 
#SBATCH --nodes=1

# The name of the job:
#SBATCH --job-name="Complexity_code"

# The project ID which this job should run under:
#SBATCH --account="punim1670"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# The amount of memory in megabytes per process in the job:
#SBATCH --mem-per-cpu=40960

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=3-1:0:00

#specify output file
#SBATCH --output=complexmodelcode.out

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

# The modules to load:
module load r/4.2.0

# The job command(s): args: runprelimfit runfinalfit focal year complexity.level
Rscript --vanilla /home/lbuche/Eco_Bayesian/Complexity_caracoles/code/Caracoles_ModelFit.R 1 1 1

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s