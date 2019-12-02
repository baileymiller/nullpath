#!/bin/bash
# The name of the job, can be anything, simply used when displaying the list of running jobs
#$ -N my_job
# Combining output/error messages into one file
#$ -j y
# Set memory request:
#$ -l vf=1G
# Set walltime request:
#$ -l h_rt=HH:MM:SS
# One needs to tell the queue system to use the current directory as the working directory
# Or else the script may fail as it will execute in your top level home directory /home/username
#$ -cwd
# then you tell it retain all environment variables (as the default is to scrub your environment)
#$ -V
# Now comes the command to be executed
qsub -N $0 -j Y -l vf=1G -l h_rt=01:00:00 -pe smp $1 -l ironfs ./run_job.sh $2 $3
exit 0
