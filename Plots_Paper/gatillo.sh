#!/bin/bash
# The name of the job, can be whatever makes sense to you
#$ -N Test

# The job should be placed into the queue 'all.q'.
#$ -q all.q

# Merge output and standard output
#$ -j y

# The batchsystem should use the current directory as working directory.
# Both files will be placed in the current
# directory. The batchsystem assumes to find the executable in this directory.
#$ -cwd

# Number of tasks
#$ -t 1-1

# request Bourne shell as shell for job.
#$ -S /bin/bash

/opt/python/bin/python2.7 Intra-Inter_Distances.py 
