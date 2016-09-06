#!/bin/bash
#
# A script to submit SGE_SparseDemand.exe
# Run using 'qsub SGE_SparseDemand.sh

#
# Specify the Working Directory
#$ -wd /data/uctpln0/IO/Lewbel/code/fortran/source
 
# Specify an output file
#$ -o /data/uctpln0/IO/Lewbel/code/fortran/output/a1.out

# Specify whether or not the standard error stream of the job i
# is merged into the standard output stream
# y = yes
# n = no
# $ -j y

# Specify a queue name to run the job on
#$ -q batch.q

# Specify the shell to use
#$ -S /bin/bash

# Specify a job name for the job
#$ -N sparse_a1

# Specify where e-mail notifications are sent
#'$ -M l.nesheim@ucl.ac.uk

# E-mail notifications sent on events
#
# b = Mail is sent at the beginning of the job
# e = Mail is sent at the end of the job
# a = Mail is sent when the job is aborted or rescheduled
# s = Mail is sent when the job is suspended
# n = No mail is sent
#'$ -m n

# Launch the job in a queue meeting the given resource request list
# -l resource=value
#$ -l mem_total=2G

# dat-time to run, format [[CC]yy]MMDDhhmm[.SS]
#'$ -a ${1}

# Specify whether resources are essential (hard) or not (soft)
# DELETE AS APPROPRIATE
#$ -soft 

# Specify job priority
#$ -p 0

# Specify a parallel environment
# In order to use this option please remove the ' between 
# the # and $ on the next line
#$ -pe openmpi 120   
# Change the value above to increase the number of slots.

. /etc/profile.d/modules.sh
module load sge/2011.11 intel/composer/64/2015.1.133 nag/f77/fll6i25dc openmpi/intel/64/1.4.5 
echo "Loaded modules"
module list

cd /data/uctpln0/IO/Lewbel/code/fortran/source
date
mpirun -np $NSLOTS SparseDemand_mpi.exe ../inputs/A1.prop
date
