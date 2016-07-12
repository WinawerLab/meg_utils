#!/bin/bash
#PBS -l nodes=1:ppn=1
#PBS -l walltime=25:00:00
#PBS -l mem=120GB
#PBS -N s_GAMMA_pipeline
#PBS -M ek99@nyu.edu
#PBS -j oe

session_num=${PBS_ARRAYID}
module load matlab/2014a

cd /scratch/ek99/meg_utils/experiments/gamma_meg

# If the files you are running are not in the same folder as this script,
# you can insert "addpath(genpath('/PATH/TO/FILES/'));" before the command
# you want to run.
# matlab -nodisplay -r "addpath(genpath('/scratch/ek99/meg_utils')); HPC_Gamma($session_num,100); exit()"
matlab -nodisplay -r "addpath(genpath('/scratch/ek99/meg_utils')); s_GAMMA_pipeline($session_num); exit()"

exit