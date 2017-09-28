#!/bin/bash#
#SBATCH -p compute* # partition (queue)
#SBATCH -N 30 # number of nodes
#SBATCH -n 120 # number of cores
#SBATCH -t 0-24:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=warnerm@sas.upenn.edu # send-to addressfor i in {1..100000}; do
Rscript collect.parallel.R