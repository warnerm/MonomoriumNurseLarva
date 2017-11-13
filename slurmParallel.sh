#!/bin/bash
#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
#SBATCH -t 0-10:00
Rscript RandomNetworks_parallel.R $name $run