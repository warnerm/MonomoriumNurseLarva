#!/bin/bash
#SBATCH -p compute # partition (queue)
#SBATCH --export=ALL
Rscript RandomNetworks_parallel.R $name