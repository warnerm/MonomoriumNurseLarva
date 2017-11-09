#!/usr/local/bin/python2.7
from subprocess import call
import subprocess
from joblib import Parallel, delayed

def InOut(argv):
	fpkm=''
	samp = ''
	boots = ''
	totGene = ''
	try:
	        opts, args = getopt.getopt(argv,"h:f:s:b:n:",["ffile=","sfile=","bfile=","nfile="])
	except getopt.GetoptError:
	        print 'OGGalign.py -t <threads> -s <seqfile>'
	        sys.exit(2)
	for opt, arg in opts:
	        if opt in ("-t","--tfile"):
	                num_cores = arg 
	        elif opt in ("-s","--sfile"):
	                seqfile = arg                    
	return num_cores,seqfile

for x in range(1000):
	call(["sbatch","src/slurmParallel.sh"])