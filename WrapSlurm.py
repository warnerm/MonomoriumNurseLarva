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
	        opts, args = getopt.getopt(argv,"h:s:b:n:",["sfile=","bfile=","nfile="])
	except getopt.GetoptError:
	        print 'WrapSlurm.py -f <fpkm> -s <samp> -b <totalboots> -n <totGene>'
	        sys.exit(2)
	for opt, arg in opts:
	        if opt in ("-s","--sfile"):
	                samp = arg   
	        elif opt in ("-b","--bfile"):
	                boots = arg   
	        elif opt in ("-n","--nfile"):
	                nGene = arg                    
	return samp,boots,nGene

def main(argv):
	samp,boots,nGene = InOut(argv)
	name = samp + '_' + boots + '_' + nGene
	for x in range(boots/500):
		call(["sbatch","slurmParallel.sh","name"])

if __name__ == "__main__":
	main(sys.argv[1:])



