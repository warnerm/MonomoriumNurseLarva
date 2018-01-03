#!/usr/local/bin/python2.7
import subprocess,sys, getopt
from subprocess import call
import os


def InOut(argv):
	samp = ''
	boots = ''
	nGene = ''
	try:
	        opts, args = getopt.getopt(argv,"h:s:b:n:",["sfile=","bfile=","nfile="])
	except getopt.GetoptError:
	        print 'WrapSlurm.py -s <samp> -b <totalboots> -n <totGene>'
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
	os.environ['name']=name
	for x in range(int(boots)/1000):
		os.environ['run'] = str(x)
		call(["sbatch","slurmParallel.sh"])

if __name__ == "__main__":
	main(sys.argv[1:])



