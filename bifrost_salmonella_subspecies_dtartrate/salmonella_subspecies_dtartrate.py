#!/usr/bin/env python3

import argparse
import subprocess as sp

component_dir= "/bifrost/components/bifrost_salmonella_subspecies_dtartrate"
DTARTRATEDB = component_dir + "/ressources/d-Tartrate/srst2_d_Tartrate.fasta"

def dtartrate(read1, read2):
	bwa = sp.Popen(
		["bwa", "mem", DTARTRATEDB, read1, read2], 
		stdout = sp.PIPE, 
		stderr = sp.PIPE).communicate()
	elprep = sp.Popen(
		["elprep","filter", "/dev/stdin", "/dev/stdout", "--filter-unmapped-reads", "--sorting-order", "coordinate", "--nr-of-threads", "1"],
		stdin=sp.PIPE, 
		stdout = sp.PIPE).communicate(input=bwa[0])
	getpos = sp.Popen(
		[component_dir + "/bifrost_salmonella_subspecies_dtartrate/getposfromsam.py", "/dev/stdin", "--pos", "10"], 
		stdin = sp.PIPE, 
		stdout = sp.PIPE).communicate(input=elprep[0])
	return [getpos[0].decode().replace('\n','')]

def subspecies(ST):
	try:
		int(ST)
		res = sp.run([component_dir + "/bifrost_salmonella_subspecies_dtartrate/NearestSubspecies.py","--st", ST], stdout=sp.PIPE)
		resultoutput = [res.stdout.decode().strip('\n')]
	except ValueError:
		""" No known st """
		resultoutput = ['NA']
	return resultoutput

def get_args():
	parser = argparse.ArgumentParser(
		description='')
	parser.add_argument("ST", help="ST of isolate")
	parser.add_argument("--id", help= "Add this if you want strain ID in output", default="NA" )
	parser.add_argument("--reads", nargs=2, help="Two .fastq.gz read files", type=argparse.FileType('r'))
	return parser.parse_args()

if __name__ == "__main__":
	args = get_args()
	dtartrate_result = dtartrate(args.reads[0].name, args.reads[1].name)
	subspecies_result = subspecies(args.ST)
	print("\t".join([args.id, args.ST, "".join(dtartrate_result[0]), subspecies_result[0]]))
