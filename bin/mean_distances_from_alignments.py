#!/usr/bin/env python

"""

Name: mean_distance_from_alignments.py

Author: Michael G. Harvey
Date: 21 July 2014

Description: Calculate mean pairwise distances between samples in a dataset from multiple
alignments. Not all samples have to be in each alignment (each distance will be calculated
only from those alignments with the respective samples). Both uncorrected (p) distances 
and jukes-cantor distances are calculated.

Usage:

python mean_distance_from_alignments.py in_dir format

Where in_dir is the location of your alignments, and format is "nexus", "fasta", "phylip", 
"clustal", "emboss", "fasta-m10", "ig", "phylip-sequential", "phylip-relaxed", or 
"stockholm".

Requirements:

BioPython (Cock et al. 2009) - http://biopython.org/wiki/Main_Page


"""

import os 
import sys
import argparse
import math
import numpy as np
import numpy.ma as ma
from Bio import AlignIO


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_dir",
			type=str,
			help="""The input directory of alignments"""
		)
	parser.add_argument(
			"format",
			type=str,
			help="""Format of input files"""
		)
	return parser.parse_args()


def get_dist(dict, all_samples, p_array, jc_array, k):
	missing = ["n", "?", "-"]
	homozyg = ["a", "g", "c", "t"]
	ambig = ["m", "r", "w", "s", "y", "k"]
	for i, sample1 in enumerate(all_samples):
		for j, sample2 in enumerate(all_samples):
			if j >= i:
				break # Stop here for symmetrical matrix
			p_distprop = float(0)
			if sample1 not in dict.keys() or sample2 not in dict.keys():
				p_distprop = -1 # Set distances to missing samples to -1
				jc_distprop = -1				
			else:
				p_dist = float(0)
				seq_a = dict[sample1]
				seq_b = dict[sample2]
				for seq_x, seq_y in zip(seq_a, seq_b):
					if seq_x not in missing and seq_y not in missing:
						if seq_x != seq_y:
							if seq_x in homozyg and seq_y in homozyg:
								#print seq_x, seq_y
								p_dist += 1
							# Deal with p_distances for heterozygous sites
							elif seq_x == "m" and seq_y in ["a", "c"]:
								p_dist += 0.5
							elif seq_x == "r" and seq_y in ["a", "g"]:
								p_dist += 0.5
							elif seq_x == "w" and seq_y in ["a", "t"]:
								p_dist += 0.5
							elif seq_x == "s" and seq_y in ["c", "g"]:
								p_dist += 0.5
							elif seq_x == "y" and seq_y in ["c", "t"]:
								p_dist += 0.5
							elif seq_x == "k" and seq_y in ["g", "t"]:
								p_dist += 0.5					
							elif seq_y == "m" and seq_x in ["a", "c"]:
								p_dist += 0.5
							elif seq_y == "r" and seq_x in ["a", "g"]:
								p_dist += 0.5
							elif seq_y == "w" and seq_x in ["a", "t"]:
								p_dist += 0.5
							elif seq_y == "s" and seq_x in ["c", "g"]:
								p_dist += 0.5
							elif seq_y == "y" and seq_x in ["c", "t"]:
								p_dist += 0.5
							elif seq_y == "k" and seq_x in ["g", "t"]:
								p_dist += 0.5			
							elif seq_x in homozyg or ambig and seq_y in homozyg or ambig:
								p_dist += 1
							else:
								print "Error: {0} or {1} not a \
								recognized nucleotide".format(seq_x, seq_y)	
				p_distprop = p_dist/len(seq_a)		
				jc_distprop = -0.75*math.log(1.0 - 4.0*p_distprop/3.0) # Jukes-Cantor eqn
			p_array[i, j, k] = p_distprop # Get dist as a proportion of alignment length
			p_array[j, i, k] = p_distprop
			jc_array[i, j, k] = jc_distprop
			jc_array[j, i, k] = jc_distprop
			p_array_ma = ma.masked_where(p_array == -1, p_array) # Mask missing values
			jc_array_ma = ma.masked_where(jc_array == -1, jc_array) 
	return p_array_ma, jc_array_ma
			

def main():
	args = get_args()
	formats = ["nexus", "fasta", "phylip", "clustal", "emboss", "fasta-m10", "ig", \
	"phylip-sequential", "phylip-relaxed", "stockholm"]
	if args.format.lower() not in formats:
		raise Exception("Unrecognized format: {0}".format(args.format))
	files = list()
	prefiles = os.listdir(args.in_dir)
	for prefile in prefiles:
		if not prefile.startswith('.'):
			files.append(prefile)
	p_distances = open("./p_distances.txt", 'wb')
	jc_distances = open("./jc_distances.txt", 'wb')
	all_samples = list()
	for file in files: # Loop to get list of all samples in dataset
		if os.stat("{0}/{1}".format(args.in_dir, file))[6] == 0:	
			print "{0} is empty".format(file)
		else:
			alignment = AlignIO.read("{0}{1}".format(args.in_dir, file), \
			args.format.lower())
			for record in alignment:
				if record.id not in all_samples:
					all_samples.append(record.id)
	p_array = ma.masked_all((len(all_samples), len(all_samples), len(files)))
	jc_array = ma.masked_all((len(all_samples), len(all_samples), len(files)))
	for k, file in enumerate(files): # Loop to generate key and call distance functions
		if os.stat("{0}/{1}".format(args.in_dir, file))[6] != 0:	
			alignment = AlignIO.read("{0}{1}".format(args.in_dir, file), \
			args.format.lower())
			samples = list()
			seqs = list()
			for record in alignment:
				samples.append(str(record.id))
				seqs.append(str(record.seq))	
			sample_key = dict(zip(samples, seqs)) # Dictionary from alignment		
			arrays = get_dist(sample_key, all_samples, p_array, jc_array, k) 
			p_array = arrays[0]
			jc_array = arrays[1]
	p_mean = p_array.mean(axis=2) # Mean for each distance across loci
	jc_mean = jc_array.mean(axis=2)
	
	# Write mean distances to output file	
	p_distances.write("\t")
	jc_distances.write("\t")
	for all_sample in all_samples:
		p_distances.write("{0}\t".format(all_sample))
		jc_distances.write("{0}\t".format(all_sample))
	p_distances.write("\n")
	jc_distances.write("\n")
	for i, all_sample1 in enumerate(all_samples):
		p_distances.write("{0}\t".format(all_sample1))
		jc_distances.write("{0}\t".format(all_sample1))
		for j, all_sample2 in enumerate(all_samples):
			p_distances.write(str(p_mean.item((i, j)))) # Get distances from matrices
			jc_distances.write(str(jc_mean.item((i, j))))
			p_distances.write("\t")
			jc_distances.write("\t")
		p_distances.write("\n")
		jc_distances.write("\n")
	p_distances.close()
	jc_distances.close()
		
if __name__ == '__main__':
    main()