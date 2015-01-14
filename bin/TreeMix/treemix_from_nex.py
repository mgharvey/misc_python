#!/usr/bin/env python

"""

Name: treemix_from_nex.py 
Author: Michael G. Harvey
Date: 12 May 2013

Description: Convert a nexus alignment of concatenated SNPs to the input file format for 
Treemix (Pickrell and Pritchard 2012).

Usage: 	python treemix_from_nex.py input_file output_directory number_of_populations \
			individuals_in_pop_1 individuals_in_pop_2 ... individuals_in_pop_n

Note that TreeMix uses only biallelic SNPs, so any sites that do not represent biallelic SNPs, will 
not be used.

"""

import os
import sys
import argparse
from Bio import AlignIO


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_file",
			type=str,
			help="""The input nexus file"""
		)
	parser.add_argument(
			"out_dir",
			type=str,
			help="""The output directory"""
		)
	parser.add_argument(
			"populations",
			type=int,
			help="""The number of populations"""
		)
	parser.add_argument(
			"pop_sizes",
			type=int,
			nargs='+',
			help="""The number of samples in each population"""
		)
	return parser.parse_args()


def main():
	args = get_args()
	alignment = AlignIO.read("{0}".format(args.in_file), "nexus")
	out = open("{0}treemix_file_out.txt".format(args.out_dir), 'wb')
	for x in range(args.populations):
		out.write("pop{0} ".format(x+1)) 
	out.write("\n")
	total_size = 0
	for w in xrange(alignment.get_alignment_length()):
		total_size = 0
		bases = alignment[:, w] 		
		uniqs = list()
		uniqs = list(set(bases))
		nmuniqs = list()
		for uniq in uniqs: # Remove sites with missing data
			if uniq not in ['?', 'n', '?', 'N']:
				nmuniqs.append(uniq)			
		if len(nmuniqs) != 2:
			print "Skipping site {0} - not a biallelic SNP".format(w)	
		else:
			allele1 = nmuniqs[0]
			allele2 = nmuniqs[1]
			total_size = 0
			for x in range(args.populations):
				pop_size = args.pop_sizes[x]
				pop_bases = alignment[total_size:(total_size+pop_size), w]
				total_size += pop_size
				a = 0
				b = 0
				for pop_base in pop_bases:
					if pop_base == allele1:
						a += 1
					elif pop_base == allele2:
						b += 1
				out.write("{0},{1} ".format(a,b))
		out.write("\n")	
		out.flush()
	out.close()
	
if __name__ == '__main__':
    main()