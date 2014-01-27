#!/usr/bin/env python

"""

Name: bayescan_from_nex.py 

Author: Michael G. Harvey
Date: 4 September 2013

Description: Convert a nexus alignment of concatenated SNPs to the input file format for 
bayescan (Foll and Gaggiotti 2008).

Usage: 	python bayescan_from_nex.py input_file output_directory number_of_populations \
			individuals_in_pop_1 individuals_in_pop_2 ... individuals_in_pop_n

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
	out = open("{0}bayescan_input.txt".format(args.out_dir), 'wb')
	out.write("[loci]={0}\n\n".format(alignment.get_alignment_length()))
	out.write("[populations]={0}\n\n".format(args.populations))
	total_size = 0
	for x in range(args.populations):
		out.write("[pop]={0}\n".format(x+1)) 
		for w in xrange(alignment.get_alignment_length()):
			out.write("{0}\t".format(w+1))
			bases = alignment[:, w]
			uniqs = list(set(bases))
			nuniqs = filter(lambda a: a != "N", uniqs) # Filter out N's to assign alleles
			if len(nuniqs) != 2:
				print "ERROR: Locus {0} not biallelic. Outfile invalid.".format(w)
			allele1 = nuniqs[0]
			allele2 = nuniqs[1]
			pop_size = args.pop_sizes[x]
			pop_bases = alignment[total_size:(total_size+pop_size), w]
			a = 0
			b = 0
			n = 0
			for pop_base in pop_bases:
				if pop_base == allele1:
					a += 1
				elif pop_base == allele2:
					b += 1
				elif pop_base == "N":
					n += 1
			if a+b+n != pop_size:
				print "ERROR: Number of alleles not equal to diploid population size."
			out.write("{0}\t".format(args.pop_sizes[x]-n))
			out.write("2\t")
			out.write("{0}\t{1}\n".format(a,b))
		out.write("\n")
		total_size += pop_size
		out.flush()
	out.close()
	
if __name__ == '__main__':
    main()
