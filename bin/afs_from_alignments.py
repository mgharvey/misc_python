#!/usr/bin/env python

"""

Name: afs_from_alignments.py 

Author: Michael G. Harvey
Date: 31 July 2013

Description: Use sites with complete data in alignments to obtain allele frequency information.
Sites with missing data are removed. All allele variants with fewer copies than the total 
number of alleles are output. The input is assumed to be both split alleles from diploid 
samples in Nexus format (one following the other for each sample). The number of heterozygotes
and homozygotes based on the sites with complete data are also output.

Usage: python afs_from_alignments.py alleles

Where alleles is the total number of chromosomes in the sample (e.g., 16 for 8 diploid 
individuals).

"""

import os
import sys
import argparse
from Bio import AlignIO


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_dir",
			type=str,
			help="""The input directory of split haplotypes"""
		)
	parser.add_argument(
			"alleles",
			type=int,
			help="""Total number of alleles in dataset (individuals*2 for diploids)"""
		)
	return parser.parse_args()


def main():
	args = get_args()
	files = list()
	prefiles = os.listdir(args.in_dir)
	for prefile in prefiles:
		if not prefile.startswith('.'):
			files.append(prefile)
	number_alleles = range(1, args.alleles)
	empty_counts = [0]*(args.alleles)
	afs_dict = dict(zip(number_alleles, empty_counts))
	homo = 0
	het = 0
	for file in files:		
		# Pull out variable sites without missing data
		alleles = list() # Make empty list of alleles
		alignment = AlignIO.read("{0}{1}".format(args.in_dir, file), "nexus")
		new_alignments = list()	
		for w in xrange(alignment.get_alignment_length()):
			bases = alignment[:, w] 
			missing = False
			for base in bases:
				if base in ["n", "N", "-", "?"]:
					missing = True
			if len(bases) < 16:
				missing = True
			if missing == False:
				uniqs = list(set(bases))
				if len(uniqs) > 1:
					new_alignments.append(bases)		
		new_alignments = zip(*new_alignments) # transpose
		alleles = list()
		for new_alignment in new_alignments:
			new_alignment = ''.join(new_alignment)
			alleles.append(new_alignment)
		for i, allele in enumerate(alleles):	
			if i % 2 == 0: # If odd
				if allele == alleles[i+1]:
					homo += 1
				elif allele != alleles[i+1]:
					het += 1
		uniq_alleles = list(set(alleles))
		for uniq_allele in uniq_alleles:
			count = alleles.count(uniq_allele)
			for i in range(args.alleles):
				if count == i:
					afs_dict[i] += 1
	
	# Print dict to stdout
	print "Number of heterozygotes: {0}".format(het)
	print "Number of homozygotes: {0}".format(homo)
	print "Allele_Frequency\tCount_of_Alleles"
	for key in afs_dict:
		print "{0}\t{1}".format(key, afs_dict[key])

if __name__ == '__main__':
    main()