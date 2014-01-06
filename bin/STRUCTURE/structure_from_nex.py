#!/usr/bin/env python

"""

Name: structure_from_nex.py 

Author: Michael G. Harvey
Date: 12 July 2013

Convert a Nexus file of concatenated SNPs to input file for STRUCTURE (Pritchard et al. 2000).

Usage: 	python structure_from_nex.py in_file out_file sample_size

Ex.: python structure_from_nex.py HapMapPairedFilt.txt \
		HapMapPairedFiltStruct.txt 73

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
			help="""The input genotype probabilities file from Tom White's scripts"""
		)
	parser.add_argument(
			"out_file",
			type=str,
			help="""The file name"""
		)
	parser.add_argument(
			"sample_size",
			type=str,
			help="""The number of samples/individuals in file"""
		)
	return parser.parse_args()


def reformat(alignment):
	array = list()
	for w in xrange(alignment.get_alignment_length()):
		j = 0
		bases = list(alignment[:, w])
		bbases = list()
		uniqs = list(set(bases))
		buniqs = list()
		for uniq in uniqs:
			if uniq not in ["-", "?", "N"]:
				buniqs.append(uniq)	
		if len(buniqs) == 2:
			allele1 = buniqs[0]
			allele2 = buniqs[1]
			for base in bases:
				if base == allele1:
					bbases.append("1")
				elif base == allele2:
					bbases.append("2")
				elif base in ["-", "?", "N"]:
					bbases.append("-9")
			a_bases = list()
			b_bases = list()
			for i in range(len(bbases)/2):
				a_bases.append(bbases[j])
				b_bases.append(bbases[j+1])
				j += 2
		if len(buniqs) == 3:
			allele1 = buniqs[0]
			allele2 = buniqs[1]
			allele3 = buniqs[2]
			for base in bases:
				if base == allele1:
					bbases.append("1")
				elif base == allele2:
					bbases.append("2")
				elif base == allele3:
					bbases.append("3")
				elif base in ["-", "?", "N"]:
					bbases.append("-9")
			a_bases = list()
			b_bases = list()
			for i in range(len(bbases)/2):
				a_bases.append(bbases[j])
				b_bases.append(bbases[j+1])
				j += 2
		if len(buniqs) == 4:
			allele1 = buniqs[0]
			allele2 = buniqs[1]
			allele3 = buniqs[2]
			allele4 = buniqs[3]
			for base in bases:
				if base == allele1:
					bbases.append("1")
				elif base == allele2:
					bbases.append("2")
				elif base == allele3:
					bbases.append("3")
				elif base == allele4:
					bbases.append("4")	
				elif base in ["-", "?", "N"]:
					bbases.append("-9")
			a_bases = list()
			b_bases = list()
			for i in range(len(bbases)/2):
				a_bases.append(bbases[j])
				b_bases.append(bbases[j+1])
				j += 2
		array.append(a_bases)
		array.append(b_bases)
	return array
	

def main():
	args = get_args()
	alignment = AlignIO.read("{0}".format(args.in_file), "nexus")
	outfile = open("{0}".format(args.out_file), 'wb')
	samples = list()
	bsamples = list()
	for record in alignment:
		samples.append(record.id)
	k = 0
	for sample in samples:
		if k % 2 == 0:
			bsamples.append(sample)
		k += 1
	ra = reformat(alignment)
	for j in range(len(ra)/2):
		outfile.write("\tloc_{0}".format(j+1))
	outfile.write("\n")
	for i, bsample in enumerate(bsamples):
		outfile.write("{0}".format(bsample))
		j = 0
		for k in range(len(ra)/2):
			count1 = ra[j][i]
			count2 = ra[j+1][i]
			outfile.write("\t{0}\t{1}".format(count1, count2))
			j += 2
		outfile.write("\n")
	outfile.close()

if __name__ == '__main__':
    main()