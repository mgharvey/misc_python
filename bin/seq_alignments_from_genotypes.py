#!/usr/bin/env python

"""

Name: seq_alignments_from_genotypes.py 

Author: Michael G. Harvey
Date: 1 October 2013

Convert genotype probabilities file output by Tom White's post-UNEAK processing scripts to nexus 
alignments containing full sequences from each individual for all loci.

Usage: 	python nexus_from_genotypes.py in_file out_dir sample_size

Example: python seq_alignments_from_genotype_probabilities.py HapMapPairedFilt_8_20_75.txt \
			./GBS_Seqs/ 8

"""

import os
import sys
import argparse
import csv
import numpy

def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_file",
			type=str,
			help="""The input genotype probabilities file from Tom White's scripts"""
		)
	parser.add_argument(
			"out_dir",
			type=str,
			help="""The output directory"""
		)
	parser.add_argument(
			"sample_size",
			type=str,
			help="""The number of samples/individuals in file"""
		)
	return parser.parse_args()


def read_samples(infile, sample_size):
	samples = list()
	first_line = infile.readline()
	parts = first_line.split()
	for i in range(int(sample_size)):
		parts2 = parts[4+i]
		parts3 = parts2.split('_')
		samples.append(parts3[0])	
	return samples


def unphase(infile, sample_size):	
	array = list()
	for line in infile:
		parts = line.split()		
		alleles = str(parts[3]).split('/')
		a1 = alleles[0]
		a2 = alleles[1]
		seq = list()		
		for i in range(int(sample_size)):
			if parts[4+i] == "NA":
				seq.append("N")
			else:
				alsp = parts[4+i]
				als = alsp.split(',')
				if a1 == "A":
					if a2 == "C":
						if als[0] == "1":
							if als[1] == "1":
								seq.append("A")
							elif als[1] == "2":
								seq.append("M")
						elif als[0] == "2":
							if als[1] == "2":
								seq.append("C")
							elif als[1] == "1":
								seq.append("M")
					if a2 == "G":
						if als[0] == "1":
							if als[1] == "1":
								seq.append("A")
							elif als[1] == "2":
								seq.append("R")
						elif als[0] == "2":
							if als[1] == "2":
								seq.append("G")
							elif als[1] == "1":
								seq.append("R")					
					if a2 == "T":
						if als[0] == "1":
							if als[1] == "1":
								seq.append("A")
							elif als[1] == "2":
								seq.append("W")
						elif als[0] == "2":
							if als[1] == "2":
								seq.append("T")
							elif als[1] == "1":
								seq.append("W")				
				if a1 == "C":
					if a2 == "A":
						if als[0] == "1":
							if als[1] == "1":
								seq.append("C")
							elif als[1] == "2":
								seq.append("M")
						elif als[0] == "2":
							if als[1] == "2":
								seq.append("A")
							elif als[1] == "1":
								seq.append("M")
					if a2 == "G":
						if als[0] == "1":
							if als[1] == "1":
								seq.append("C")
							elif als[1] == "2":
								seq.append("S")
						elif als[0] == "2":
							if als[1] == "2":
								seq.append("G")
							elif als[1] == "1":
								seq.append("S")					
					if a2 == "T":
						if als[0] == "1":
							if als[1] == "1":
								seq.append("C")
							elif als[1] == "2":
								seq.append("Y")
						elif als[0] == "2":
							if als[1] == "2":
								seq.append("T")
							elif als[1] == "1":
								seq.append("Y")			
				if a1 == "G":
					if a2 == "A":
						if als[0] == "1":
							if als[1] == "1":
								seq.append("G")
							elif als[1] == "2":
								seq.append("R")
						elif als[0] == "2":
							if als[1] == "2":
								seq.append("A")
							elif als[1] == "1":
								seq.append("R")
					if a2 == "C":
						if als[0] == "1":
							if als[1] == "1":
								seq.append("G")
							elif als[1] == "2":
								seq.append("S")
						elif als[0] == "2":
							if als[1] == "2":
								seq.append("C")
							elif als[1] == "1":
								seq.append("S")					
					if a2 == "T":
						if als[0] == "1":
							if als[1] == "1":
								seq.append("G")
							elif als[1] == "2":
								seq.append("K")
						elif als[0] == "2":
							if als[1] == "2":
								seq.append("T")
							elif als[1] == "1":
								seq.append("K")	
				if a1 == "T":
					if a2 == "A":
						if als[0] == "1":
							if als[1] == "1":
								seq.append("T")
							elif als[1] == "2":
								seq.append("W")
						elif als[0] == "2":
							if als[1] == "2":
								seq.append("A")
							elif als[1] == "1":
								seq.append("W")
					if a2 == "C":
						if als[0] == "1":
							if als[1] == "1":
								seq.append("T")
							elif als[1] == "2":
								seq.append("Y")
						elif als[0] == "2":
							if als[1] == "2":
								seq.append("C")
							elif als[1] == "1":
								seq.append("Y")					
					if a2 == "G":
						if als[0] == "1":
							if als[1] == "1":
								seq.append("T")
							elif als[1] == "2":
								seq.append("K")
						elif als[0] == "2":
							if als[1] == "2":
								seq.append("G")
							elif als[1] == "1":
								seq.append("K")	
		print seq
		seq1 = list(parts[1])
		seq2 = list(parts[2])
		for i in range(len(seq1)):
			if seq1[i] != seq2[i]:
				preseq = ''.join(seq1[0:i])
				print preseq
				postseq = ''.join(seq1[i:])
				print postseq
		locus_list = list()
		for i in range(len(seq)):
			seqlist = list()
			seqlist.append(preseq)
			seqlist.append(seq[i])
			seqlist.append(postseq)
			finalseq = ''.join(seqlist)
			locus_list.append(finalseq)
		array.append(locus_list)
	print array
	return array


def main():
	args = get_args()
	infile = open("{0}".format(args.in_file), 'r')
	samples = read_samples(infile, args.sample_size)
	array = unphase(infile, args.sample_size)
	for i in range(len(array)):		
		j = 0
		outfile = open("{0}locus{1}.nex".format(args.out_dir, i+1), 'wb')
		outfile.write("#NEXUS\n")
		outfile.write("begin data;\n")
		outfile.write("\tdimensions ntax={0} nchar={1};\n".format(args.sample_size, len(array[i][0])))
		outfile.write("\tformat datatype=dna missing=? gap=-;\n")
		outfile.write("\tmatrix\n")
		for sample in samples:
			outfile.write("{0}\t{1}\n".format(sample, array[i][j]))
			j += 1
		outfile.write(";\n")
		outfile.write("end;\n")
		infile.close()
		outfile.close()

if __name__ == '__main__':
    main()