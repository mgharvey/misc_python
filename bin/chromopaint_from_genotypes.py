#!/usr/bin/env python

"""

Name: chromopaint_from_genotypes.py 

Author: Michael G. Harvey
Date: 12 July 2013

Convert genotype probabilities file output by Tom White's post-UNEAK processing scripts to input 
file for fineSTRUCTURE (Lawson et al. 2012).

Usage: 	python chromopaint_from_genotypes.py in_file out_file sample_size

Ex.: python chromopaint_from_genotypes.py HapMapPairedFilt.txt \
		HapMapPairedFiltStruct.txt 73
		
Note that if you are using a donor list input file for fineSTRUCTURE, your samples have to be 
ordered by population before running fineSTRUCTURE. 

"""


import os
import sys
import argparse


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
	snps = 0
	missing = 0
	for line in infile:
		parts = line.split()
		alleles = str(parts[3]).split('/')
		a1 = alleles[0]
		a2 = alleles[1]
		seq = list()
		for i in range(int(sample_size)):
			if parts[4+i] == "NA":
				seq.append("-9,-9")
			else:
				counts = str(parts[4+i]).split(',')
				if counts[0] == "1":
					counts1 = "0"
				elif counts[0] == "2":
					counts1 = "1"
				else:
					print "ERROR: SNP {0} contains value that is not '1' or '2'.".format(i)
				if counts[1] == "1":
					counts2 = "0"
				elif counts[0] == "2":
					counts2 = "1"
				else:
					print "ERROR: SNP {0} contains value that is not '1' or '2'.".format(i)			
				seq.append("{0},{1}".format(counts1, counts2))
		snps += 1
		if "-9,-9" in seq:
			print "SNP {0} contains missing data and is being omitted.".format(snps)
			missing += 1
		else:
			array.append(seq)
	print "{0} of {1} SNPs contain complete data and are being written to output.".format((snps-missing), snps)
	return array


def main():
	args = get_args()
	infile = open("{0}".format(args.in_file), 'r')
	outfile = open("{0}".format(args.out_file), 'wb')
	samples = read_samples(infile, args.sample_size)
	array = unphase(infile, args.sample_size)
	alignment = zip(*array)
	i = 0
	outfile.write("0\n{0}\n{1}\nP".format(args.sample_size, len(alignment[0])))	
	for j in range(len(alignment[0])):
		outfile.write("\s{0}".format(j+1))
	outfile.write("\n")
	for j in range(len(alignment[0])):
		outfile.write("S".format(j))
	outfile.write("\n")
	for sample in samples:
		seq1 = list()
		seq2 = list()
		for j in range(len(alignment[0])):
			counts = str(alignment[i][j]).split(',')
			count1 = counts[0]
			count2 = counts[1]
			seq1.append("{0}".format(count1))
			seq2.append("{0}".format(count2))
		outfile.write("{0}\n".format(''.join(seq1)))
		outfile.write("{0}\n".format(''.join(seq2)))
		i += 1
	infile.close()
	outfile.close()

if __name__ == '__main__':
    main()