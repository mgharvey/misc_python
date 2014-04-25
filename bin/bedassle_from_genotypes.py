#!/usr/bin/env python

"""

Name: bedassle_from_genotypes.py 

Author: Michael G. Harvey
Date: 26 February 2014

Convert genotype probabilities file output by Tom White's post-UNEAK processing scripts to allale
and sample size input files for BEDASSLE (Bradburd et al. 2013).

Usage: 	python bedassle_from_genotypes.py in_file a_outfile s_outfile sample_size

Ex.: python bedassle_from_genotypes.py HapMapPairedFilt.txt \
		HapMapPairedFiltBedassleA.txt HapMapPairedFiltBedassleB.txt 73


"""


import os
import sys
import argparse


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"genotype_file",
			type=str,
			help="""The input genotype probabilities file from Tom White's scripts"""
		)
	parser.add_argument(
			"a_outfile",
			type=str,
			help="""The file name"""
		)
	parser.add_argument(
			"s_outfile",
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


def process_alleles(infile, sample_size):	
	aarray = list()
	sarray = list()
	for line in infile:
		parts = line.split()
		alleles = str(parts[3]).split('/')
		a1 = alleles[0]
		a2 = alleles[1]
		aseq = list() # allele list		
		sseq = list() # sample size list
		for i in range(int(sample_size)):
			if parts[4+i] == "NA":
				aseq.append("0")
				sseq.append("0")
			else:
				counts = str(parts[4+i]).split(',')
				aseq.append("{0}".format(counts.count('1')))
				sseq.append('2')
		aarray.append(aseq)
		sarray.append(sseq)
	return aarray, sarray


def main():
	args = get_args()
	infile = open("{0}".format(args.genotype_file), 'r')
	a_outfile = open("{0}".format(args.a_outfile), 'wb')
	s_outfile = open("{0}".format(args.s_outfile), 'wb')
	samples = read_samples(infile, args.sample_size)
	arrays = process_alleles(infile, args.sample_size)
	aarray = arrays[0]
	sarray = arrays[1]
	aalignment = zip(*aarray)
	salignment = zip(*sarray)
	i = 0
	a_outfile.write("Ind")	
	for j in range(len(aalignment[0])):
		a_outfile.write("\tloc{0}".format(j+1))
	a_outfile.write("\n")
	s_outfile.write("Ind")	
	for j in range(len(salignment[0])):
		s_outfile.write("\tloc{0}".format(j+1))
	s_outfile.write("\n")
	for k, sample in enumerate(samples): # Loop through samples
		a_outfile.write("{0}".format(sample))
		s_outfile.write("{0}".format(sample))
		for j in range(len(aalignment[0])):
			acount = str(aalignment[i][j])
			scount = str(salignment[i][j])
			a_outfile.write("\t{0}".format(acount))
			s_outfile.write("\t{0}".format(scount))
		a_outfile.write("\n")
		s_outfile.write("\n")
		i += 1
	infile.close()
	a_outfile.close()
	s_outfile.close()

if __name__ == '__main__':
    main()