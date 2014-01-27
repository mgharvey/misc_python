#!/usr/bin/env python

"""

Name: structurama_from_genotypes.py 

Author: Michael G. Harvey
Date: 11 July 2013

Convert genotype probabilities file output by Tom White's post-UNEAK processing scripts to input 
file for structurama (Huelsenbeck and Andolfatto 2007).

Usage: 	python structurama_from_genotypes.py in_file out_file sample_size

Ex.: python structurama_from_genotypes.py HapMapPairedFilt.txt \
		HapMapPairedFiltStruct.txt 73

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
	for line in infile:
		parts = line.split()
		alleles = str(parts[3]).split('/')
		a1 = alleles[0]
		a2 = alleles[1]
		seq = list()		
		for i in range(int(sample_size)):
			if parts[4+i] == "NA":
				seq.append("?,?")
			else:
				counts = str(parts[4+i]).split(',')
				seq.append("{0},{1}".format(counts[0], counts[1]))
		array.append(seq)
	return array


def main():
	args = get_args()
	infile = open("{0}".format(args.in_file), 'r')
	outfile = open("{0}".format(args.out_file), 'wb')
	samples = read_samples(infile, args.sample_size)
	array = unphase(infile, args.sample_size)
	alignment = zip(*array)
	outfile.write("begin data;\n")
	outfile.write("\tdimensions nind={0} nloci={1};\n".format(args.sample_size, len(alignment[0])))
	outfile.write("\tinfo\n")
	i = 0
	for sample in samples:
		outfile.write("\t{0}".format(sample))
		for j in range(len(alignment[0])):
			counts = str(alignment[i][j]).split(',')
			count1 = counts[0]
			count2 = counts[1]
			outfile.write(" ( {0} , {1} )".format(count1, count2))
		if i < (len(samples)-1):
			outfile.write(",")
		outfile.write("\n")
		i += 1
	outfile.write("\t;\n")
	outfile.write("end;\n")
	infile.close()
	outfile.close()

if __name__ == '__main__':
    main()