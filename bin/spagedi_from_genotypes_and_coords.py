#!/usr/bin/env python

"""

Name: spagedi_from_genotypes_and_coords.py 

Author: Michael G. Harvey
Date: 12 July 2013

Convert genotype probabilities file output by Tom White's post-UNEAK processing scripts to input 
file for SPAGeDi (Hardy and Vekemans 2002).

Usage: 	python spagedi_from_genotypes_and_coords.py in_file out_file sample_size

Ex.: python spagedi_from_genotypes_and_coords.py HapMapPairedFilt.txt \
		HapMapPairedFiltspage.txt 73

This script generates an input file for SPAGeDi specifically for the purpose of generating distance
matrices of genetic coefficients and geographic distances. Different analyses in SPAGeDi may require
different information in the input file than is included here.
	
The input file of coordinates should be a tab-delimited text file with three columns: the leftmost 
with the sample names (identical to the portion of the names before the first underscore in the 
genotypes file), the middle column with the latitude, and the rightmost column with the longitude. 
Latitude and longitude should be in decimal degrees. Example:

XM1 	-20.3461	0.4755
XM2		-21.8171	1.4889
...etc.

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
			"coord_file",
			type=str,
			help="""A tab-delimited file of coordinates"""
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


def read_coords(coords):
	coord_array = list()
	for line in coords:
		ind = line.split()
		coord_array.append(ind)
	return coord_array	
		
		
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
	array = list()
	for line in infile:
		parts = line.split()
		alleles = str(parts[3]).split('/')
		a1 = alleles[0]
		a2 = alleles[1]
		seq = list()		
		for i in range(int(sample_size)):
			if parts[4+i] == "NA":
				seq.append("0,0")
			else:
				counts = str(parts[4+i]).split(',')
				seq.append("{0},{1}".format(counts[0], counts[1]))
		array.append(seq)
	return array


def main():
	args = get_args()
	infile = open("{0}".format(args.genotype_file), 'r')
	coords = open("{0}".format(args.coord_file), 'r')
	coord_array = read_coords(coords)
	outfile = open("{0}".format(args.out_file), 'wb')
	samples = read_samples(infile, args.sample_size)
	if len(coord_array) != len(samples):
		print "ERROR: Numbers of samples ({0}) and coordinates ({1}) differ. Check file formats.".format(len(samples), len(coord_array))
	else:
		array = process_alleles(infile, args.sample_size)
		alignment = zip(*array)
		i = 0
		outfile.write("{0}\t0\t-2\t{1}\t1\t2\n".format(args.sample_size, len(alignment[0])))	
		outfile.write("0\n")	
		outfile.write("Ind\tLat\tLong")	
		for j in range(len(alignment[0])):
			outfile.write("\tloc{0}".format(j+1))
		outfile.write("\n")
		for k, sample in enumerate(samples): # Loop through samples
			outfile.write("{0}".format(sample))
			for l in range(len(coord_array)): # Loop through samples
				if sample == coord_array[l][0]:
					print sample, coord_array[l][0]
					outfile.write("\t{0}\t{1}".format(coord_array[l][1], coord_array[l][2]))
			for j in range(len(alignment[0])):
				counts = str(alignment[i][j]).split(',')
				count1 = counts[0]
				count2 = counts[1]
				outfile.write("\t{0},{1}".format(count1, count2))
			outfile.write("\n")
			i += 1
		outfile.write("END")	
	coords.close()
	infile.close()
	outfile.close()

if __name__ == '__main__':
    main()