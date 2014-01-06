#!/usr/bin/env python

"""

Name: split_diplotypes_from_nexus.py 

Author: Michael G. Harvey
Date: 17 June 2013

Description: Split diploid genotypes into haplotypes assuming no linkage.

Usage: python split_diplotypes_from_nexus.py in_file out_dir


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
			help="""The input directory"""
		)	
	parser.add_argument(
			"out_dir",
			type=str,
			help="""The output directory"""
		)
	return parser.parse_args()


def phase(alignment):
	new_alignment = list()
	for w in xrange(alignment.get_alignment_length()):
		bases = alignment[:,w]
		new_bases = list()
		for base in bases:
			if base in ["A","a"]:
				new_bases.append("A")
				new_bases.append("A")
			if base in ["C","c"]:
				new_bases.append("C")
				new_bases.append("C")
			if base in ["G","g"]:
				new_bases.append("G")
				new_bases.append("G")
			if base in ["T","t"]:
				new_bases.append("T")
				new_bases.append("T")
			if base in ["M","m"]:
				new_bases.append("A")
				new_bases.append("C")
			if base in ["R","r"]:
				new_bases.append("A")
				new_bases.append("G")
			if base in ["W","w"]:
				new_bases.append("A")
				new_bases.append("T")
			if base in ["S","s"]:
				new_bases.append("C")
				new_bases.append("G")
			if base in ["Y","y"]:
				new_bases.append("C")
				new_bases.append("T")
			if base in ["K","k"]:
				new_bases.append("G")
				new_bases.append("T")
			if base in ["N","n","-","?"]:
				new_bases.append("N")
				new_bases.append("N")
		new_alignment.append(new_bases)
	return new_alignment


def main():
	args = get_args()	
	files = list()
	prefiles = os.listdir(args.in_dir)
	for prefile in prefiles: # Remove hidden files
		if not prefile.startswith('.'):
			files.append(prefile)
	os.chdir(args.in_dir)
	for file in files:	
		alignment = AlignIO.read("{0}{1}".format(args.in_dir, file), "nexus")
		samples = list()
		for record in alignment:
			samples.append(record.id)
		new_alignment = phase(alignment)
		final_alignment = zip(*new_alignment)
		ff = open("{0}{1}".format(args.out_dir,file), 'wb')
		ff.write("#NEXUS\n\n")
		ff.write("Begin data;\n")
		ff.write("\tDimensions ntax={0} nchar={1};\n".format(len(final_alignment), len(new_alignment)))
		ff.write("\tFormat datatype=dna gap=-;\n")
		ff.write("\tMatrix\n")
		j = 0
		for i in range((len(final_alignment)/2)):
			ff.write("{0}a\t".format(samples[i]))
			ff.writelines(final_alignment[:][j])
			ff.write("\n")
			ff.write("{0}b\t".format(samples[i]))
			ff.writelines(final_alignment[:][j+1])
			ff.write("\n")
			j += 2	
		ff.write(";\n")
		ff.write("End;\n")
		ff.close()	

if __name__ == '__main__':
    main()