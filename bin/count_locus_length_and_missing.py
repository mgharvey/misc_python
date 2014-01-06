#!/usr/bin/env python

"""

Name: count_locus_length_and_missing.py 

Author: Michael G. Harvey
Date: 29 September 2013

Description: Count length of alignments and number of sites with missing data.

Usage: python count_locus_length_and_missing.py in_dir out_file



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
			help="""The input directory of nexus files"""
		)
	parser.add_argument(
			"out_file",
			type=str,
			help="""The desired output file name"""
		)
	return parser.parse_args()


def main():
	args = get_args()
	os.chdir(args.in_dir)
	files = list()
	prefiles = os.listdir(args.in_dir)
	for prefile in prefiles: # Remove hidden files
		if not prefile.startswith('.'):
			files.append(prefile)
	out = open("{0}".format(args.out_file), 'wb')
	out.write("Alignments: {0}\n".format(len(files)))
	for i, file in enumerate(files):
		missing = 0
		polymorphic = 0
		alignment = AlignIO.read("{0}{1}".format(args.in_dir, file), "nexus")
		for record in alignment:
			bases = record.seq
			for base in bases:
				if base in ["N","n","-","?"]:
					missing += 1
		rows = len(alignment)
		for w in xrange(alignment.get_alignment_length()):
			bases = alignment[:,w] 
			nm_bases = list()
			for base in bases:
				if base not in ["N","n","-","?"]:
					nm_bases.append(base)
			unique = set(nm_bases)
			if len(unique) > 1:
				polymorphic += 1
				print unique
		out.write("{0}\t{1}\t{2}\t{3}\n".format(file, alignment.get_alignment_length(), missing, polymorphic))
		out.flush()
	out.close()
	
if __name__ == '__main__':
    main()
