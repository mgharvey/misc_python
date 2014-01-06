#!/usr/bin/env python

"""

Name: fasta_from_nex.py 

Author: Michael G. Harvey
Date: 14 October 2013

Description: Convert Nexus file to fasta file.

Usage: python fasta_from_nex.py in_file out_dir


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
			help="""The input nexus file"""
		)
	parser.add_argument(
			"out_dir",
			type=str,
			help="""The output directory"""
		)
	return parser.parse_args()


def main():
	args = get_args()
	files = list()
	prefiles = os.listdir(args.in_dir)
	for prefile in prefiles: # Remove hidden files
		if not prefile.startswith('.'):
			files.append(prefile)
	os.chdir(args.in_dir)
	for file in files:
		print file
		alignment = AlignIO.read("{0}{1}".format(args.in_dir, file), "nexus")
		AlignIO.write(alignment, "{0}{1}".format(args.out_dir, file), "fasta")

if __name__ == '__main__':
    main()
