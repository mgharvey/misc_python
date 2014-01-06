#!/usr/bin/env python

"""

Name: split_fasta_to_two.py 

Author: Michael G. Harvey
Date: 14 October 2013

Description: split fasta files into two to separate trans- and cis-Andean individuals for separate
analyses in COMPUTE


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
			"out_dir1",
			type=str,
			help="""The output directory 1"""
		)
	parser.add_argument(
			"out_dir2",
			type=str,
			help="""The output directory 2"""
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
		alignment = AlignIO.read("{0}{1}".format(args.in_dir, file), "fasta")		
		alignment1 = alignment[:6,:]
		alignment1.append(alignment[14,:])
		alignment1.append(alignment[15,:])
		alignment2 = alignment[6:14,:]
		print alignment2
		AlignIO.write(alignment1, "{0}trans_{1}".format(args.out_dir1, file), "fasta")
		AlignIO.write(alignment2, "{0}cis_{1}".format(args.out_dir2, file), "fasta")


if __name__ == '__main__':
    main()
