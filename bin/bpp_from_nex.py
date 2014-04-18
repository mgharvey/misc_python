#!/usr/bin/env python

"""

Name: bpp_from_nex.py 

Author: Michael G. Harvey
Date: 14 April 2014

Description: Convert folder of Nexus files to BP&P input.

Usage: python bpp_from_nex.py in_dir out_file


"""

import os
import sys
import argparse
from subprocess import check_call, call, PIPE
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
			"out_file",
			type=str,
			help="""The output directory"""
		)
	return parser.parse_args()


def main():
	args = get_args()
	prefiles = list()
	files = list()
	prefiles = os.listdir(args.in_dir)
	for prefile in prefiles: # Remove hidden files
		if not prefile.startswith('.'):
			if not prefile.startswith('tmp'):
				files.append(prefile)
	os.chdir(args.in_dir)
	out = open("{0}".format(args.out_file), 'wb')
	for file in files:
		alignment = AlignIO.read("{0}{1}".format(args.in_dir, file), "nexus")
		out.write("{0}  {1}\n".format(len(alignment), alignment.get_alignment_length()))
		i = 1
		for record in alignment:
			out.write("{0}^{1}  {2}\n".format(record.id, i, record.seq))
			i += 1
	out.close()


if __name__ == '__main__':
    main()