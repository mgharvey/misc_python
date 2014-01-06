#!/usr/bin/env python

"""

Name: run_mrbayes.py 

Author: Michael G. Harvey
Date: 6 October 2013

Description: Run MrBayes on directory of Nexus files.

Usage: python run_mrbayes.py in_file


"""

import os
import sys
import argparse


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_dir",
			type=str,
			help="""The input directory of nexus files"""
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
		os.system("mb {0}{1}".format(args.in_dir, file))

if __name__ == '__main__':
    main()
