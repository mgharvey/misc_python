#!/usr/bin/env python

"""

Name: run_mbsum.py 

Author: Michael G. Harvey
Date: 8 October 2013

Description: Run mbsum (BUCKy) on a directory of MrBayes output files.

Usage: python run_mbsum.py in_dir out_dir burnin


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
			help="""The input directory of MrBayes output files"""
		)
	parser.add_argument(
			"out_dir",
			type=str,
			help="""The output directory"""
		)
	parser.add_argument(
			"burnin",
			type=str,
			help="""The number of trees to discard as burn-in"""
		)
	return parser.parse_args()


def main():
	args = get_args()
	files = list()
	prefiles = os.listdir(args.in_dir)
	for prefile in prefiles:
		if not prefile.startswith('.'): # Remove hidden files
			if prefile.endswith('1.t'): # Only take MrBayes treefiles
				files.append(prefile)
	os.chdir(args.in_dir)
	for file in files:
		parts = file.split('.')
		name = list()
		for part in parts[:-2]:
			name.append("{0}".format(part))
			name.append(".")
		middlename = ''.join(name)
		outname = middlename[:-1]
		newname = "{0}run?.t".format(middlename)
		os.system("mbsum -n {0} -o {1}{2}.in {3}".format(args.burnin, args.out_dir, outname, newname))
	
if __name__ == '__main__':
    main()
