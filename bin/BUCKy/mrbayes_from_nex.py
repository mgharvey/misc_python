#!/usr/bin/env python

"""

Name: mrbayes_from_nex.py 

Author: Michael G. Harvey
Date: 6 October 2013

Description: Add MrBayes commands to Nexus files.

Usage: python mrbayes_from_nex.py in_dir out_dir

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
	parser.add_argument(
			"out_dir",
			type=str,
			help="""The output directory of nexus files with MrBayes commands"""
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
		infile = open("{0}{1}".format(args.in_dir, file), 'r')
		outfile = open("{0}{1}".format(args.out_dir, file), 'wb')
		for line in infile:
			outfile.write(line)
		outfile.write("\n\nBegin mrbayes;\n")
		outfile.write("\tset autoclose=yes;\n")
		outfile.write("\tlset nst=6 rates=gamma;\n")
		outfile.write("\tmcmc ngen=1100000 nruns=2 nchains=4 burnin=200 starttree=random;\n")
		outfile.write("End;\n")
		infile.close()
		outfile.close()
		

if __name__ == '__main__':
    main()
