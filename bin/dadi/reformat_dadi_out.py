#!/usr/bin/env python

"""

Name: reformat_dadi_out.py 

Author: Michael G. Harvey
Date: 8 May 2013

Description: Reformat a custom output txt file from dadi into a spreadsheet-compatible delimited 
textfile.

Usage: python reformat_dadi_out.py in_dir out_file

"""


import os, sys
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
			"out_file",
			type=str,
			help="""The desired output file name"""
		)
	return parser.parse_args()

def main():
	args = get_args()
	dir = args.in_dir
	os.chdir(dir)
	prefiles = os.listdir(dir)
	files = list()
	for prefile in prefiles:
		if prefile.startswith("dadi"):
			files.append(prefile)
	outfile = open("{0}".format(args.out_file), 'wb')
	i = 0
	print files
	for file in files:
		infile = open("{0}{1}".format(dir, file), 'r')
		elements = list()
		for line in infile:
			if "Likelihood:" in line:
				line = line.rstrip()
				elements.append(line.replace("Likelihood: ",""))
			if "Optimized Parameters:" in line:
				line = line.rstrip()
				line = line.replace("Optimized Parameters: array([ ","")
				parts = line.split(', ')
				for part in parts:
					elements.append(part.replace(",",""))
			if "        " in line:
				line = line.rstrip()
				line = line.replace("         ","")
				line = line.replace("])","")
				parts = line.split(', ')
				for part in parts:
					elements.append(part.replace(",",""))
		for element in elements:
			outfile.write("{0}\t".format(element))
		outfile.write("\n")
		i+=1
		infile.close()
		outfile.flush()
	outfile.close()
		
if __name__ == '__main__':
    main()

