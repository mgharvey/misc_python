#!/usr/bin/env python

"""

Name: run_mraic.py

Author: Michael G. Harvey
Date: 5 July 2013

Description: Run mraic.pl (Nylanderb 2004) on a folder of alignments in phylip/phyml format.

Usage: python run_mraic.py mraic_dir in_dir out_dir

python run_mraic.py /Users/michaelharvey/Applications/mraic /Users/michaelharvey/Desktop/pic/beast/deep_UCEs/77_loci_phylip ~/Desktop/mraic_out
python run_mraic.py /Users/michaelharvey/Applications/mraic /Users/michaelharvey/Desktop/pic/beast/shallow_UCEs/Xm/orthologs/phylip ~/Desktop/mraic_UCE_shallow_out

"""



import os
import sys
import argparse

def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"mraic_dir",
			type=str,
			help="""The directory for mraic.pl"""
		)
	parser.add_argument(
			"in_dir",
			type=str,
			help="""The output directory"""
		)
	parser.add_argument(
			"out_dir",
			type=str,
			help="""The output directory"""
		)
	return parser.parse_args()
	

def main():
	args = get_args()
	outfile = open("{0}/mraic_out.txt".format(args.out_dir), 'wb')
	files = list()
	prefiles = os.listdir("{0}".format(args.in_dir))
	for prefile in prefiles: # Remove hidden files
		if not prefile.startswith('.'):
			files.append(prefile)
	os.chdir("{0}".format(args.mraic_dir))
	for file in files:
		os.system("perl mraic.pl {0}/{1}".format(args.in_dir, file))
		infile = open("{0}/{1}.MrAIC.txt".format(args.in_dir, file), 'r')
		for line in infile:
			if line.startswith("Minimum AICc model:"):
				parts = line.split()
				outfile.write("{0}\t{1}\n".format(file, parts[3]))				
		infile.close()
		outfile.flush()
	outfile.close()

if __name__ == '__main__':
	main()