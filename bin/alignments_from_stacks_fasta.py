#!/usr/bin/env python

"""

Name: alignments_from_stacks_fasta.py

Author: Michael G. Harvey
Date: 6 June 2014

Description: Convert the fasta output from the populations module of the Stacks (Catchen et al. 
2013) pipeline to alignments for coalescent analyses. Designed for use with diploid data. For 
convenience, files can be output in fasta, phylip-like (no limits on sample name length), or nexus 
format. There is an option to remove those alignments containing only a single sequence.

Usage: python alignments_from_stacks_fasta.py in_file [--fasta] [--phylip] [--nexus] [--remove_singletons]


"""

import os
import sys
import random
import argparse
import dendropy
from dendropy import popgenstat
from Bio import AlignIO


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_file",
			type=str,
			help="""The input file of fasta sequences from stacks populations module"""
		)
	parser.add_argument(
			"--fasta",
			action='store_true',
			help="""Output alignments in fasta format"""
		)
	parser.add_argument(
			"--phylip",
			action='store_true',
			help="""Output alignments in phylip-like format"""
		)
	parser.add_argument(
			"--nexus",
			action='store_true',
			help="""Output alignments in nexus format"""
		)
	parser.add_argument(
			"--remove_singletons",
			action='store_true',
			help="""Remove alignments containing only a single sequence"""
		)
	return parser.parse_args()


def fasta(align_rows, locus):
	out = open("./out/locus_{0}.fa".format(locus), 'wb')
	samples = list()
	for align_row in align_rows: # Make a list of samples in this locus
		samples.append(str(align_row[0]))		
	second = False
	for align_row in align_rows:
		if samples.count(str(align_row[0])) == 1: # If sample is only present 1x
			out.write(">Sample_{0}a\n{1}\n".format(align_row[0], align_row[1]))
		elif samples.count(str(align_row[0])) == 2: # If sample is present 2x (het)
			if second == False:
				out.write(">Sample_{0}a\n{1}\n".format(align_row[0], align_row[1]))
				second = True
			elif second == True:
				out.write(">Sample_{0}b\n{1}\n".format(align_row[0], align_row[1]))
				second = False
		elif samples.count(str(align_row[0])) > 2:
			print "Sample {0} contains >2 alleles. Exiting.".format(align_row[0])
			exit()
	out.close()			
	align_rows = list()
	
	
def phylip(align_rows, locus):
	out = open("./out/locus_{0}.phy".format(locus), 'wb')
	out.write(" {0} {1}\n".format(len(align_rows), len(align_rows[0][1])))
	samples = list()
	for align_row in align_rows: # Make a list of samples in this locus
		samples.append(str(align_row[0]))		
	second = False
	for align_row in align_rows:
		if samples.count(str(align_row[0])) == 1: # If sample is only present 1x
			out.write("Sample_{0}a  {1}\n".format(align_row[0], align_row[1]))
		elif samples.count(str(align_row[0])) == 2: # If sample is present 2x (het)
			if second == False:
				out.write("Sample_{0}a  {1}\n".format(align_row[0], align_row[1]))
				second = True
			elif second == True:
				out.write("Sample_{0}b  {1}\n".format(align_row[0], align_row[1]))
				second = False
		elif samples.count(str(align_row[0])) > 2:
			print "Sample {0} contains >2 alleles. Exiting.".format(align_row[0])
			exit()
	out.close()			
	align_rows = list()


def nexus(align_rows, locus):
	out = open("./out/locus_{0}.nex".format(locus), 'wb')
	out.write("#NEXUS\n")
	out.write("Begin data;\n")
	out.write("Dimensions ntax={0} nchar={1}\n".format(len(align_rows), len(align_rows[0][1])))
	out.write("Format datatype=dna symbols=\"ACTG\" missing=? gap=-;samples = list()\n")
	out.write("Matrix;\n")
	samples = list()
	for align_row in align_rows: # Make a list of samples in this locus
		samples.append(str(align_row[0]))		
	second = False
	for align_row in align_rows:
		if samples.count(str(align_row[0])) == 1: # If sample is only present 1x
			out.write("Sample_{0}a	{1}\n".format(align_row[0], align_row[1]))
		elif samples.count(str(align_row[0])) == 2: # If sample is present 2x (het)
			if second == False:
				out.write("Sample_{0}a	{1}\n".format(align_row[0], align_row[1]))
				second = True
			elif second == True:
				out.write("Sample_{0}b	{1}\n".format(align_row[0], align_row[1]))
				second = False
		elif samples.count(str(align_row[0])) > 2:
			print "Sample {0} contains >2 alleles. Exiting.".format(align_row[0])
			exit()
	out.write(";\n")
	out.write("End;")
	out.close()			
	align_rows = list()


def main():
	args = get_args()
	infile = open("{0}".format(args.in_file), 'r')
	os.system("mkdir out")
	prev_parts = infile.readline().split('_') # Get number of first locus in infile
	prev_locus = int(prev_parts[1])
	infile.seek(0)
	align_rows = list() # Create an empty list for sample and sequence info
	i = 1
	for line in infile:
		if line.startswith('>'):
			parts = line.split('_')
			locus = int(parts[1])
			if locus != prev_locus: # If this upcoming locus is a new locus, write the current one
				if args.remove_singletons: # Remove loci with only one sequence
					if len(align_rows) > 1: 
						if args.fasta:
							fasta(align_rows, prev_locus)
						if args.phylip:
							phylip(align_rows, prev_locus)
						if args.nexus:
							nexus(align_rows, prev_locus)					
				else: # Output all loci (including those with only one sequence)
					if args.fasta:
						fasta(align_rows, prev_locus)
					if args.phylip:
						phylip(align_rows, prev_locus)
					if args.nexus:
						nexus(align_rows, prev_locus)					
				i += 1
				align_rows = list()
			sample = int(parts[3])
			seq = next(infile).rstrip()
			align_rows.append([sample, seq])
			prev_locus = locus
	if args.fasta:
		fasta(align_rows, prev_locus)
	if args.phylip:
		phylip(align_rows, prev_locus)
	if args.nexus:
		nexus(align_rows, prev_locus)					
	i += 1
	infile.close()
	print "Loci processed: {0}".format(i) # Total number of loci examined

if __name__ == '__main__':
    main()