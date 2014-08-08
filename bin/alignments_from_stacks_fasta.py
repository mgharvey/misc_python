#!/usr/bin/env python

"""

Name: alignments_from_stacks_fasta.py

Author: Michael G. Harvey
Date: 6 June 2014

Description: Convert the fasta output from the populations module of the Stacks (Catchen et al. 
2013) pipeline to alignments for coalescent analyses. Warning messages will be printed if the
number of alleles at a locus exceeds the expected ploidy level of the samples (which must be
specified by the user), but these loci will still be printed to output. Note that Stacks will 
often output more alleles than the ploidy level if the arguments "--max_locus_stacks" and "-H" 
are not used to restrict the number of alleles in the ustacks program. For convenience, files 
can be output in fasta, phylip-like (no limits on sample name length), or nexus format. There 
is an option to remove those alignments containing only a single sequence. Haploid and diploid 
samples can be further processed with process_stacks_alignments.py (also in this repository).

Usage: python alignments_from_stacks_fasta.py in_file ploidy [--fasta] [--phylip] [--nexus] [--remove_singletons]

Where in_file is the location of the input batch*.fa file from Stacks, ploidy is an integer,
and the other options are set by typing the string in the brackets.

"""

import os
import sys
import argparse


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_file",
			type=str,
			help="""The input file of fasta sequences from stacks populations module"""
		)
	parser.add_argument(
			"ploidy",
			type=int,
			help="""The ploidy of the study samples (an integer value)"""
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


def fasta(align_rows, locus, ploidy):
	out = open("./out_fasta/locus_{0}.fa".format(locus), 'wb')
	samples = list()
	for align_row in align_rows: # Make a list of samples in this locus
		samples.append(str(align_row[0]))		

	i = 0
	allele_labs = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p']
	for align_row in align_rows:
		if samples.count(str(align_row[0])) > ploidy:
			print  "WARNING: Number of alleles > ploidy at locus {0}, sample {0}".format(locus, align_row[0]) 
		if i < (samples.count(str(align_row[0]))-1):		
			out.write(">Sample_{0}{1}\n{2}\n".format(align_row[0], allele_labs[i], align_row[1]))
			i += 1
		elif i == (samples.count(str(align_row[0]))-1):
			out.write(">Sample_{0}{1}\n{2}\n".format(align_row[0], allele_labs[i], align_row[1]))
			i = 0
	out.close()			
	align_rows = list()


def phylip(align_rows, locus, ploidy):
	out = open("./out_phylip/locus_{0}.phy".format(locus), 'wb')
	out.write(" {0} {1}\n".format(len(align_rows), len(align_rows[0][1])))
	samples = list()
	for align_row in align_rows: # Make a list of samples in this locus
		samples.append(str(align_row[0]))		
	i = 0
	allele_labs = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p']
	for align_row in align_rows:
		if samples.count(str(align_row[0])) > ploidy:
			print  "WARNING: Number of alleles > ploidy at locus {0}, sample {0}".format(locus, align_row[0]) 
		if i < (samples.count(str(align_row[0]))-1):		
			out.write("Sample_{0}{1}  {2}\n".format(align_row[0], allele_labs[i], align_row[1]))
			i += 1
		elif i == (samples.count(str(align_row[0]))-1):
			out.write("Sample_{0}{1}  {2}\n".format(align_row[0], allele_labs[i], align_row[1]))
			i = 0
	out.close()			
	align_rows = list()


def nexus(align_rows, locus, ploidy):
	out = open("./out_nexus/locus_{0}.nex".format(locus), 'wb')
	out.write("#NEXUS\n")
	out.write("Begin data;\n")
	out.write("Dimensions ntax={0} nchar={1}\n".format(len(align_rows), len(align_rows[0][1])))
	out.write("Format datatype=dna symbols=\"ACTG\" missing=? gap=-;samples = list()\n")
	out.write("Matrix;\n")
	samples = list()
	for align_row in align_rows: # Make a list of samples in this locus
		samples.append(str(align_row[0]))		
	i = 0
	allele_labs = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p']
	for align_row in align_rows:
		if samples.count(str(align_row[0])) > ploidy:
			print  "WARNING: Number of alleles > ploidy at locus {0}, sample {0}".format(locus, align_row[0]) 
		if i < (samples.count(str(align_row[0]))-1):		
			out.write("Sample_{0}{1}  {2}\n".format(align_row[0], allele_labs[i], align_row[1]))
			i += 1
		elif i == (samples.count(str(align_row[0]))-1):
			out.write("Sample_{0}{1}  {2}\n".format(align_row[0], allele_labs[i], align_row[1]))
			i = 0
	out.write(";\n")
	out.write("End;")
	out.close()			
	align_rows = list()


def main():
	args = get_args()
	infile = open("{0}".format(args.in_file), 'r')
	if args.fasta:
		os.system("mkdir out_fasta")
	if args.phylip:
		os.system("mkdir out_phylip")
	if args.nexus:
		os.system("mkdir out_nexus")
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
							fasta(align_rows, prev_locus, args.ploidy)
						if args.phylip:
							phylip(align_rows, prev_locus, args.ploidy)
						if args.nexus:
							nexus(align_rows, prev_locus, args.ploidy)					
				else: # Output all loci (including those with only one sequence)
					if args.fasta:
						fasta(align_rows, prev_locus, args.ploidy)
					if args.phylip:
						phylip(align_rows, prev_locus, args.ploidy)
					if args.nexus:
						nexus(align_rows, prev_locus, args.ploidy)					
				i += 1
				align_rows = list()
			sample = int(parts[3])
			seq = next(infile).rstrip()
			align_rows.append([sample, seq])
			prev_locus = locus
	if args.fasta:
		fasta(align_rows, prev_locus, args.ploidy)
	if args.phylip:
		phylip(align_rows, prev_locus, args.ploidy)
	if args.nexus:
		nexus(align_rows, prev_locus, args.ploidy)					
	i += 1
	infile.close()
	print "Loci processed: {0}".format(i) # Total number of loci examined

if __name__ == '__main__':
    main()