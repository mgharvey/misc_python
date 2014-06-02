#!/usr/bin/env python

"""

Name: add_phased_snps_to_seqs.py

Author: Michael G. Harvey
Date: 1 June 2014

Using a table of phasing information and consensus sequences from GATK, insert phased heterozygous 
sites and output both haplotype sequences for an individual.

Usage: 	python add_phased_snps_to_seqs.py fasta_file phased_file desired_output_file [ --resolve]

Note: The resolve argument randomly resolves unphased SNPs (forces phasing). Generally not a good 
idea.

"""

import os
import sys
import random
import argparse


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"seq_file",
			type=str,
			help="""The input file of fasta sequences"""
		)
	parser.add_argument(
			"phase_file",
			type=str,
			help="""The input file of phase data"""
		)
	parser.add_argument(
			"out_file",
			type=str,
			help="""The output file of phased sequences"""
		)
	parser.add_argument(
			"--resolve",
			action="""store_true""",
			help="""Randomly resolve unphased SNPs?"""
		)
	return parser.parse_args()


def unphase(alleles): # Finds ambiguity code for each unresolveable heterozygous site	
	a1 = alleles[0]
	a2 = alleles[1]
	if a1 == "A":
		if a2 == "C":
			code = "M"
		elif a2 == "G":
			code = "R"
		elif a2 == "T":
			code = "W"					
		elif a2 == "A":
			code = "A"					
	elif a1 == "C":
		if a2 == "A":
			code = "M"
		elif a2 == "G":
			code = "S"
		elif a2 == "T":
			code = "Y"
		elif a2 == "C":
			code = "C"					
	elif a1 == "G":
		if a2 == "A":
			code = "R"
		elif a2 == "C":
			code = "S"				
		elif a2 == "T":
			code = "K"
		elif a2 == "G":
			code = "G"					
	elif a1 == "T":
		if a2 == "A":
			code = "W"
		elif a2 == "C":
			code = "Y"					
		elif a2 == "G":
			code = "K"	
		elif a2 == "T":
			code = "T"					
	elif a1 == ".":
		code = "N"
	return str(code)


def main():
	args = get_args()
	seqfile = open("{0}".format(args.seq_file))
	phasefile = open("{0}".format(args.phase_file))
	sequences = seqfile.readlines()
	phasings = phasefile.readlines()
	outfile = open("{0}".format(args.out_file), 'wb')
	i = 0
	for sequence in sequences: # For each line in the sequence file
		if sequence.startswith('>'):			
			parts = sequence.split('_')
			firstpart = parts[0].split('>')
			locus = str(firstpart[1]) # Get locus name
			seq = sequences[i+1].lstrip().rstrip()
			seq1list = list(seq) # Set up haplotype 1
			seq2list = list(seq) # Set up haplotype 2
			prev = 'FALSE' 
			for phasing in phasings: # For each line in the phasing file
				parts = phasing.split()
				firstpart = parts[0].split('_')
				phaselocus = firstpart[0] # Get locus name
				if phaselocus == locus: # If locus name from phasefile matches that from seqfile
					parts = phasing.split()
					#snplocs.append(int(parts[1]))
					info = parts[3]
					if '/' in info: # Get the alleles
						alleles = info.split('/')
					elif '|' in info:
						alleles = info.split('|')
					if len(set(alleles)) > 1: # Is the site biallelic for this individual?
						if prev == 'TRUE': # If a biallelic site already exists at this locus
							if '|' in info: # If it was accurately phased
								seq1list[int(parts[1])-1] = alleles[0] # Just insert alleles
								seq2list[int(parts[1])-1] = alleles[1]	
							elif '/' in info: # If it was not phased
								if args.resolve: # Are we randomly resolving?
									rand = random.randint(0,1) # Random number
									if rand == 0:
										other = 1
									elif rand == 1:
										other = 0
									seq1list[int(parts[1])-1] = alleles[rand] # Randomly selected allele
									seq2list[int(parts[1])-1] = alleles[other] # Other allele	
								else: # If not randomly resolving, insert ambiguity code
									base = unphase(alleles)
									seq1list[int(parts[1])-1] = base
									seq2list[int(parts[1])-1] = base			
						elif prev == 'FALSE': # If no biallelic site already exists
							seq1list[int(parts[1])-1] = alleles[0] # Just insert alleles
							seq2list[int(parts[1])-1] = alleles[1]	
						prev = 'TRUE' # Register that a biallelic SNP has been found at this locus
					else: # If the site is monomorphic for this individual				
						seq1list[int(parts[1])-1] = alleles[0].replace('.','N') # Just insert alleles (replace errors with 'N')
						seq2list[int(parts[1])-1] = alleles[1].replace('.','N')							
			name = sequence.rstrip()
			outfile.write("{0}a".format(name))
			outfile.write("\n")
			outfile.write(''.join(seq1list))
			outfile.write("\n")		
			outfile.write("{0}b".format(name))
			outfile.write("\n")
			outfile.write(''.join(seq2list))
			outfile.write("\n")
		i += 1
	seqfile.close()
	phasefile.close()
	outfile.close()
									
if __name__ == '__main__':
    main()