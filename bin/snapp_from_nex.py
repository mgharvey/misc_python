#!/usr/bin/env python

"""

Name: snapp_from_nex.py 

Author: Michael G. Harvey
Date: 12 May 2013

Description: Convert a nexus alignment of concatenated SNPs to the input file format for 
SNAPP (Bryant et al. 2012).

Usage: 	python snapp_from_nex.py input_file output_directory number_of_populations \
			haplotypes_in_pop_1 haplotypes_in_pop_2 ... haplotypes_in_pop_n

Note: Haplotypes must occur in the alignment in the order that you input the number of haplotypes 
per population in the above command. For example, if there are 4 haplotypes (2 diploid individuals) 
in population 1, those should be the top 4 haplotypes in the alignment.

Note 2: The code will prompt you about whether or not you want to include sites with missing data. 
The original version of SNAPP could not handle missing data. In order to use missing data, you need 
a SNAPP add-on. To get it, start BEAUti, click menu File/Manage Add-ons, select SNAPP from the list, 
and click the install button. You may need to restart BEAUti. 

"""

import os
import sys
import argparse
from Bio import AlignIO


def get_args():
	parser = argparse.ArgumentParser(
			description="""Program description""")
	parser.add_argument(
			"in_file",
			type=str,
			help="""The input nexus file"""
		)
	parser.add_argument(
			"out_dir",
			type=str,
			help="""The output directory"""
		)
	parser.add_argument(
			"populations",
			type=int,
			help="""The number of populations"""
		)
	parser.add_argument(
			"pop_sizes",
			type=int,
			nargs='+',
			help="""The number of samples in each population"""
		)
	return parser.parse_args()


def query_yes_no(question, default="yes"):
    valid = {"yes":True,   "y":True,  "ye":True,
             "no":False,     "n":False}
    if default == None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)
    while True:
        sys.stdout.write(question + prompt)
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "\
                             "(or 'y' or 'n').\n")


def reformat_alignment(infile, new_alignment, missing):	
	alignment = AlignIO.read("{0}".format(infile), "nexus")
	for w in xrange(alignment.get_alignment_length()):
		bases = alignment[:, w] 		
		uniqs = list()
		uniqs = list(set(bases))
		if missing == "True": 
			nuniqs = filter(lambda a: a != "N", uniqs) # Filter out N's to check for biallelic-ness
			if len(nuniqs) != 2: # Check for biallelic-ness
				print "Skipping site {0} - not a biallelic SNP".format(w)	
			else:
				allele1 = nuniqs[0]
				allele2 = nuniqs[1]
				new_snp = [0]*(len(bases))				
				for x in range(len(bases)):
					if bases[x] == allele1:
						new_snp[x] = 0
					elif bases[x] == allele2:
						new_snp[x] = 1
					elif bases[x] == "N":
						new_snp[x] = "?"			
				new_alignment.append(new_snp)
		elif missing == "False":
			if "N" in uniqs:
				print "Skipping site {0} - contains missing data".format(w)	 
			else:
				nuniqs = filter(lambda a: a != "N", uniqs) # Filter out N's to check for biallelic-ness
				if len(nuniqs) != 2: # Check for biallelic-ness
					print "Skipping site {0} - not a biallelic SNP".format(w)	
				else:
					allele1 = nuniqs[0]
					allele2 = nuniqs[1]
					new_snp = [0]*(len(bases))				
					for x in range(len(bases)-1):
						if bases[x] == allele1:
							new_snp[x] = 0
						elif bases[x] == allele2:
							new_snp[x] = 1
						elif bases[x] == "N":
							new_snp[x] = "?"
					new_alignment.append(new_snp)
	return new_alignment


def make_names(populations, pop_sizes):
	names = list()
	i = 0
	j = 0
	for y in range(populations): # For each pop
		for z in range(pop_sizes[i]): # For each individual in that pop
			names.append("{0}_{1}\t".format(i+1,j+1)) # Write name
			j += 1 # Plus one individual
		i += 1 # Plus one pop
	return names


def write_outfile(infile, outdir, new_alignment, names):
	alignment = AlignIO.read("{0}".format(infile), "nexus")
	out = open("{0}snapp_input.nex".format(outdir), 'wb')
	out.write("#NEXUS\n\n")
	out.write("Begin data;\n")
	out.write("\tDimensions ntax={0} nchar={1};\n".format(len(alignment[:,1]), len(new_alignment)))
	out.write("\tFormat datatype=binary symbols=\"01\" gap=- missing=?;\n")
	out.write("\tMatrix\n")
	for x in range(len(alignment[:,1])): # Loop through individuals
		out.write(names[x])				
		for y in range(len(new_alignment)): # Loop over alignment columns
			out.write("{0}".format(new_alignment[y][x])) # Write nucleotides
			out.flush()
		out.write("\n")	
		out.flush()
	out.write("\t;\n")
	out.write("End;")
	out.close()
	
	
def main():
	args = get_args()
	missing = str(query_yes_no("\nInclude sites with missing data?: "))
	new_align = list()
	write_outfile(args.in_file, args.out_dir, reformat_alignment(args.in_file, new_align, missing), make_names(args.populations, args.pop_sizes))

if __name__ == '__main__':
    main()
