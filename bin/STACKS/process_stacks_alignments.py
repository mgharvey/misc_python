"""

Name: process_stacks_alignments.py

Author: Michael G. Harvey
Date: 30 July 2014

Description: To be used following processing with alignments_from_stacks_fasta.py. 
Collapses alleles into single sequences with ambiguity codes for heterozygous sites, 
removes alignments with only one individual, and relabels sampling using a sample map.

Usage:

python process_stacks_alignments.py in_dir sample_map out_dir

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
			help="""The input directory of alignments"""
		)
	parser.add_argument(
			"sample_map",
			type=str,
			help="""The sample map"""
		)
	parser.add_argument(
			"out_dir",
			type=str,
			help="""The output directory"""
		)
	return parser.parse_args()


def main():
	args = get_args()
	files = list()
	prefiles = os.listdir(args.in_dir)
	for prefile in prefiles:
		if not prefile.startswith('.'):
			files.append(prefile)
	map_file = open(args.sample_map, 'r')
	stacks_samples = list()
	my_samples = list()
	for line in map_file:
		parts = line.split(':')
		stacks_samples.append(parts[0])
		my_samples.append(parts[1].rstrip())
	sample_key = dict(zip(stacks_samples, my_samples))
	for file in files:
		infile = open("{0}/{1}".format(args.in_dir, file))
		new_align = {}
		names = list()
		for line in infile:
			if line.startswith(" "):
				parts = line.split()
				length = parts[1]
			elif line.startswith("Sample"):
				parts = line.split()
				prename = parts[0]
				name = prename[:-1]
				bases = parts[1]
				if name in names:
					new_bases = list()
					previous_bases = new_align[sample_key[name]]
					for i, base in enumerate(bases):
						if base == previous_bases[i]:
							new_base = base
						elif base in ["A","C"] and previous_bases[i] in ["A","C"]:
							new_base = "M"
						elif base in ["A","G"] and previous_bases[i] in ["A","G"]:
							new_base = "R"
						elif base in ["A","T"] and previous_bases[i] in ["A","T"]:
							new_base = "W"
						elif base in ["C","G"] and previous_bases[i] in ["C","G"]:
							new_base = "S"
						elif base in ["C","T"] and previous_bases[i] in ["C","T"]:
							new_base = "Y"
						elif base in ["G","T"] and previous_bases[i] in ["G","T"]:
							new_base = "K"
						else:
							print "Error!: {0} or {1} not a recognized nucleotide".format(base, previous_bases[i])
						new_bases.append(new_base)
					new_align[sample_key[name]] = ''.join(new_bases)						
				else:
					new_align[sample_key[name]] = bases
				names.append(name)
		if len(new_align) > 1: # Only print loci with >1 individual
			out = open("{0}/{1}".format(args.out_dir, file), 'wb')
			out.write(" {0} {1}\n".format(len(new_align), length))
			for key in new_align:
				out.write("{0}  {1}\n".format(key, new_align[key]))
		infile.close()
				
if __name__ == '__main__':
    main()