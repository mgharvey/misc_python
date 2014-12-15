"""

Name: process_stacks_alignments.py

Author: Michael G. Harvey
Date: 30 July 2014

Description: To be used following processing with alignments_from_stacks_fasta.py. 
Collapses alleles into single sequences with ambiguity codes for heterozygous sites, 
removes alignments with only one individual, and relabels samples using a sample map.

Usage:

python process_stacks_alignments.py in_dir sample_map out_dir

The sample map is a separate file that looks like this:

Sample_1:Rallus_crepitans_63400
Sample_2:Rallus_crepitans_63404
Sample_3:Rallus_crepitans_63475
Sample_4:Rallus_crepitans_63477
Sample_5:Rallus_crepitans_63464
Sample_6:Rallus_crepitans_63467
Sample_7:Rallus_crepitans_63471
Sample_8:Rallus_crepitans_63473

The sample numbers correspond to the Sample ID's assigned by Stacks. You can find these in
the .tsv files (e.g. see the second column in the .tags.tsv files) for each individual 
that are output by Stacks.

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
		print file
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
					bases_b = bases
					new_align[sample_key[name]] = [bases_a, bases_b]
				else:
					bases_a = bases
					new_align[sample_key[name]] = [bases_a, bases_a]
				names.append(name)
		if len(new_align) > 1: # Only print loci with >1 individual
			out = open("{0}/{1}".format(args.out_dir, file), 'wb')
			out.write("#NEXUS\n\n")
			out.write("Begin data;\n")
			out.write("\tDimensions ntax={0} nchar={1};\n".format((2*len(new_align)), length))
			out.write("\tFormat datatype=dna gap=-;\n")
			out.write("\tMatrix\n")
			for key in new_align:
				out.write("{0}a\t{1}\n".format(key, new_align[key][0]))
				out.write("{0}b\t{1}\n".format(key, new_align[key][1]))
			out.write(";\n")
			out.write("End;\n")
			out.close()
		infile.close()
	map_file.close()
				
if __name__ == '__main__':
    main()