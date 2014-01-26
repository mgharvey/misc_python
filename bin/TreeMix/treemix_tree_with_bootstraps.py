#!/usr/bin/env python

"""

Name: treemix_tree_with_bootreps.py 

Author: Michael G. Harvey
Date: 13 May 2013

Dependencies:

treemix (https://code.google.com/p/treemix/)
sumtrees package in dendropy (http://pythonhosted.org/DendroPy/scripts/sumtrees.html)

Usage: 	python treemix_tree_with_bootstraps.py input_file output_directory outgroup_name
			number_of_bootstraps_to_run size_of_bootstrap_snp_blocks

Make sure Treemix is in your path. If, when installing treemix, there is an issue locating boost, 
try using CPPFLAGS to specify the location of the boost libraries (replacing "PATH_TO_BOOST" with 
the correct location):

>./configure CPPFLAGS=-I/PATH_TO_BOOST/boost_1_53_0/

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
			"outgroup",
			type=str,
			help="""The name of the outgroup"""
		)
	parser.add_argument(
			"bootreps",
			type=int,
			help="""The number of bootstrap replicates to run"""
		)
	parser.add_argument(
			"block_size",
			type=int,
			help="""The number of SNPs to include in bootstrap blocks"""
		)	
	return parser.parse_args()

def run_bootreps(in_file, out_dir, bootreps, block_size, outgroup):
	infile = "{0}.gz".format(in_file) # SNP data in treemix format
	reps2 = bootreps
	for i in range(reps2):
		print "\nBootRep: {0}\n".format(i+1)
		outfile = "{0}bootstraps/treemix_bootrep_{1}".format(out_dir, i+1)
		os.system("treemix -i {0} -bootstrap -k {1} -root {2} -o {3}".format(infile, block_size, outgroup, outfile))
	
def combine_bootreps(out_dir, bootreps):
	outfile = open("{0}cat_trees.tre".format(out_dir), 'wb')
	for i in range(bootreps):
		infile = open("{0}bootstraps/treemix_bootrep_{1}.treeout.gz".format(out_dir, i+1), 'r')
		for line in infile:
			outfile.write(line)
		infile.close()
	outfile.close()

def main():
	args = get_args()
	os.system("gzip {0}".format(args.in_file)) # gzip input file
	os.system("treemix -i {0}.gz -root {1} -o {2}out_stem".format(args.in_file, args.outgroup, args.out_dir)) # build ML tree
	os.system("mkdir {0}bootstraps".format(args.out_dir)) # make dir for bootreps
	run_bootreps(args.in_file, args.out_dir, args.bootreps, args.block_size, args.outgroup) # run bootreps
	combine_bootreps(args.out_dir, args.bootreps) # combine bootreps into one file
	# os.system("sumtrees.py --rooted -t {0}out_stem.treeout.tre -o {1}boottree.txt {2}cat_trees.tre".format(args.out_dir, args.out_dir, args.out_dir))
	# sumtrees.py --rooted -t ./out_stem.treeout.tre -o ./boottree.txt ./cat_trees.tre

"""
I had problems with file encoding, you may need to paste the trees from the cat_trees.tre and 
out_stem.treeout.tre files into existing .tre files and save as to get them to format correctly.

The resulting tree is in NEXUS by default, and can be visualized in e.g. FigTree. You can also 
visualize the tree in R using the commands (replacing "PATH_TO_TREEMIX" with the correct path):

>source("/PATH_TO_TREEMIX/treemix-1.11/src/plotting_funcs.R")
>plot_tree("/PATH_TO_TREEMIX/treemix-1.11/out_stem")

"""


if __name__ == '__main__':
    main()
