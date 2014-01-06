#!/usr/bin/env python

"""

Name: fill_in_missing_taxa.py 

Author: Michael G. Harvey
Date: 24 July 2013

Description: Take directory of Nexus files containing some alignments that are missing some taxa, 
and fill in missing taxa (with missing data for entire sequence) before printing out new Nexus 
files.

Usage: python fill_in_missing_taxa.py in_dir out_dir


"""

import os
import sys
import argparse
from Bio import AlignIO


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
			help="""The output directory"""
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
	samples = list()
	for file in files:
		alignment = AlignIO.read("{0}{1}".format(args.in_dir, file), "nexus")
		for record in alignment:
			name = record.id
			if name not in samples:
				samples.append(name)
	print "Samples detected across alignments:\n{0}".format(samples)
	for file in files:
		alignment = AlignIO.read("{0}{1}".format(args.in_dir, file), "nexus")
		filenameparts = file.split('.')
		filename = filenameparts[0]
		out = open("{0}{1}.nex".format(args.out_dir, filename), 'wb')
		out.write("#NEXUS\n\n")
		out.write("Begin data;\n")
		out.write("\tDimensions ntax={0} nchar={1};\n".format(len(alignment), alignment.get_alignment_length()))
		out.write("\tFormat datatype=DNA gap=- missing=?;\n")
		out.write("\tMatrix\n")
		recs = list()
		for record in alignment:
			if record.id not in recs:
				recs.append(record.id)
		for sample in samples:
			if sample in recs:
				for record in alignment:
					if sample == record.id:
						out.write("{0}\t".format(record.id))
						out.write("{0}\n".format(record.seq))
			elif sample not in recs:
				out.write("{0}\t".format(sample))
				for w in xrange(alignment.get_alignment_length()):
					out.write("{0}".format("n"))
				out.write("\n")
			out.flush()
		out.write("\t;\n")
		out.write("End;")
		out.close()

if __name__ == '__main__':
    main()
