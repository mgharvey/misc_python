#!/usr/bin/env python

"""

Name: make_controlfiles_SeqCap.py 

Author: Michael G. Harvey
Date: 13 October 2013

Description: An example Python script for generating control files for sequence capture datasets 
for G-PhoCS.

Usage: python make_controlfiles_SeqCap.py


"""

import os
import sys
import random
import re
from sys import stdout
from time import sleep

outdir = "/Users/michaelharvey/Documents/Research/My_Dissertation/Chapter_2__UCEs_systematics/UCEs_vs_GBS/10.1.13/GPhoCS/SeqCap_control/"
os.chdir(outdir)

for z in range(1000):
	ff = open("SeqCap_{0}.ctl".format(z+1), 'w')
	ff.write("GENERAL-INFO-START\n\n")
	ff.write("\tseq-file\t/scratch/mgharvey/SysBio/gphocs/gphocs_input_SeqCap.txt\n")
	ff.write("\ttrace-file\t/scratch/mgharvey/SysBio/gphocs/SeqCap_control/SeqCap_{0}.log\n".format(z+1))
	ff.write("\tlocus-mut-rate\tCONST\n\n")
	ff.write("\tmcmc-iterations\t1000000\t# Sample iterations after burn-in (Gronau et al. did 200,000)\n")
	ff.write("\tburn-in\t0\t # (Gronau et al. did 100,000)\n")
	ff.write("\titerations-per-log\t10\t# Sample every 10 iterations (Gronau et al. did 10)\n")
	ff.write("\tlogs-per-line\t10\n\n")
	ff.write("\tfind-finetunes\tTRUE\n\n")
	ff.write("\ttau-theta-print\t1000\t# tau and theta values are low, but this should be more than enough (Gronau et al. did 1)\n")
	ff.write("\ttau-theta-alpha\t1.0\t# for STD/mean ratio of 100%\n")
	ff.write("\ttau-theta-beta\t5000\t# mean near expected value\n\n")
	ff.write("\tmig-rate-print\t1.0\n")
	ff.write("\tmig-rate-alpha\t1.0\t# Wide variance (Gronau et al. used wider prior of 0.002, 0.00001)\n")
	ff.write("\tmig-rate-beta\t3\t# mean of 1 \n\n")
	ff.write("GENERAL-INFO-END\n\n")
	ff.write("CURRENT-POPS-START\n\n")
	ff.write("\tPOP-START\n")
	ff.write("\t\tname\tPopA\n")
	ff.write("\t\tsamples\tXm_119_CA d Xm_104_CA d Xm_110_CH d Xm_025_CH d\n")
	ff.write("\tPOP-END\n\n")
	ff.write("\tPOP-START\n")
	ff.write("\t\tname\tPopB\n")
	ff.write("\t\tsamples\tXm_103_NA d Xm_118_NA d Xm_111_SA d Xm_102_SA d\n")
	ff.write("\tPOP-END\n\n")
	ff.write("CURRENT-POPS-END\n\n")
	ff.write("ANCESTRAL-POPS-START\n\n")
	ff.write("\tPOP-START\n")
	ff.write("\t\tname\tPopAB\n")
	ff.write("\t\tchildren\tPopA\tPopB\n")
	ff.write("\t\ttau-initial\t0.0001\n")
	ff.write("\t\ttau-alpha\t1.0\t# Wide variance\n")
	ff.write("\t\ttau-beta\t5000\t# Mean near expected value, similar to priors in Gronau et al.\n")
	ff.write("\tPOP-END\n\n")
	ff.write("ANCESTRAL-POPS-END\n\n")
	ff.write("MIG-BANDS-START\n\n")
	ff.write("\tBAND-START\n")
	ff.write("\t\tsource PopA\n")
	ff.write("\t\ttarget PopB\n")
	ff.write("\tBAND-END\n")
	ff.write("\tBAND-START\n")
	ff.write("\t\tsource PopB\n")
	ff.write("\t\ttarget PopA\n")
	ff.write("\tBAND-END\n")
	ff.write("MIG-BANDS-END\n")					
	ff.close()
