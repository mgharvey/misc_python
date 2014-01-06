#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

"""
File: summary_stats_GBS.py
Author: Michael G. Harvey

Created by Michael G. Harvey on 21 September 2012
Copyright (c) 2012 Michael G. Harvey. All rights reserved.

Description: generate multi-population allele frequency spectra from SNPs from GBS dataset and 
estimates summary statistics using dadi

Usage: python summary_stats_GBS.py

NOTE: This script depends on a modified version of the custom script "NeXus.py" written by 
R. Gutenkunst specifically for use with alignments in nexus format.

"""

import dadi
import NeXus_mod


"""
Generate frequency spectrum for 2 populations for xAndes.

"""    
    
    
if __name__ == '__main__':
    pop_assignments = {'XM11a':'p1',
                       'XM11b':'p1',
                       'XM13a':'p1',
                       'XM13b':'p1',
                       'XM2a':'p1',
                       'XM2b':'p1',
                       'XM6a':'p1',
                       'XM6b':'p1',
                       'XM30a':'p2',
                       'XM30b':'p2',
                       'XM31a':'p2',
                       'XM31b':'p2',
                       'XM33a':'p2',
                       'XM33b':'p2',
                       'XM40a':'p2',
                       'XM40b':'p2'
                       }


    # Generate data dictionary from NeXus alignment

    dd = NeXus_mod.data_dict_from_file('GBS_SNPs_split.nex', pop_assignments)
    
    # Generate frequency spectrum from data dictionary
    
    fs1 = dadi.Spectrum.from_data_dict(dd, ['p1'], [8], polarized=False)
    fs2 = dadi.Spectrum.from_data_dict(dd, ['p2'], [8], polarized=False)
    
    print "Population 1"
    theta1 = fs1.Watterson_theta()
    print theta1
    pi1 = fs1.pi()
    print pi1
    D1 = fs1.Tajima_D()
    print D1

    print "Population 2"
    theta1 = fs2.Watterson_theta()
    print theta1
    pi1 = fs2.pi()
    print pi1
    D1 = fs2.Tajima_D()
    print D1
