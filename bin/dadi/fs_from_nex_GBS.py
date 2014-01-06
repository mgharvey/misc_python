#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

"""
File: fs_from_nex_GBS.py
Author: Michael G. Harvey

Created by Michael G. Harvey on 21 September 2012
Copyright (c) 2012 Michael G. Harvey. All rights reserved.

Description: generate multi-population allele frequency spectra from SNPs from GBS dataset

Usage: python fs_from_nex_GBS.py

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
    
    fs = dadi.Spectrum.from_data_dict(dd, ['p1', 'p2'], [8,8], polarized=False)
    
    # Plot the fs
    
    import pylab
    dadi.Plotting.plot_single_2d_sfs(fs, vmin=0.1)
    pylab.show()
    
    # Print fs to file
    
    fs.tofile('GBS.fs')
    
    
