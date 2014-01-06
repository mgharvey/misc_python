#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

"""
File: fs_from_nex_SeqCap.py
Author: Michael G. Harvey

Created by Michael G. Harvey on 21 September 2012
Copyright (c) 2012 Michael G. Harvey. All rights reserved.

Description: generate multi-population allele frequency spectra from SNPs from sequence capture
dataset

Usage: python fs_from_nex_SeqCap.py

NOTE: This script depends on a modified version of the custom script "NeXus.py" written by 
R. Gutenkunst specifically for use with alignments in nexus format.

"""

import dadi
import NeXus_mod


"""
Generate frequency spectrum for 2 populations for xAndes.

"""    
    
    
if __name__ == '__main__':
    pop_assignments = {'Xenopsminutus025_CHa':'p1',
                       'Xenopsminutus025_CHb':'p1',
                       'Xenopsminutus102_SAa':'p2',
                       'Xenopsminutus102_SAb':'p2',
                       'Xenopsminutus103_NAa':'p2',
                       'Xenopsminutus103_NAb':'p2',
                       'Xenopsminutus104_CAa':'p1',
                       'Xenopsminutus104_CAb':'p1',
                       'Xenopsminutus110_CHa':'p1',
                       'Xenopsminutus110_CHb':'p1',
                       'Xenopsminutus111_SAa':'p2',
                       'Xenopsminutus111_SAb':'p2',
                       'Xenopsminutus118_NAa':'p2',
                       'Xenopsminutus118_NAb':'p2',
                       'Xenopsminutus119_CAa':'p1',
                       'Xenopsminutus119_CAb':'p1'
                       }


    # Generate data dictionary from NeXus alignment

    dd = NeXus_mod.data_dict_from_file('SeqCap_SNPs.nex', pop_assignments)
    
    # Generate frequency spectrum from data dictionary
    
    fs = dadi.Spectrum.from_data_dict(dd, ['p1', 'p2'], [8,8], polarized=False)
    
    # Plot the fs
    
    import pylab
    dadi.Plotting.plot_single_2d_sfs(fs, vmin=0.1)
    pylab.show()
    
    # Print fs to file
    
    fs.tofile('SeqCap.fs')
    
    
