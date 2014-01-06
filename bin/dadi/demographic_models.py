#!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python

"""

File: demographic_models.py
Author: Michael G. Harvey

Description: models used for analyses of sequence capture and GBS data in dadi


"""

import numpy
import dadi

def twopop_model((nu1, nu2, nu3, T1, T2, m1, m2), (n1,n2), pts):

    # Define the grid we'll use
    xx = yy = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    # Integration with one pop
    phi = dadi.Integration.one_pop(phi, xx, T1, nu1)

    # The first divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # Integration with two pops
    phi = dadi.Integration.two_pops(phi, xx, T2, nu2, nu3, m12=m1, m21=m2)
                                    
    # Finally, calculate the spectrum.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return sfs

def twopop_nomig_model((nu1, nu2, nu3, T1, T2), (n1,n2), pts):

    # Define the grid we'll use
    xx = yy = dadi.Numerics.default_grid(pts)

    # phi for the equilibrium ancestral population
    phi = dadi.PhiManip.phi_1D(xx)
    
    # Integration with one pop
    phi = dadi.Integration.one_pop(phi, xx, T1, nu1)

    # The first divergence
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)
    
    # Integration with two pops
    phi = dadi.Integration.two_pops(phi, xx, T2, nu2, nu3, m12=0, m21=0)
                                    
    # Finally, calculate the spectrum.
    sfs = dadi.Spectrum.from_phi(phi, (n1,n2), (xx,xx))
    return sfs
