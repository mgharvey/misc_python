
"""
File: xAndes.py
Author: Michael G. Harvey

Description: Run dadi 5x on fs from GBS data.

Usage: python xAndes_GBS.py 


"""



import os, sys
import numpy
from numpy import array
import dadi
import matplotlib
import matplotlib.pyplot as plt
import demographic_models



for i in range(5):

	dir = ("/scratch/mgharvey/SysBio/dadi/GBS")
	os.chdir(dir)
	outfile = open("./optimized_output_GBS_{0}.txt".format(i+1), 'wb')

	data = dadi.Spectrum.from_file('GBS.fs')
	ns = data.sample_sizes
	pts_l = [30,40,50]


	func = demographic_models.twopop_model
	params = array([1, 1, 1, 1, 1, 1, 1])
	upper_bound = [20, 20, 20, 10, 10, 20, 20]
	lower_bound = [0.01, 0.01, 0.01, 0, 0, 0, 0]


	func_ex = dadi.Numerics.make_extrap_log_func(func)
	model = func_ex(params, ns, pts_l)
	ll_model = dadi.Inference.ll_multinom(model, data)
	theta = dadi.Inference.optimal_sfs_scaling(model, data)


	p0 = dadi.Misc.perturb_params(params, fold=1, upper_bound=upper_bound)
	popt = dadi.Inference.optimize_log_lbfgsb(p0, data, func_ex, pts_l,
				lower_bound=lower_bound,
				upper_bound=upper_bound,
				verbose=len(params), maxiter=20)

	model = func_ex(popt, ns, pts_l)
	ll_opt = dadi.Inference.ll_multinom(model, data)
	optpar = repr(popt)

	outfile.write("Likelihood: {0}\nOptimized Parameters: {1}\nTheta (Scaling Factor): {2}\nConverted thetas: {3} {4} {5} {6} {7}".format(ll_opt, optpar, theta, (popt[0]*theta), (popt[1]*theta), (popt[2]*theta), (popt[3]*theta), (popt[4]*theta)))
	outfile.close()