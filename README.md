# BRieszk

BRIESZK ('brisk') is a collection of (Matlab for now) routines that minimize the Riesz energy of N-point configurations, and do so fairly quickly. It was written for personal and internal use at the Vanderbilt University MinEnergy group.
It works on Matlab R2016b; earlier versions may pose compatibility issues for now.

To see details about a specific function, type **help _name_of_the_function_** in Matlab prompt while in the BRieszk folder.

DISCLAIMER: if any of these routines performs fast, it is only due to all the calculations being done in a (possibly, very) sloppy manner. The end result is by no means guaranteed to be a minimizer, even a local one. If you're lucky, it may be _somewhat close_ to a minimizer, whatever 'close' could mean here.

## List of functions:
* riesz_shell(cnf, n, N, r, R, dim, s, plotit, silent)
* riesz_sphere(cnf, N, dim, s, plotit, silent)
* riesz_torus(cnf, N, s, r, R)
* torus_inversion(x, y, z, r, R)
