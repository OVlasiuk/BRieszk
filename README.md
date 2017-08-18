# BRieszk
<p align="center">
<img src="https://raw.githubusercontent.com/OVlasiuk/BRieszk/master/sph_paisley.png" width="400">
</p>

---
BRIESZK ('brisk') is a collection of (Matlab for now) routines that minimize the Riesz energy of N-point configurations, and do so fairly quickly, using nearest neighbor truncation. It was written for personal and internal use at the Vanderbilt University MinEnergy group.
It works on Matlab R2016b+; earlier versions may pose compatibility issues.

Apart from the optimization routines, there are also functions to produce nearest neighbor statistics up to several nearest neighbors for a given node set. To see details about a specific function, type **help _name_of_the_function_** in Matlab prompt while in the BRieszk folder. 

DISCLAIMER: if any of these routines performs fast, it is only due to all the calculations being done in a (possibly, very) sloppy manner. The end result is by no means guaranteed to be a full-energy minimizer, even a local one. If you're lucky, it may be _somewhat close_ to a minimizer, whatever 'close' could mean here.

## List of functions:
* dcompare(cnf, densityF, plotit) <br>
Compares node separation distances and values of density at these nodes.
* pt_analyzer(cnf, in_domainF) <br>
Nearest neighbor statistics and domain inclusion.
* riesz_shell(cnf, n, N, r, R, dim, s, plotit, analyzeit, silent) <br>
* riesz_sphere(cnf, N, dim, s, plotit, silent)
* riesz_torus(cnf, N, s, r, R)
* torus_inversion(x, y, z, r, R) <br>
Transforms Cartesian ambient coordinates to the intrinsic angles.
