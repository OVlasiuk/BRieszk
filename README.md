# BRieszk
<p align="center">
<img src="https://raw.githubusercontent.com/OVlasiuk/BRieszk/master/img/sph_paisley.png" width="800">
</p>

---
BRIESZK ('brisk') is a collection of routines for (currently, only) Matlab that allow to approximate gradient flows of the Riesz energy in several contexts, and to work with the output of such flows. It was written for personal and internal use at the Vanderbilt University MinEnergy group.<br>
BRieszk works on Matlab R2016b+; earlier versions may pose compatibility issues.

Apart from the optimization routines, there are also functions to produce nearest neighbor statistics for a given node set. To see details about a specific function, type **help _name_of_the_function_** in Matlab prompt while in the BRieszk folder. 

DISCLAIMER: if any of these routines performs fast, it is only due to all the calculations being done in a (possibly, very) sloppy manner. The end result is by no means guaranteed to be a full-energy minimizer, even a local one. If you're lucky, it may be _somewhat close_ to a minimizer, whatever 'close' could mean here.

## Function/file listing:

* f_analyzer.m -- Nearest neighbor statistics and domain inclusion.
* f_cnfinit.m -- An example of producing a random configuration on an implicit surface.
* f_dcompare.m -- Compares node separation distances and values of density at these nodes.
* f_repel.m -- Compute a naive repulsive flow using the Riesz kernel.
* f_torinv.m -- Transform Cartesian ambient coordinates to the intrinsic angles on a 2-torus.
* f_vorsurf.m -- Produce a Voronoi tessellation of a surface.
* riesz_shell.m -- Uniform node distribution in a spherical shell.
* riesz_sphere.m -- ... on a sphere.
* riesz_surf.m -- ... on an implicit surface.
* riesz_torus.m -- ... on (the surface of) a 2-torus.

