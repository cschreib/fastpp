# FAST++

This is a C++ version of the popular SED fitting code [FAST](http://w.astro.berkeley.edu/~mariska/FAST.html) (Kriek et al. 2009). Below is a list of the major differences and selling points:

 - FAST++ is free software and does not require an IDL license
 - FAST++ is faster
 - FAST++ can handle larger parameter grids
 - FAST++ can work with FITS tables
 - FAST++ can be used as a drop-in replacement for FAST (same input and output formats)

If you use this code for your own work, please cite this repository and the 2009 paper of M. Kriek where FAST was first introduced.

# Install instructions

To install FAST++, you first need to build it. For this, you will need a recent C++ compiler and [CMake](https://cmake.org/). You also need to have installed the C/C++ libraries [cfitsio](http://heasarc.gsfc.nasa.gov/fitsio/fitsio.html) and [phy++](http://cschreib.github.io/phypp/). Then, download the latest FAST++ release
