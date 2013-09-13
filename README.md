INTRODUCTION
============

This source package is a **Finite Elements Method** implementation for the **Finite Strains theory**. 
It also contains a bunch of MATLAB/Octave prototype implementations, auxulary utilities/parsers,
exact solutions for specific cases(uniaxial deformation, uniform pressure of the cylinder -
*Lame problem*).


DIRECTORY STRUCTURE
===================
 * **solver-large** - The C-language solver for finite-strains (*large deformations*) problems with displacements boundary conditions (for now). Material models supported: *Neo-Hookean compressible* material model; *A5 compressible*.
 * **solver-prototype** - a bunch of MATLAB/Octave prototypes for different FEA problems
 * **exact-solutions** - contains exact solutions for the following problems:
   * Uniaxial tension of the block with different material models
   * Uniform pressure of the cylinder (Lame problem)
 * **utilities** - set of file format converters/results analysers

REQUIREMENTS
============
In order to compile the solver the following dependencies needed:
- Logger library: https://github.com/fourier/liblogger
- SEXP parser library: https://github.com/fourier/libsexp
- Sparse Matrix library: https://github.com/fourier/libspmatrix
