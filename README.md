This package finds the eigenvalue and eigenfunctions of a one dimensional time independent Schrodinger equation (Sturm-Liouville problem) using the Numerov algorithm. The core is implemented in C++ and through Rcpp can be called from inside R.

### Example
The workflow is like this:

```r
library(numerov)
setPotential(px, py)
computeSpectrum(10)
getEnergies()
getWavefunctions()
```

The arguments of the function *setPotential()* are the *x* and *y* coordinates of the potential. The numerical boundaries used by Numerov are deduced from this (so be careful, if you are finding the first 100 eigenfunctions of the harmonic oscillator make sure you provide a potential range such that indeed the higher order eigenfunctions are small enough at its ends).

*computeSpectrum(N)* computes the first N eigenvalues and eigenfunctions. The routine starts by looking at the global minimum of the potential and uses it as starting point to find the zeros of an internal spectral curve definition. The scanning is controlled by the second parameter of this function, which should be always smaller than the distance between two eigenvalues.

*getEnergies()* and *getWavefunctions()* does the obvious thing.

You can use Numerov straight in your C++ code, just include the *numerov.hpp* file in your project.

