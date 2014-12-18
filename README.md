
#markuspline

Fortran module which contains many useful tools to caculate and manage splines as in [1].
Please refer to the documentation for more informations.

[1] Markus Holzmann and Bernard Bernu.
    Optimized periodic 1/r Coulomb potential in two dimensions.
    Journal of Computational Physics, 206(1):111 – 121, 2005.
    
#Build a library

Use the following instructions:

   `gfortran -O3 -march=native -c module_markuspline.f90`

   `ar rcv libmarkuspline.a *.o`
   
   `ranlib libmarkuspline.a`
