# Active-Model-N
This is a Fortran 2003 code for numerically solving the PDEs of 2D Active Model N incorporating nonreciprocal force interactions and torque (see [Active pattern formation emergent from single-species nonreciprocity, arXiv:2404.10093](https://doi.org/10.48550/arXiv.2404.10093))

## **Usage**
The main source code is given in ActiveModelN_2D.F03, while ran3.f90 contains a function of random number generator available in the book "Numerical Recipes". To use and compile the codes, the FFTW3 library is needed, as well as a file fftw3.f03 for the related Fortran 2003 interfaces (which is provided by FFTW; see the FFTW manual at [fftw.org](https://fftw.org/)).

## **Contributing**
This code was developed by [Zhi-Feng Huang](https://s.wayne.edu/huang/).

## **Citation**
If using the code in research, please cite:\
Z.-F. Huang, M. te Vrugt, J. Mayer Martins, R. Wittkowski, and H. Löwen, "Active pattern formation emergent from single-species nonreciprocity", [arXiv:2404.10093](https://doi.org/10.48550/arXiv.2404.10093).
