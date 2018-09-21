# cauchy-kernel-interpolation-accuracy

Code for checking accuracy of interpolation scheme underlying FIt-SNE.

Compile:


`g++ nbodyfft.cpp main.cpp -O3 -o nbody  -lfftw3 -lm -lfftw3_omp -pthread  -std=c++11` -fopenmp
