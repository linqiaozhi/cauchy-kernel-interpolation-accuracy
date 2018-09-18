# cauchy-kernel-interpolation-accuracy

Code for checking accuracy of interpolation scheme underlying FIt-SNE. The number of threads is hardcoded as `nthreads`. Currently, only 2D code works.

Compile:


`g++ nbodyfft.cpp main.cpp -O3 -o nbody  -lfftw3 -lm  -std=c++11 -fopenmp`

or 

`g++ nbodyfft.cpp main.cpp -O3 -o nbody  -lfftw3 -lm  -std=c++11 -pthread`

Depending on whether you want to use OpeMP or C++11 threading
