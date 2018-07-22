#include "nbodyfft.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h> // pulls in declaration of malloc, free
#include <fftw3.h>

//g++ nbody.cpp main.cpp -O3 -o 1dnbody
double squared_cauchy(double x, double y) {
	return pow(1.0 + pow(x - y, 2), -2);
}


double squared_cauchy_2d(double x1, double x2, double y1, double y2) {
	return pow(1.0 + pow(x1 - y1, 2) + pow(x2 - y2, 2), -2);
}
int test2d(int ndim) {

	//Initialize some charges and locations
	int n = 1E4;
	int n_interpolation_points = 5;
	auto n_boxes_per_dim = 600;
	double xmin = -50;
	double xmax = 50;
	double ymin = -50;
	double ymax = 50;

	double * xs =(double*) malloc(n* sizeof(double));
	double * ys =(double*) malloc(n* sizeof(double));
	double * charges =(double*) malloc(ndim*n* sizeof(double));

	srand(0);
	for (int i=0; i< n; i++){
		xs[i] = (rand() /(double) RAND_MAX - 0.5)*(xmax - xmin);
		ys[i] = (rand() /(double) RAND_MAX - 0.5)*(ymax - ymin);
		for (int idim = 0; idim<ndim; idim++){
			charges[i*ndim+idim] = rand() /(double) RAND_MAX;
		}

	}

	int n_boxes = n_boxes_per_dim * n_boxes_per_dim;

	auto *box_lower_bounds = new double[2 * n_boxes];
	auto *box_upper_bounds = new double[2 * n_boxes];
	auto *y_tilde_spacings = new double[n_interpolation_points];
	int n_interpolation_points_1d = n_interpolation_points * n_boxes_per_dim;
	auto *x_tilde = new double[n_interpolation_points_1d]();
	auto *y_tilde = new double[n_interpolation_points_1d]();
	auto *fft_kernel_tilde = new complex<double>[2 * n_interpolation_points_1d * 2 * n_interpolation_points_1d];

	auto *pot = new double[n * ndim]();
	clock_t begin = clock();
	precompute_2d(xmax, xmin, ymax, ymin, n_boxes_per_dim, n_interpolation_points,
			&squared_cauchy_2d,
			box_lower_bounds, box_upper_bounds, y_tilde_spacings, x_tilde, y_tilde, fft_kernel_tilde);
	clock_t end1 = clock();
	n_body_fft_2d(n, ndim, xs, ys, charges, n_boxes_per_dim, n_interpolation_points, box_lower_bounds,
			box_upper_bounds, y_tilde_spacings, fft_kernel_tilde, pot);
	clock_t end2 = clock();
	double time_spent1 = (double)(end1 - begin) / CLOCKS_PER_SEC;
	double time_spent2 = (double)(end2 - end1) / CLOCKS_PER_SEC;
	printf("Precompute: %2e\nFast: %.2e seconds, %.2e per second\n", time_spent1, time_spent2, n/time_spent2);

	kernel_type_2d kernel = &squared_cauchy_2d;

	//Brute force compute for at most 100 points
	int ntest = fmin(100,n);
	double * truepot =  (double*) calloc(ntest*ndim,sizeof(double));
	for (int idim=0; idim<ndim;idim++){
		for (int i = 0; i< ntest; i++){
			for (int j = 0; j< n; j++){
				truepot[i*ndim+idim] += charges[j*ndim+idim]*kernel( xs[i], ys[i], xs[j], ys[j]);
			}
		}
	}
	double diffnorm = 0;
	double norm = 0;
	for (int i=0; i<ntest;i++){
		for (int j=0; j<ndim;j++){
			norm += fabs(truepot[i*ndim+j]);
			diffnorm += fabs(truepot[i*ndim+j]- pot[i*ndim+j]);
			if (i<10) {
		//	printf("%d,%d:pot  %f, true pot %f \n", i,j, pot[i*ndim+j], truepot[i*ndim+j]);
			}
		}
	}
	norm /= (ntest*ndim);
	diffnorm /=(ntest*ndim);

	printf("%.2e L1 norm of difference \n", diffnorm);
	printf("%.2e L1 norm of difference divided by L1 norm of the true matrix\n", diffnorm/norm);




	return 0;

}

int main() {
	printf("Computing one potential:\n");
	test2d(1);
	printf("Computing five potentials:\n");
	test2d(5);

	return 1;
}
