#include "nbodyfft.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h> // pulls in declaration of malloc, free
#include <fftw3.h>

//g++ nbody.cpp main.cpp -O3 -o 1dnbody

double cauchy2d(double x1,double x2,  double y1,double y2, double bandx, double bandy ){
	return pow(1.0/(double) (1.0+pow(x1-y1,2) + pow(x2-y2,2)),2);
	//return 1.0/(double) (1.0+pow(x-y,2));
}


int test1d(int ndim) {
//
//	//Initialize some charges and locations
//	int n = 1E4;
//	int nterms = 5;
//	auto n_boxes_per_dim = 600;
//	double xmin = -50;
//	double xmax = 50;
//
//	double * xs =(double*) malloc(n* sizeof(double));
//	double * charges =(double*) malloc(ndim*n* sizeof(double));
//	srand(0);
//	for (int i=0; i< n; i++){
//		xs[i] = (rand() /(double) RAND_MAX - 0.5)*(xmax - xmin);
//		for (int idim = 0; idim<ndim; idim++){
//			charges[i*ndim+idim] = rand() /(double) RAND_MAX;
//		}
//
//	}
//
//    	int nlat =n_boxes_per_dim;
//	int nboxes = nlat*nlat;
//
//        double * band = (double *) calloc(N,sizeof(double));
//	double * boxl =(double*) malloc(2*nboxes* sizeof(double));
//	double * boxr =(double*) malloc(2*nboxes* sizeof(double));
//	double *prods =(double*) malloc(nterms* sizeof(double));
//	double *xpts =(double*) malloc(nterms* sizeof(double));
//	int nfourh = nterms*nlat;
//	double *xptsall =(double*) calloc(nfourh*nfourh, sizeof(double));
//	double *yptsall =(double*) calloc(nfourh*nfourh, sizeof(double));
//	int *irearr =(int*) calloc(nfourh*nfourh, sizeof(int));
//
//	clock_t startTime	=	clock();
//	fftw_complex * zkvalf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*nfourh*2*nfourh);
//	precompute2(maxloc, minloc, maxloc, minloc,nlat, nterms, &cauchy2d, band,boxl, boxr,  prods, xpts, xptsall,yptsall,irearr,zkvalf );
//	clock_t end1 = clock();
//	nbodyfft2(N,ndim,xs, ys,chargesQij, nlat, nterms,boxl, boxr,  prods, xpts, xptsall,yptsall, irearr, zkvalf,potentialQij);
//	clock_t end2 = clock();
//	double time_spent1 = (double)(end1 - begin) / CLOCKS_PER_SEC;
//	double time_spent2 = (double)(end2 - end1) / CLOCKS_PER_SEC;
//	printf("Precompute: %2e\nFast: %.2e seconds, %.2e per second\n", time_spent1, time_spent2, n/time_spent2);
//
//	kernel_type kernel = &squared_cauchy;
//
//	//Brute force compute for at most 100 points
//	int ntest = fmin(100,n);
//	double * truepot =  (double*) calloc(ntest*ndim,sizeof(double));
//	for (int idim=0; idim<ndim;idim++){
//		for (int i = 0; i< ntest; i++){
//			for (int j = 0; j< n; j++){
//				truepot[i*ndim+idim] += charges[j*ndim+idim]*kernel( xs[i],  xs[j]);
//			}
//		}
//	}
//	double diffnorm = 0;
//	double norm = 0;
//	for (int i=0; i<ntest;i++){
//		for (int j=0; j<ndim;j++){
//			norm += fabs(truepot[i*ndim+j]);
//			diffnorm += fabs(truepot[i*ndim+j]- pot[i*ndim+j]);
//			if (i<10) {
//		//	printf("%d,%d:pot  %f, true pot %f \n", i,j, pot[i*ndim+j], truepot[i*ndim+j]);
//			}
//		}
//	}
//	norm /= (ntest*ndim);
//	diffnorm /=(ntest*ndim);
//
//	printf("%.2e L1 norm of difference \n", diffnorm);
//	printf("%.2e L1 norm of difference divided by L1 norm of the true matrix\n", diffnorm/norm);
//
//
//
//
	return 0;

}
int test2d(int ndim) {

	//Initialize some charges and locations
	int N = 1E6;
	int nterms = 3;
	auto nlat = 100;
	double xmin = -50;
	double xmax = 50;
	double ymin = -50;
	double ymax = 50;
        printf("N: %d, nlat: %d, nterms: %d, xmin,ymin: %lf, xmax,ymax:%lf \n", N,nlat, nterms, xmin,xmax);

	double * xs =(double*) malloc(N* sizeof(double));
	double * ys =(double*) malloc(N* sizeof(double));
	double * chargesQij =(double*) malloc(ndim*N* sizeof(double));
	double * potentialQij =(double*) malloc(ndim*N* sizeof(double));



	srand(0);
	for (int i=0; i< N; i++){
		xs[i] = (rand() /(double) RAND_MAX - 0.5)*(xmax - xmin);
		ys[i] = (rand() /(double) RAND_MAX - 0.5)*(ymax - ymin);
		for (int idim = 0; idim<ndim; idim++){
			chargesQij[idim*N+i] = rand() /(double) RAND_MAX;
		}

	}


        double minloc = ymin;
	double maxloc = ymax;


	int nboxes = nlat*nlat;

	double * band = (double *) calloc(N,sizeof(double));
	double * boxl =(double*) malloc(2*nboxes* sizeof(double));
	double * boxr =(double*) malloc(2*nboxes* sizeof(double));
	double *prods =(double*) malloc(nterms* sizeof(double));
	double *xpts =(double*) malloc(nterms* sizeof(double));
	int nfourh = nterms*nlat;
	double *xptsall =(double*) calloc(nfourh*nfourh, sizeof(double));
	double *yptsall =(double*) calloc(nfourh*nfourh, sizeof(double));
	int *irearr =(int*) calloc(nfourh*nfourh, sizeof(int));

        struct timespec start, end, start2, end2;
        clock_gettime(CLOCK_MONOTONIC, &start);
	fftw_complex * zkvalf = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * 2*nfourh*2*nfourh);
	precompute2(maxloc, minloc, maxloc, minloc,nlat, nterms, &cauchy2d, band,boxl, boxr,  prods, xpts, xptsall,yptsall,irearr,zkvalf );
        clock_gettime(CLOCK_MONOTONIC, &end);
	double time_spent1 = (end.tv_nsec - start.tv_nsec)/(double)1E6;
        printf("Precomputation took %.2lf ms\n", time_spent1);
        clock_gettime(CLOCK_MONOTONIC, &start2);
        unsigned int nthreads = 8;
	nbodyfft2(N,ndim,xs, ys,chargesQij, nlat, nterms,boxl, boxr,  prods, xpts, xptsall,yptsall, irearr, zkvalf,potentialQij,nthreads);
        clock_gettime(CLOCK_MONOTONIC, &end2);
	double time_spent2 = (diff(start2,end2))/(double)1E6;
	printf("Fast: %.2lf ms, %.2lf points per ms\n",  time_spent2, N/time_spent2);

	//Brute force compute for at most 100 points
	int ntest = fmin(100,N);
	double * truepot =  (double*) calloc(ntest*ndim,sizeof(double));
	for (int idim=0; idim<ndim;idim++){
		for (int i = 0; i< ntest; i++){
			for (int j = 0; j< N; j++){
				truepot[idim*ntest+i] += chargesQij[idim*N+j]*cauchy2d( xs[i], ys[i], xs[j], ys[j],0,0);
			}
		}
	}
	double diffnorm = 0;
	double norm = 0;
	for (int i=0; i<ntest;i++){
		for (int j=0; j<ndim;j++){
			norm += fabs(truepot[i*ndim+j]);
			diffnorm += fabs(truepot[j*ntest+i]- potentialQij[j+i*ndim]);
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


	//printf("Computing one potential in 1D:\n");
	//test1d(1);
	//printf("Computing five potentials in 1D:\n");
	//test1d(5);

	//printf("Computing one potential in 2D:\n");
	//test2d(1);
	test2d(4);

	return 1;
}
