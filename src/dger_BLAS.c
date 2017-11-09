#include<R.h>
#include<Rmath.h>
#include <R_ext/BLAS.h>

//R CMD SHLIB -o slash.dll -LC:/ARCHIV~1/R/R-30~1.1/bin/i386 -lRblas dger_BLAS.c

void Kronecker_Product(double *, double *, int, int, double *, int, int);
void transpose(double *, int , int );

void dger( double *A, double *x, double *y, int *m, int *n  ){

	//A := alpha*x*y**T + A,
	int inc=1, lda=*m;
	double alpha=1;

	F77_NAME(dger)( m, n, &alpha, x, &inc, y, &inc, A, &lda );

}

void Dmatrix( double *Dt, double *I, int *k, double *T, int *n ){
	int i, j, ind, p=(*n)*(*n+1)/2, nsq=(*n)*(*n);
	double *x = (double *)Calloc( p, double );
	double *y = (double *)Calloc( nsq, double );
	
	for ( j=0; j<*n; j++){
		for( i=j; i<*n; i++ ){
			for( ind=0; ind<p; ind++)	x[ind] = I[  p * (k[ (*n)*j + i ]-1) + ind  ];
			for( ind=0; ind<nsq; ind++)	y[ind] = T[ (j*(*n) + i)*nsq + ind ];
			
			dger( Dt, x, y, &p, &nsq );		
		}	
	}
	
	Free(x);
	Free(y);

}

void Kmatrix( double *K, double *H, int *r, int *c ){
	int p = (*r) * (*c);
	int i, j, ind, inc=1;
	double *Hij = (double *) Calloc( p , double);
	double *tHij = (double *) Calloc( p , double);
	double *HxH = (double *) Calloc( p*p, double );
	

	for( i=0; i<*r; i++){
		for( j=0; j<*c; j++){
			for( ind=0; ind < p; ind++ ) Hij[ind] = H[ (i + j*(*r))*p  + ind ];
			F77_NAME(dcopy)( &p, Hij, &inc, tHij, &inc);
			transpose(tHij, *r, *c);
			Kronecker_Product( HxH, Hij, *r, *c, tHij, *c, *r);
			for( ind=0; ind < p*p ; ind ++ ) K[ind] += HxH[ind];
		}	
	}
	
	Free(Hij);
	Free(tHij);
	Free(HxH);
	
}