////////////////////////////////////////////////////////////////////////////////
// File: kronecker_product.c                                                  //
// Routines:                                                                  //
//    Kronecker_Product                                                       //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Kronecker_Product(double *C, double *A, int nrows, int ncols,        //
//                                           double *B, int mrows, int mcols) //
//                                                                            //
//  Description:                                                              //
//     This routine initializes the (nrows x mrows) x (ncols x mcols) matrix  //
//     C as the Kronecker product of the nrows x ncols matrix A and the       //
//     mrows x mcols matrix B.                                                //
//                                                                            //
//     Def:  The Kronecker product of the n x m matrix A and the n' x m'      //
//           matrix B is the n*n' x m*m' matrix C where C in block matrix     //
//           form is given by                                                 //
//                 -                                              -           //
//                 |  A[0][0] B     A[0][1] B   ...   A[0][m-1] B |           //
//            C =  |  A[1][0] B     A[1][1] B   ...   A[1][m-1] B |           //
//                 | ............................................ |           //
//                 | A[n-1][0] B   A[n-1][1] B  ...  A[n-1][m-1] B|           //
//                 -                                              -           //
//                                                                            //
//     The number of rows of C is the product of the number of rows of A and  //
//     the number of rows of B, similarly, the number of columns of C is the  //
//     product of the number of columns of A and the number of columns of B.  //
//                                                                            //
//     Def: Given vector spaces V and W, the tensor product of V and W is the //
//     vector space V (x) W and bilinear map i: VxW -> V (x) W,               //
//     i(v,w)=v (x) w (unique up to isomorphism) such that given              //
//     a third vector space U and a bilinear map h:VxW->U there exists a      //
//     unique linear map k: (V (x) W) -> U such that h = k i.                 //
//        Given bases {e1[i]} of V1, {e2[i]} of V2, {f1[i]} of W1 and {f2[i]} //
//     of W2, then if A[i][j] is the matrix of A with respect to the bases    //
//     {e1[i]} and {f1[i]} and B[i][j] is the matrix of B with respect to the //
//     bases {e2[i]} and {f2[i]}, then the Kronecker product of the matrices  //
//     A and B, is the matrix of the linear map (A (x) B) defined by the      //
//     equation                                                               //
//              ( A (x) B)( v1 (x) v2 ) =  (Av1) (x) (Bv2)                    //
//     with respect to the basis vectors {e1[i] (x) e2[j]} and                //
//     {f1[i] (x) f2[j]} ordered lexicographically, i.e. <e1[0] (x) e2[0],    //
//     e1[0] (x) e2[1], ... , e1[0] (x) e2[m'-1], e1[1] (x) e2[0], ... ,      //
//     e1[m-1] (x) e2[m'-1]> and <f1[0] (x) f2[0], ..., f1[0] (x) f2[n'-1],   //
//     ..., f1[n-1] (x) f2[n'-1]>.                                            //
//                                                                            //
//  Arguments:                                                                //
//     double *C       Pointer to the first element of the matrix C, the      //
//                     Kronecker product of the two matrix arguments in the   //
//                     order they appear. The matrix C should be declared     //
//                     double C[nrows*mrows][ncols*mcols] in the calling      //
//                     routine.                                               //
//     double *A       Pointer to the first matrix in the Kronecker product.  //
//     int    nrows    The number of rows of A.                               //
//     int    ncols    The number of cols of A.                               //
//     double *B       Pointer to the second matrix in the Kronecker product. //
//     int    mrows    The number of rows of B.                               //
//     int    mcols    The number of cols of B.                               //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N1                                                             //
//     #define M1                                                             //
//     #define N2                                                             //
//     #define M2                                                             //
//                                                                            //
//                                                                            //
//     double A[N1][M1], B[N2][M2], C[N1*N2][M1*M2];                          //
//                                                                            //
//     (your code to initialize the matrices A and B)                         //
//                                                                            //
//     Kronecker_Product(&C[0][0], (double *)A, N1, M1, (double*)B, N2, M2);  //
//     printf(" The Kronecker product C of A and B is  \n");                  //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////

void Kronecker_Product(double *C, double *A, int nrows, int ncols, 
                                               double *B, int mrows, int mcols)
{
   int ccols, i, j, k, l;
   int block_increment;
   double *pB;
   double *pC, *p_C;
 
   ccols = ncols * mcols;
   block_increment = mrows * ccols;
   for (i = 0; i < nrows; C += block_increment, i++)
      for (p_C = C, j = 0; j < ncols; p_C += mcols, A++, j++) 
         for (pC = p_C, pB = B, k = 0; k < mrows; pC += ccols, k++)
            for (l = 0; l < mcols; pB++, l++) *(pC+l) = *A * *pB; 

}

void transpose(double *m, int w, int h)
{
	int start, next, i;
	double tmp;
 
	for (start = 0; start <= w * h - 1; start++) {
		next = start;
		i = 0;
		do {	i++;
			next = (next % h) * w + next / h;
		} while (next > start);
		if (next < start || i == 1) continue;
 
		tmp = m[next = start];
		do {
			i = (next % h) * w + next / h;
			m[next] = (i == start) ? tmp : m[i];
			next = i;
		} while (next > start);
	}
}
