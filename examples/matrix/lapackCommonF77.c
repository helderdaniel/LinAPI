/* ATLAS F77BLAS + ATLAS F77LAPACK from C++  demo
   @2009  hdaniel
*/

/*
typedef long int        __CLPK_integer;
typedef long int        __CLPK_logical;
typedef float           __CLPK_real;
typedef double          __CLPK_doublereal;

typedef struct { __CLPK_real r, i; } __CLPK_complex;
typedef struct { __CLPK_doublereal r, i; } __CLPK_doublecomplex;
*/

#include <iostream>
#include <iomanip>
using namespace std;

#include "../matrixutl.h"

extern "C" {
  #include <cblas.h>
  #include <clapack.h>  

  char TRANSA, TRANSB;
  long M, N, K, NRHS, INCX, INCY, LDA, LDB, LDC, INFO;
  double ALPHA, BETA;

  void daxpy_(const long *Np, double *DAp, double *DX, const long *INCXp, 
              double *DY, const long *INCYp);
  void dgemv_(const char *TRANSA, const long *Mp, const long *Np, double *ALPHAp,
              double *A, const long *LDAp, double *X, const long *INCXp, 
              double *BETAp, double *Y, const long *INCYp);
  void dgemm_(const char *TRANSA, const char *TRANSB, const long *Mp, 
              const long *Np, const long *Kp, double *ALPHAp, double *A, 
              const long *LDAp, double *B, const long *LDBp, double *BETAp, 
              double *C, const long *LDCp);
  void dgesv_(const long *Np, const long *NRHSp, double *A, const long *LDAp,
              long *IPIV, double *B, const long *LDBp, long *INFOp );
  void dgtsv_(const long *Np, const long *NRHSp, double *DL, double *D, 
              double *DU, double *B, const long *LDBp, long *INFOp);
}


int main() {
/*type A[] = {
  3, 1, 3,
  1, 5, 9,
  2, 6, 5
};

Fortran alocates matrices on contiguous memory positions by columns and C by rows,
so to use Fortran subroutines matrices must be defined by rows in C, which can be 
done easly defining a transpose matrix in C.
(Since vectors have only one column this is not required)
On C print routine it is required to exchange column and indices to show not
the transpose but the actual matrix stored by columns.
*/
type A[] = {
  3, 1, 2,
  1, 5, 6,
  3, 9, 5
};

type C[] = { 
  0, 0, 0, 
  0, 0, 0, 
  0, 0, 0
};

type x[] = { -1, -1, 1 };

type y[] = { 1, 2, 3 };

type I[] = { 1, 1, 1 };

type I3[] = { 
  1, 0, 0, 
  0, 1, 0,
  0, 0, 1
};

  //operators
  printm ("A", A, 3, 3, 1); //Since stored by columns transpose to show actual matrix in a C program
  printm ("x", x, 3);
  printm ("y", y, 3);  

  /* Level 1 BLAS:vector operations
     man DAXPY, DAXPY, ... (Double A*X+Y)
     xAXPY perform:
     
     y=alpha*x+y
   */

  // y = x + y (ATLAS F77BLAS interface) 
  N=3; ALPHA=1.0; INCX=INCY=1;
  daxpy_(&N, &ALPHA, x, &INCX, y, &INCY);
  printm ("y=x+y", y, 3);

  /* Level 2 BLAS: matrix-vector operations
     man SGEMV, DGEMV, .... (Double GEneral Matrix Vector)
     xGEMV perform:

    y = alpha*A*x + beta*y or 
    y = alpha*A'*x + beta*y
  */      

  // y = A*x    
  TRANSA='N'; M=N=3; ALPHA=1.0; LDA=3; INCX=INCY=1; BETA=0.0;
  dgemv_(&TRANSA, &M, &N, &ALPHA, A, &LDA, x, &INCX, &BETA, y, &INCY);
  printm ("y=A.x", y, 3);

  //y=A(:,1)+...+A(:,3)+y
  TRANSA='N'; M=N=3; ALPHA=1.0; LDA=3; INCX=INCY=1; BETA=1.0;
  dgemv_(&TRANSA, &M, &N, &ALPHA, A, &LDA, I, &INCX, &BETA, y, &INCY);
  printm ("y=A(:,1)+...+A(:,3)+y", y, 3);

  /* Level 3 BLAS: matrix-matrix operations
     man SGEMM, DGEMM, .... (Double GEneral Matrix Matrix)
     xGEMM perform:

    C = alpha*op(A)*op(B) + beta*C, where op(X) can be the tranposed matrix
    specified by 2nd and 3rd parameters
  */      

  // C = A * A
  TRANSA=TRANSB='N'; M=N=K=3; ALPHA=1.0; LDA=LDB=LDC=3; BETA=1.0;
  dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, A, &LDB, &BETA, C, &LDC);
  printm ("C=A*A", C, 3, 3, 1); //Since stored by columns transpose to show actual matrix in a C program
  clear(C, 3, 3);

  // A = A + A
  TRANSA=TRANSB='N'; M=N=K=3; ALPHA=1.0; LDA=LDB=LDC=3; BETA=1.0;
  dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, I3, &LDB, &BETA, A, &LDC);
  printm ("A=A+A", A, 3, 3, 1); //Since stored by columns transpose to show actual matrix in a C program

  //restore previous values of A. Since A was added to itself A/2 will do the trick
  TRANSA=TRANSB='N'; M=N=K=3; ALPHA=0.5; LDA=LDB=LDC=3; BETA=0.0;
  dgemm_(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, I3, &LDB, &BETA, A, &LDC);

  //restore original values of y
  y[0]=1; y[1]=2; y[2]=3; 
    
  printm ("original A", A, 3, 3, 1); //Since stored by columns transpose to show actual matrix in a C program
  printm ("original x", x, 3);        
  printm ("original y", y, 3);

  /* Lapack linear system solver for general matrices:
     man SGESV, DGESV, ... (Double GEneral SolVer)
     
     y = A.x <=> x = inv(A) * y <=> x' = 0.244, 0.477, -0.070
  */
  long ipiv [3]; //permutation matrix  
  N=3; NRHS=1; LDA=LDB=3; 
  dgesv_(&N, &NRHS, A, &LDA, ipiv, y, &LDB, &INFO);
  if (INFO != 0) cerr << "Error = " << INFO << endl;
  else printm ("y=A.x <=> x=inv(A).y", y, 3);

  /* Lapack linear system solver for square tridiagonal matrix (obtained by Gaussian 
     elimination with partial pivoting.
     A tridiagonal matrix is a matrix that is "almost" a diagonal matrix:
        has nonzero elements only in the main diagonal, 
        the first diagonal below this, and 
        the first diagonal above the main diagonal.
       
     man DGTSV, DGTSV, ... (Double General Tridiagonal SolVer)
     (this function has no prototype on atals CLAPACK interface, so
      must be acessed from fortram only)
     
     z = TD.x <=> x = inv(TD) * z <=> x' =  4.8, -3.8, 1.6, 1.2
  */
  /*
  type TD[] = {  
    1, 1, 0, 0,
    1, 2, 3, 0,
    0, 1, 2, 3,
    0, 0, 1, 2
  };
  Store TD by colums for Fortran SUBROUTINE menas transpose matrix in C program:
  */
  type TD[] = {  
    1, 1, 0, 0,
    1, 2, 1, 0,
    0, 3, 2, 1,
    0, 0, 3, 2
  };
  type DU[]= { 1, 3, 3 };
  type D[]= { 1, 2, 2, 2 };  
  type DL[] = { 1, 1, 1 };


  type z[]= { 1, 2, 3, 4 };
  printm ("TD", TD, 4, 4, 1); //Since stored by columns transpose to show actual matrix in a C program
  printm ("z", z, 4);        

  long ipiv2 [4]; //permutation matrix  
  N=4; NRHS=1; LDA=LDB=4;  
  dgesv_(&N, &NRHS, TD, &LDA, ipiv2, z, &LDB, &INFO);
  if (INFO != 0) cerr << "Error = " << INFO << endl;
  else printm ("z=TD.x <=> x=inv(TD).z (using F77 interface DGESV)", z, 4);

  //restore original z vector
  z[0]=1; z[1]=2; z[2]=3; z[3]=4; 

  N = 4; NRHS = 1;
  dgtsv_(&N, &NRHS, DL, D, DU, z, &N, &INFO); 
  if (INFO != 0) cerr << "Error = " << INFO << endl;
  else printm ("z=TD.x <=> x=inv(TD).z (using F77 interface DGTSV)", z, 4);
 
  return 0;
}
