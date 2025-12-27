/* LinSolver class for LinAPI
	solves systems of equations
	Lapack version
   @2009 hdaniel mmadeira
 */



#include "../linsolver.h"
#include "../linexception.h"

//Fix no export keyword for separate templates in *. h and *.cpp
//This allows compiler that implement export keyword or not to
//compile this source code.
//However it can be removed when support arrives at gcc (remove also from *.h)
//and remove this filename from exclusion (with patsubst) in makefile 
//so that makefile compiles all *.cpp on this dir
#ifndef _LINSOLVER_LAPACK_
	#define _LINSOLVER_LAPACK_

	#ifndef USE_EXPORT_KEYWORD
   	#define export 
	#endif


/* External LAPACK functions prototypes
*/
extern "C" {
	void sgesv_(const INTEGER *Np, const INTEGER *NRHSp, SINGLE *A, const INTEGER *LDAp,
              INTEGER *IPIV, SINGLE *B, const INTEGER *LDBp, INTEGER *INFOp );
	void dgesv_(const INTEGER *Np, const INTEGER *NRHSp, DOUBLE *A, const INTEGER *LDAp,
              INTEGER *IPIV, DOUBLE *B, const INTEGER *LDBp, INTEGER *INFOp );
	void cgesv_(const INTEGER *Np, const INTEGER *NRHSp, COMPLEX *A, const INTEGER *LDAp,
              INTEGER *IPIV, COMPLEX *B, const INTEGER *LDBp, INTEGER *INFOp );
	void zgesv_(const INTEGER *Np, const INTEGER *NRHSp, DBLCOMPLEX *A, const INTEGER *LDAp,
              INTEGER *IPIV, DBLCOMPLEX *B, const INTEGER *LDBp, INTEGER *INFOp );

	void dgtsv_(const INTEGER *Np, const INTEGER *NRHSp, DOUBLE *DL, DOUBLE *D,
              DOUBLE *DU, DOUBLE *B, const INTEGER *LDBp, INTEGER *INFOp);
	//missing single complex, double trapezoidal functions
}



/* System solver
 */
/* Specializing a wrapper for each Sca/Lapack function that operates on different data types  	allows a template call too the wrapper, avoiding code duplication
	(need to have default template since it is not declared anywhere (e.g. as member function)
 */
template <class T> 
void gesvw (INTEGER N, INTEGER NRHS, T *A, INTEGER LDA, INTEGER *ipiv,
				T *Y, INTEGER LDB, INTEGER *INFO) {
	//this should not happen since Operator2D contructor prevents
	//creation of unsupported types
	throw LinTypeNotSupportedException<T, LinSolver> ("gesv");
}
template <>
void gesvw<SINGLE> (INTEGER N, INTEGER NRHS, SINGLE *A, INTEGER LDA, INTEGER *ipiv,
							SINGLE *Y, INTEGER LDB, INTEGER *INFO) {
	sgesv_(&N, &NRHS, A, &LDA, ipiv, Y, &LDB, INFO);
}
template <>
void gesvw<DOUBLE> (INTEGER N, INTEGER NRHS, DOUBLE *A, INTEGER LDA, INTEGER *ipiv,
							DOUBLE *Y, INTEGER LDB, INTEGER *INFO) {
	dgesv_(&N, &NRHS, A, &LDA, ipiv, Y, &LDB, INFO);
}
template <>
void gesvw<COMPLEX> (INTEGER N, INTEGER NRHS, COMPLEX *A, INTEGER LDA, INTEGER *ipiv,
							COMPLEX *Y, INTEGER LDB, INTEGER *INFO) {
	cgesv_(&N, &NRHS, A, &LDA, ipiv, Y, &LDB, INFO);
}
template <>
void gesvw<DBLCOMPLEX> (INTEGER N, INTEGER NRHS, DBLCOMPLEX *A, INTEGER LDA, INTEGER *ipiv,
								DBLCOMPLEX *Y, INTEGER LDB, INTEGER *INFO) {
	zgesv_(&N, &NRHS, A, &LDA, ipiv, Y, &LDB, INFO);
}

export template <class T> 
const LinVector<T>& LinSolver::gesv (const LinMatrix<T>& A, const LinVector<T>& op) {
const INTEGER N = A.rows();
const INTEGER NRHS = op.cols();
const INTEGER LDA=max((INTEGER) 1, N);
const INTEGER LDB=max((INTEGER) 1, N);
INTEGER *ipiv = new INTEGER [N]; //permutation matrix
INTEGER INFO;
LinVector<T> *y = new LinVector<T>(op); //do not destruct op

	gesvw(N, NRHS, A._data, LDA, ipiv, (*y)._data, LDB, &INFO);
	delete [] ipiv;
	if (INFO != 0) throw LinException("Lapack solver ?gesv error ", INFO);

	return *y;
}


#endif