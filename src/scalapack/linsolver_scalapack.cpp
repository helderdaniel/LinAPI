/* LinSolver class for LinAPI
	solves systems of equations
	ScaLapack version
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
#ifndef _LINSOLVER_SCALAPACK_
	#define _LINSOLVER_SCALAPACK_

	#ifndef USE_EXPORT_KEYWORD
   	#define export 
	#endif


/* External ScaLAPACK functions prototypes
*/
extern "C" {
	//SCALAPACK driver functions
	void psgesv_(const INTEGER *Np, const INTEGER *NRHSp, SINGLE *A, const INTEGER *IAp,
					const INTEGER *JAp, INTEGER *DESCA, INTEGER *IPIV, SINGLE *Y, const INTEGER *IYp, const INTEGER *JYp, INTEGER *DESCY, INTEGER *INFOp);
	void pdgesv_(const INTEGER *Np, const INTEGER *NRHSp, DOUBLE *A, const INTEGER *IAp,
					const INTEGER *JAp, INTEGER *DESCA, INTEGER *IPIV, DOUBLE *Y, const INTEGER *IYp, const INTEGER *JYp, INTEGER *DESCY, INTEGER *INFOp);
	void pcgesv_(const INTEGER *Np, const INTEGER *NRHSp, COMPLEX *A, const INTEGER *IAp,
					const INTEGER *JAp, INTEGER *DESCA, INTEGER *IPIV, COMPLEX *Y, const INTEGER *IYp, const INTEGER *JYp, INTEGER *DESCY, INTEGER *INFOp);
	void pzgesv_(const INTEGER *Np, const INTEGER *NRHSp, DBLCOMPLEX *A, const INTEGER *IAp,
					const INTEGER *JAp, INTEGER *DESCA, INTEGER *IPIV, DBLCOMPLEX *Y,
					const INTEGER *IYp, const INTEGER *Jyp, INTEGER *DESCY, INTEGER *INFOp);

	//missing all trapezoidal solver functions ..
	
}


/* System solver
 */
/* Specializing a wrapper for each Sca/Lapack function that operates on different data types  	allows a template call too the wrapper, avoiding code duplication
	(need to have default template since it is not declared anywhere (e.g. as member function)
 */
template <class T> 
void gesvw (INTEGER N, INTEGER NRHS, T *A, INTEGER IAp, INTEGER JAp, INTEGER *DESCA,
				INTEGER *ipiv, T *Y, INTEGER IYp, INTEGER JYp, INTEGER *DESCY,
				INTEGER *INFO) {
	//this should not happen since Operator2D contructor prevents
	//creation of unsupported types
	throw LinTypeNotSupportedException<T, LinSolver> ("gesv");
}
template <>
void gesvw<SINGLE> (INTEGER N, INTEGER NRHS, SINGLE *A, INTEGER IAp, INTEGER JAp,
						   INTEGER *DESCA, INTEGER *ipiv, SINGLE *Y, INTEGER IYp,
							INTEGER JYp, INTEGER *DESCY, INTEGER *INFO) {
	psgesv_(&N, &NRHS, A, &IAp, &JAp, DESCA, ipiv, Y, &IYp, &JYp, DESCY, INFO);
}
template <>
void gesvw<DOUBLE> (INTEGER N, INTEGER NRHS, DOUBLE *A, INTEGER IAp, INTEGER JAp,
						   INTEGER *DESCA, INTEGER *ipiv, DOUBLE *Y, INTEGER IYp,
							INTEGER JYp, INTEGER *DESCY, INTEGER *INFO) {
	pdgesv_(&N, &NRHS, A, &IAp, &JAp, DESCA, ipiv, Y, &IYp, &JYp, DESCY, INFO);
}
template <>
void gesvw<COMPLEX> (INTEGER N, INTEGER NRHS, COMPLEX *A, INTEGER IAp, INTEGER JAp,
						   INTEGER *DESCA, INTEGER *ipiv, COMPLEX *Y, INTEGER IYp,
							INTEGER JYp, INTEGER *DESCY, INTEGER *INFO) {
	pcgesv_(&N, &NRHS, A, &IAp, &JAp, DESCA, ipiv, Y, &IYp, &JYp, DESCY, INFO);
}
template <>
void gesvw<DBLCOMPLEX> (INTEGER N, INTEGER NRHS, DBLCOMPLEX *A, INTEGER IAp, INTEGER JAp,
						   INTEGER *DESCA, INTEGER *ipiv, DBLCOMPLEX *Y, INTEGER IYp,
							INTEGER JYp, INTEGER *DESCY, INTEGER *INFO) {
	pzgesv_(&N, &NRHS, A, &IAp, &JAp, DESCA, ipiv, Y, &IYp, &JYp, DESCY, INFO);
}


export template <class T> 
const LinVector<T>& LinSolver::gesv (const LinMatrix<T>& A, const LinVector<T>& op) {
const INTEGER N = A.rows();
const INTEGER NRHS = op.cols();
const INTEGER LDA=max((INTEGER) 1, N);
const INTEGER LDB=max((INTEGER) 1, N);
const INTEGER subRow=Operator2D<T>::subRow;
const INTEGER subCol=Operator2D<T>::subCol;
INTEGER *ipiv = new INTEGER [N]; //permutation matrix
INTEGER INFO;
LinVector<T> *y = new LinVector<T>(op); //do not destruct op

	gesvw(N, NRHS, A._data, subRow, subCol, A._gridDesc, ipiv, (*y)._data,
			subRow, subCol, (*y)._gridDesc, &INFO);
	delete [] ipiv;
	if (INFO != 0) throw LinException("Lapack solver ?gesv error ", INFO);

	return *y;
}


#endif