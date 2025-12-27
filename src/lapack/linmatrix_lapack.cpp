/* LinAPI matrix class v1.0
   Lapack version
   @2009 hdaniel mmadeira
 */

#include "../linmatrix.h"
#include "../linexception.h"

//Fix no export keyword for separate templates in *. h and *.cpp
//This allows compiler that implement export keyword or not to
//compile this source code.
//However it can be removed when support arrives at gcc (remove also from *.h)
//and remove this filename from exclusion (with patsubst) in makefile 
//so that makefile compiles all *.cpp on this dir 
#ifndef _LINMATRIX_LAPACK_
	#define _LINMATRIX_LAPACK_

	#ifndef USE_EXPORT_KEYWORD
   #define export 
#endif


/* External LAPACK functions prototypes
*/
extern "C" {
	//Level 2 BLAS
	//y = alpha*op(A)*x + beta*y, where op(A) can be the tranposed matrix
	void sgemv_(const CHAR *TRANSA, const INTEGER *Mp, const INTEGER *Np,
	            SINGLE *ALPHAp, SINGLE *A, const INTEGER *LDAp, SINGLE *X,
	            const INTEGER *INCXp, SINGLE *BETAp, SINGLE *Y, const INTEGER *INCYp);
	void dgemv_(const CHAR *TRANSA, const INTEGER *Mp, const INTEGER *Np,
	            DOUBLE *ALPHAp, DOUBLE *A, const INTEGER *LDAp, DOUBLE *X,
	            const INTEGER *INCXp, DOUBLE *BETAp, DOUBLE *Y, const INTEGER *INCYp);
	void cgemv_(const CHAR *TRANSA, const INTEGER *Mp, const INTEGER *Np,
	            COMPLEX *ALPHAp, COMPLEX *A, const INTEGER *LDAp, COMPLEX *X,
	            const INTEGER *INCXp, COMPLEX *BETAp, COMPLEX *Y, const INTEGER *INCYp);
	void zgemv_(const CHAR *TRANSA, const INTEGER *Mp, const INTEGER *Np,
	            DBLCOMPLEX *ALPHAp, DBLCOMPLEX *A, const INTEGER *LDAp, DBLCOMPLEX *X,
	            const INTEGER *INCXp, DBLCOMPLEX *BETAp, DBLCOMPLEX *Y, const INTEGER *INCYp);

	//Level 3 BLAS
	//C = alpha*op(A)*op(B) + beta*C, where op(X) can be the tranposed matrix
	void sgemm_(const CHAR *TRANSA, const CHAR *TRANSB, const INTEGER *Mp,
	            const INTEGER *Np, const INTEGER *Kp, const SINGLE *ALPHAp, SINGLE *A,
	            const INTEGER *LDAp, SINGLE *B, const INTEGER *LDBp, const SINGLE *BETAp,
	            SINGLE *C, const INTEGER *LDCp);
	void dgemm_(const CHAR *TRANSA, const CHAR *TRANSB, const INTEGER *Mp,
	            const INTEGER *Np, const INTEGER *Kp, const DOUBLE *ALPHAp, DOUBLE *A,
	            const INTEGER *LDAp, DOUBLE *B, const INTEGER *LDBp, const DOUBLE *BETAp,
	            DOUBLE *C, const INTEGER *LDCp);
	void cgemm_(const CHAR *TRANSA, const CHAR *TRANSB, const INTEGER *Mp,
	            const INTEGER *Np, const INTEGER *Kp, const COMPLEX *ALPHAp, COMPLEX *A,
	            const INTEGER *LDAp, COMPLEX *B, const INTEGER *LDBp, const COMPLEX *BETAp,
	            COMPLEX *C, const INTEGER *LDCp);
	void zgemm_(const CHAR *TRANSA, const CHAR *TRANSB, const INTEGER *Mp,
	            const INTEGER *Np, const INTEGER *Kp, const DBLCOMPLEX *ALPHAp, DBLCOMPLEX *A,
	            const INTEGER *LDAp, DBLCOMPLEX *B, const INTEGER *LDBp, const
	            DBLCOMPLEX *BETAp, DBLCOMPLEX *C, const INTEGER *LDCp);
}

/* Wrapper functions
 */
//GEMV
template <class T> //need to have default template
void gemvw (CHAR TRANSA, INTEGER Mp, INTEGER Np, T ALPHAp, T *A,
				INTEGER LDAp, T *X, INTEGER INCXp, T BETAp, T *Y,
				INTEGER INCYp) {
	//this should not happen since Operator2D contructor prevents
	//creation of unsupported types
	throw LinTypeNotSupportedException<T, LinMatrix<T> > ("gemv");
}
template <>
void gemvw<SINGLE> (CHAR TRANSA, INTEGER Mp, INTEGER Np, SINGLE ALPHAp, SINGLE *A,
							INTEGER LDAp, SINGLE *X, INTEGER INCXp, SINGLE BETAp, SINGLE *Y,
							INTEGER INCYp) {
	sgemv_(&TRANSA, &Mp, &Np, &ALPHAp, A, &LDAp, X, &INCXp, &BETAp, Y, &INCYp);
}
template <>
void gemvw<DOUBLE> (CHAR TRANSA, INTEGER Mp, INTEGER Np, DOUBLE ALPHAp, DOUBLE *A,
                    INTEGER LDAp, DOUBLE *X, INTEGER INCXp, DOUBLE BETAp, DOUBLE *Y,
                    INTEGER INCYp) {
	dgemv_(&TRANSA, &Mp, &Np, &ALPHAp, A, &LDAp, X, &INCXp, &BETAp, Y, &INCYp);
}
template <>
void gemvw<COMPLEX> (CHAR TRANSA, INTEGER Mp, INTEGER Np, COMPLEX ALPHAp, COMPLEX *A,
                     INTEGER LDAp, COMPLEX *X, INTEGER INCXp, COMPLEX BETAp, COMPLEX *Y,
                     INTEGER INCYp) {
	cgemv_(&TRANSA, &Mp, &Np, &ALPHAp, A, &LDAp, X, &INCXp, &BETAp, Y, &INCYp);
}
template <>
void gemvw<DBLCOMPLEX> (CHAR TRANSA, INTEGER Mp, INTEGER Np, DBLCOMPLEX ALPHAp,
                        DBLCOMPLEX *A, INTEGER LDAp, DBLCOMPLEX *X, INTEGER INCXp,
                        DBLCOMPLEX BETAp, DBLCOMPLEX *Y, INTEGER INCYp) {
	zgemv_(&TRANSA, &Mp, &Np, &ALPHAp, A, &LDAp, X, &INCXp, &BETAp, Y, &INCYp);
}


//GEMM
template <class T> //need to have default template
void gemmw (CHAR TRANSA, CHAR TRANSB, INTEGER Mp, INTEGER Np, INTEGER Kp,
				T ALPHAp, T *A, INTEGER LDAp, T *B, INTEGER LDBp, T BETAp,
				T *C, INTEGER LDCp) {
	throw LinTypeNotSupportedException<T, LinMatrix<T> > ("gemm");
}
template <>
void gemmw<SINGLE> (CHAR TRANSA, CHAR TRANSB, INTEGER Mp, INTEGER Np, INTEGER Kp,
                    SINGLE ALPHAp, SINGLE *A, INTEGER LDAp, SINGLE *B, INTEGER LDBp,
                    SINGLE BETAp, SINGLE *C, INTEGER LDCp) {
	sgemm_ (&TRANSA, &TRANSB, &Mp, &Np, &Kp, &ALPHAp, A, &LDAp, B, &LDBp, &BETAp, C, &LDCp);
}
template <>
void gemmw<DOUBLE> (CHAR TRANSA, CHAR TRANSB, INTEGER Mp, INTEGER Np, INTEGER Kp,
							DOUBLE ALPHAp, DOUBLE *A, INTEGER LDAp, DOUBLE *B, INTEGER LDBp,
							DOUBLE BETAp, DOUBLE *C, INTEGER LDCp) {
	dgemm_ (&TRANSA, &TRANSB, &Mp, &Np, &Kp, &ALPHAp, A, &LDAp, B, &LDBp, &BETAp, C, &LDCp);
}
template <>
void gemmw<COMPLEX> (CHAR TRANSA, CHAR TRANSB, INTEGER Mp, INTEGER Np, INTEGER Kp,
                     COMPLEX ALPHAp, COMPLEX *A, INTEGER LDAp, COMPLEX *B, INTEGER LDBp,
                     COMPLEX BETAp, COMPLEX *C, INTEGER LDCp) {
	cgemm_ (&TRANSA, &TRANSB, &Mp, &Np, &Kp, &ALPHAp, A, &LDAp, B, &LDBp, &BETAp, C, &LDCp);
}
template <>
void gemmw<DBLCOMPLEX> (CHAR TRANSA, CHAR TRANSB, INTEGER Mp, INTEGER Np, INTEGER Kp,
                        DBLCOMPLEX ALPHAp, DBLCOMPLEX *A, INTEGER LDAp, DBLCOMPLEX *B,
                        INTEGER LDBp, DBLCOMPLEX BETAp, DBLCOMPLEX *C, INTEGER LDCp) {
	zgemm_ (&TRANSA, &TRANSB, &Mp, &Np, &Kp, &ALPHAp, A, &LDAp, B, &LDBp, &BETAp, C, &LDCp);
}

// End wrapper functions




export template <class T>
void LinMatrix<T>::identity() {
//const int elements = _localRows*_localCols;
	for (int i=0; i<this->_rows; ++i)
		for (int j=0; j<this->_cols; ++j)	
			if (i==j)	//Set only diagonal (all others are already zero
				this->_data [j*this->_rows+i] = (T) 1; //Transpose (C / F77 matrix alloc)
}


/* Linear algebra
 */
/* Matrix, Matrix add
 */
export template <class T>
Operator2D<T>& LinMatrix<T>::operator += (const Operator2D<T>& op) {
const LinMatrix<T>& y = static_cast<const LinMatrix<T>&> (op);
	for (int i=0; i<this->_localElements; ++i) (*this)._data[i] += y._data[i];
	return *this;
}

//implementar função para usar dgemv, como nos vectores para
export template <>
void LinMatrix<SINGLE>::gemv (const SINGLE& alpha, const SINGLE& beta, LinVector<SINGLE>& y) {
	
};


/* Matrix, Vector product
	y = A*x
 */
//usando dgemv numa funcção separada garante que
//o retorno é um vector e num uma matrix


/* Matrix, MAtrix product
 */
//C = A*B 
//Operator *= não pode ser sobrecarregado, pq a multiplicação de A por B 
//poderá resultar numa matriz com dimensão diferente do receptor: A
//Assim o operador * deverá ser não membro

/* Pre this->cols()==op.rows();
 */
export template <class T>
Operator2D<T>& LinMatrix<T>::operator *= (const Operator2D<T>& op) {
const LinMatrix<T>& B = static_cast<const LinMatrix<T>&> (op);

const CHAR trans='N';
const INTEGER M = (*this).rows();
const INTEGER K = (*this).cols(); // or B.rows();
const INTEGER N = B.cols();
const INTEGER LDA=max((INTEGER) 1, M); //or max(1,K) if this transposed
const INTEGER LDB=max((INTEGER) 1, K); //or max(1,N) if op transposed
const INTEGER LDC=max((INTEGER) 1, M); //always
const T alpha=1;
const T beta=1;
LinMatrix C(M, N); //Create result

	gemmw(trans, trans, M, N, K, alpha, (*this)._data, LDA, B._data, LDB, beta, C._data, LDC);
		
	this->swap(C); //swap result with this
	return *this;
}

#endif