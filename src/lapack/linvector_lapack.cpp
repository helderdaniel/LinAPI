/* LinAPI vector class v1.0
   Lapack version
   @2009 hdaniel mmadeira
 */


#include "../linvector.h"
#include "../linexception.h"


//Fix no export keyword for separate templates in *. h and *.cpp
//This allows compiler that implement export keyword or not to
//compile this source code.
//However it can be removed when support arrives at gcc (remove also from *.h)
//and remove this filename from exclusion (with patsubst) in makefile 
//so that makefile compiles all *.cpp on this dir
#ifndef _LINVECTOR_LAPACK_
	#define _LINVECTOR_LAPACK_

	#ifndef USE_EXPORT_KEYWORD
   	#define export 
	#endif


/* External LAPACK functions prototypes
*/
extern "C" {
	//Level 1 BLAS
	//y=alpha*x+y
	void saxpy_(const INTEGER *Np, const SINGLE *DAp, SINGLE *DX,
					const INTEGER *INCXp, SINGLE *DY, const INTEGER *INCYp);
	void daxpy_(const INTEGER *Np, const DOUBLE *DAp, DOUBLE *DX,
					const INTEGER *INCXp, DOUBLE *DY, const INTEGER *INCYp);
	void caxpy_(const INTEGER *Np, const COMPLEX *DAp, COMPLEX *DX,
					const INTEGER *INCXp, COMPLEX *DY, const INTEGER *INCYp);
	void zaxpy_(const INTEGER *Np, const DBLCOMPLEX *DAp, DBLCOMPLEX *DX,
					const INTEGER *INCXp, DBLCOMPLEX *DY, const INTEGER *INCYp);
}

/* Linear algebra
 */

/* Vector add
 */
export template <class T>
Operator2D<T>& LinVector<T>::operator += (const Operator2D<T>& op) {
const LinVector<T>& y = static_cast<const LinVector<T>&> (op);

	for (int i=0; i<this->_localElements; ++i) (*this)._data[i] += y._data[i];
	/* OR using LAPACK functions (must be splited in 4 function):
		T alpha = 1;
		?axpy_((INTEGER *) &this->_rows, &alpha, y._data, (INTEGER *) &y->interleave, (*this)._data, (INTEGER *) &this->interleave);
	 */
	return *this;
}


export template <>
void LinVector<SINGLE>::axpy (const SINGLE& alpha, LinVector<SINGLE>& y) {
	saxpy_((INTEGER *) &this->_rows, &alpha, (*this)._data, (INTEGER *) &this->interleave, y._data, (INTEGER *) &this->interleave);
}

export template <>
void LinVector<DOUBLE>::axpy (const DOUBLE& alpha, LinVector<DOUBLE>& y) {
	daxpy_((INTEGER *) &this->_rows, &alpha, (*this)._data, (INTEGER *) &this->interleave, y._data, (INTEGER *) &this->interleave);
}

export template <>
void LinVector<COMPLEX>::axpy (const COMPLEX& alpha, LinVector<COMPLEX>& y) {
	caxpy_((INTEGER *) &this->_rows, &alpha, (*this)._data, (INTEGER *) &this->interleave, y._data, (INTEGER *) &this->interleave);
}

export template <>
void LinVector<DBLCOMPLEX>::axpy (const DBLCOMPLEX& alpha, LinVector<DBLCOMPLEX>& y) {
	zaxpy_((INTEGER *) &this->_rows, &alpha, (*this)._data, (INTEGER *) &this->interleave, y._data, (INTEGER *) &this->interleave);
}

/* Dot product
 */

/* Norm
 */

// ...

#endif