/* LinAPI vector class v1.0
	Scalapack version
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
#ifndef _LINVECTOR_SCALAPACK_
	#define _LINVECTOR_SCALAPACK_

	#ifndef USE_EXPORT_KEYWORD
   	#define export 
	#endif


/* External SCALAPACK functions prototypes
*/
extern "C" {
	//Level 1 PBLAS
	//y=alpha*x+y
	void psaxpy_(const INTEGER *Np, const SINGLE *DAp, SINGLE *DX, const INTEGER *IXp,
					 const INTEGER *JXp, INTEGER *DESCX, const INTEGER *INCXp, SINGLE *DY,
					 const INTEGER *IYp, const INTEGER *JYp, INTEGER *DESCY, const INTEGER *INCYp);
	void pdaxpy_(const INTEGER *Np, const DOUBLE *DAp, DOUBLE *DX, const INTEGER *IXp,
					 const INTEGER *JXp, INTEGER *DESCX, const INTEGER *INCXp, DOUBLE *DY,
					 const INTEGER *IYp, const INTEGER *JYp, INTEGER *DESCY, const INTEGER *INCYp);
	void pcaxpy_(const INTEGER *Np, const COMPLEX *DAp, COMPLEX *DX, const INTEGER *IXp,
					 const INTEGER *JXp, INTEGER *DESCX, const INTEGER *INCXp, COMPLEX *DY,
					 const INTEGER *IYp, const INTEGER *JYp, INTEGER *DESCY, const INTEGER *INCYp);
	void pzaxpy_(const INTEGER *Np, const DBLCOMPLEX *DAp, DBLCOMPLEX *DX, const INTEGER *IXp,
					 const INTEGER *JXp, INTEGER *DESCX, const INTEGER *INCXp, DBLCOMPLEX *DY,
					 const INTEGER *IYp, const INTEGER *JYp, INTEGER *DESCY, const INTEGER *INCYp);
}


/* Linear algebra
 */

/* p?axpy_ wrapper
   (need to have default template since it is not declared anywhere (e.g. as member function)
 */
template <class T> 
void axpyw (const INTEGER Np, const T DAp, T *DX, const INTEGER IXp,
				const INTEGER JXp, INTEGER *DESCX, const INTEGER INCXp, T *DY,
				const INTEGER IYp, const INTEGER JYp, INTEGER *DESCY, const INTEGER INCYp) {
	//this should not happen since Operator2D contructor prevents
	//creation of unsupported types
	throw LinTypeNotSupportedException<T, LinVector<T> > ("axpy");
}
template <>
void axpyw<SINGLE> (const INTEGER Np, const SINGLE DAp, SINGLE *DX, const INTEGER IXp,
				const INTEGER JXp, INTEGER *DESCX, const INTEGER INCXp, SINGLE *DY,
				const INTEGER IYp, const INTEGER JYp, INTEGER *DESCY, const INTEGER INCYp) {
	psaxpy_(&Np, &DAp, DX, &IXp, &JXp, DESCX, &INCXp, DY, &IYp, &JYp, DESCY, &INCYp);
}
template <>
void axpyw<DOUBLE> (const INTEGER Np, const DOUBLE DAp, DOUBLE *DX, const INTEGER IXp,
				const INTEGER JXp, INTEGER *DESCX, const INTEGER INCXp, DOUBLE *DY,
				const INTEGER IYp, const INTEGER JYp, INTEGER *DESCY, const INTEGER INCYp) {
	pdaxpy_(&Np, &DAp, DX, &IXp, &JXp, DESCX, &INCXp, DY, &IYp, &JYp, DESCY, &INCYp);
}
template <>
void axpyw<COMPLEX> (const INTEGER Np, const COMPLEX DAp, COMPLEX *DX, const INTEGER IXp,
				const INTEGER JXp, INTEGER *DESCX, const INTEGER INCXp, COMPLEX *DY,
				const INTEGER IYp, const INTEGER JYp, INTEGER *DESCY, const INTEGER INCYp) {
	pcaxpy_(&Np, &DAp, DX, &IXp, &JXp, DESCX, &INCXp, DY, &IYp, &JYp, DESCY, &INCYp);
}
template <>
void axpyw<DBLCOMPLEX> (const INTEGER Np, const DBLCOMPLEX DAp, DBLCOMPLEX *DX,
				const INTEGER IXp, const INTEGER JXp, INTEGER *DESCX, const INTEGER INCXp,
				DBLCOMPLEX *DY, const INTEGER IYp, const INTEGER JYp, INTEGER *DESCY,
				const INTEGER INCYp) {
	pzaxpy_(&Np, &DAp, DX, &IXp, &JXp, DESCX, &INCXp, DY, &IYp, &JYp, DESCY, &INCYp);
}

/* Vector add

   Simple add element by element (as Lapack version):

	for (int i=0; i<this->_localElements; ++i) (*this)._data[i] += y._data[i];

	do not work if _rowsBlockFactor or _colsBlockFactor differ in receptor and argument
 	Must use specialized scalapack function, such as p?daxpy for vector addition
 */
export template <class T>
Operator2D<T>& LinVector<T>::operator += (const Operator2D<T>& op) {
const LinVector<T>& y = static_cast<const LinVector<T>&> (op);
T alpha = 1;

	axpyw((INTEGER) this->_rows, alpha, y._data, (INTEGER) this->subRow,
			(INTEGER) this->subCol, y._gridDesc, (INTEGER) this->interleave, (*this)._data,
			(INTEGER) this->subRow, (INTEGER) this->subCol, (*this)._gridDesc,
			(INTEGER) this->interleave);
	return *this;
}

export template <class T>
void LinVector<T>::axpy (const T& alpha, LinVector<T>& y) {
	axpyw((INTEGER) this->_rows, alpha, (*this)._data, (INTEGER ) this->subRow,
			(INTEGER) this->subCol, (*this)._gridDesc, (INTEGER) this->interleave, y._data, (INTEGER)  this->subRow, (INTEGER) this->subCol, y._gridDesc,
			(INTEGER) this->interleave);
}


/* Dot product
 */

/* Norm
 */

// ...

#endif