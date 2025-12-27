/* LinAPI vector class v1.0
   @2009 hdaniel, mmadeira
 */

#ifndef _LINVECTOR_H_
#define _LINVECTOR_H_

#include "operator2d.h"

template <class T>
class LinVector : public Operator2D<T> {
	void abstract() {}; //make this derived class NOT abstract
	void identity() { setAll((T) 1); };

public:
	LinVector () : Operator2D<T> () {}
	LinVector (int m, const char* name=NULL, bool ident=false, int mb=2,int nb=1)
					: Operator2D<T> (m, 1, name, mb,nb) { if (ident) identity(); };
	LinVector (int m, T* data, const char* name=NULL, int mb=2,int nb=1)
					: Operator2D<T> (m, 1, data, name, mb,nb) {};
	//Not needed (by default base constructor is called
	//LinVector (const LinVector<T>& op) : Operator2D<T> (op) {};
	//needed to convert base constructor to this class 
	//LinVector (const Operator2D<T>& op) : Operator2D<T> (op) {};
	LinVector<T>& operator= (const LinVector<T>& op);

	//linear algebra
	Operator2D<T>& operator += (const Operator2D<T>& op);
	//y=alpha*x+y
	//x.axpy(alpha,y) (optimized Sca/Lapack to accumulate on argument y) 
	void axpy (const T& alpha, LinVector<T>& y);
};


/* This form of assignment operator must be implemented on concrete class since
   base class Operator2D<T> tmp(op) can not be instantiated (is abstract)
 */
template <class T>
LinVector<T>& LinVector<T>::operator= (const LinVector<T>& op)  {
	LinVector<T> tmp(op);
	this->swap(tmp); //call base Operator2D<T> tmpswap function;
	//if exists concrete class data it must be swaped here
	return *this;
}


/* Linear Algebra
 */
template <class T>
const LinVector<T> operator + (const LinVector<T>& op1, const LinVector<T>& op2) {
	LinVector<T> tmp (op1);
	tmp += op2;
	return tmp;
}


//Fix no export keyword for separate templates in *. h and *.cpp
//This allows compiler that implement export keyword or not to
//compile this source code.
//However it can be removed when support arrives at gcc (remove also from *.cpp)
//Makefile define 2 diferent MACROS _SCALAPACK_ and _LAPACK_
//depending on version lapack/scalapack
#ifndef USE_EXPORT_KEYWORD
	#ifdef _SCALAPACK_
   	#include "scalapack/linvector_scalapack.cpp"
	#endif
	#ifdef _LAPACK_
		#include "lapack/linvector_lapack.cpp"
	#endif
#endif



#endif