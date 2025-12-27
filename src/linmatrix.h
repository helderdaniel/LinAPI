/* LinAPI matrix class v1.0
   @2009 hdaniel, mmadeira
 */

#ifndef _LINMATRIX_H_
#define _LINMATRIX_H_

#include "operator2d.h"
#include "linvector.h"

template <class T>
class LinMatrix : public Operator2D<T> {
	void abstract() {}; //make this derived class NOT abstract
	void identity();

public:
	LinMatrix () : Operator2D<T> () {}
	LinMatrix (int m, int n, const char* name=NULL, bool ident=false, int mb=2,int nb=2)
					: Operator2D<T> (m,n, name, mb,nb) { if (ident) identity(); };
	LinMatrix (int m, int n, T* data, const char* name=NULL, int mb=2,int nb=2)
					: Operator2D<T> (m,n,data,name, mb,nb) {};
	//Not needed (by default base constructor is called
	//LinMatrix (const LinMatrix<T>& op) : Operator2D<T> (op) {};
	//needed to convert base constructor to this class (there must be a better way)
	//LinMatrix (const Operator2D<T>& op) : Operator2D<T> (op) {};
	LinMatrix<T>& operator= (const LinMatrix<T>& op);

	//Convert vector to more general matrix type
	//Not needed (at least yet)
	//LinVectors when passed as function args are cast Operator2D and then downcast LinMatrix
	//LinMatrix (const LinVector<T>& op);

	//linear algebra
	Operator2D<T>& operator += (const Operator2D<T>& op);
	Operator2D<T>& operator *= (const Operator2D<T>& op);
	//y=alpha*A*x+beta*y
	//y=alpha*A'*x+beta*y
	//x.demv(alpha,beta,y) (optimized Sca/Lapack to accumulate on argument y)
	void gemv (const T& alpha, const T& beta, LinVector<T>& y);
};


/* This form of assignment operator must be implemented on concrete class since
   base class Operator2D<T> tmp(op) can not be instantiated (is abstract)
 */
template <class T>
LinMatrix<T>& LinMatrix<T>::operator= (const LinMatrix<T>& op) {
	LinMatrix<T> tmp(op);
	this->swap(tmp); //call base Operator2D<T> swap function;
	//if exists concrete class data it must be swaped here
	return *this;
}


/* Linear Algebra
 */
template <class T>
const LinMatrix<T> operator + (const LinMatrix<T>& op1, const LinMatrix<T>& op2) {
	LinMatrix<T> tmp (op1);
	tmp += op2;
	return tmp;
}

//This way object is created in function and returned only reference(more efficient?) than
//above operator +, where created object is returned

//Two function defined to allow operation with Matrix or Vector as 2nd operator
//Devia ter dupla distribuição?? para operações com vectores e matrizes?
template <class T>
const LinMatrix<T>& operator * (const LinMatrix<T>& op1, const LinMatrix<T>& op2) {
	LinMatrix<T> *tmp = new LinMatrix<T>(op1);
	(*tmp) *= op2;
	return *tmp;
}
template <class T>
const LinVector<T>& operator * (const LinMatrix<T>& op1, const LinVector<T>& op2) {
	LinMatrix<T> *tmp = new LinMatrix<T>(op1);
	(*tmp) *= op2;
	return *((LinVector<T>*) tmp);
}





//Fix no export keyword for separate templates in *. h and *.cpp
//This allows compiler that implement export keyword or not to
//compile this source code.
//However it can be removed when support arrives at gcc (remove also from *.cpp)
//Makefile define 2 diferent MACROS _SCALAPACK_ and _LAPACK_
//depending on version lapack/scalapack
#ifndef USE_EXPORT_KEYWORD
	#ifdef _SCALAPACK_
   	#include "scalapack/linmatrix_scalapack.cpp"
	#endif
	#ifdef _LAPACK_
		#include "lapack/linmatrix_lapack.cpp"
	#endif
#endif


#endif