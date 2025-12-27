/* LinSolver class for LinAPI
	solves systems of equations
   @2009 hdaniel mmadeira
 */

#ifndef _LINSOLVER_H_
#define _LINSOLVER_H_

#include "linmatrix.h"
#include "linvector.h"

//LinSolver static class
//template <class T>
class LinSolver {
public:
	template <class T>
	static const LinVector<T>& gesv (const LinMatrix<T>& A, const LinVector<T>& y);
	//template <class T>
	//static const LinMatrix<T>& gesv (const LinMatrix<T>& A, const LinMatrix<T>& Y);
};


//Fix no export keyword for separate templates in *. h and *.cpp
//This allows compiler that implement export keyword or not to
//compile this source code.
//However it can be removed when support arrives at gcc (remove also from *.cpp)
//Makefile define 2 diferent MACROS _SCALAPACK_ and _LAPACK_
//depending on version lapack/scalapack
#ifndef USE_EXPORT_KEYWORD
	#ifdef _SCALAPACK_
    #include "scalapack/linsolver_scalapack.cpp"
	#endif
	#ifdef _LAPACK_
		#include "lapack/linsolver_lapack.cpp"
	#endif
#endif


#endif
