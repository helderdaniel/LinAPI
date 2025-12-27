/* LinAPI 2D operators base class v1.0
   Lapack version
   @2009 hdaniel mmadeira
 */

#include "../operator2d.h" 
/*#include "mpiobj.h"
#include <cmath>
#include <cctype>
*/
#include <sstream>
#include <iomanip>
#include <cstdlib>

//Fix no export keyword for separate templates in *. h and *.cpp
//This allows compiler that implement export keyword or not to
//compile this source code.
//However it can be removed when support arrives at gcc (remove also from *.h)
//and remove this filename from exclusion (with patsubst) in makefile 
//so that makefile compiles all *.cpp on this dir 
#ifndef _OPERATOR2D_LAPACK_
	#define _OPERATOR2D_LAPACK_

	#ifndef USE_EXPORT_KEYWORD
   	#define export 
	#endif

/* External LAPACK functions prototypes
*/
extern "C" {

}

export template <class T>
Operator2Ddesc* Operator2D<T>::initDesc (int MB, int NB, int ICTXT, int LLD) {
	return NULL; //scalapack matrix descriptos not needed in Lapack version
}

export template <class T>
void Operator2D<T>::init(int m, int n, const char* name, int mb, int nb) {

	//check valid Sca/Lapack type
	if (!isValidLinType<T>())
		throw LinTypeNotSupportedException<T, typeof(*this)> ("constructor");

	//print format
	_width=defaultWidth;
	_precision=defaultPrecision;

	//Global operator info
	_rows=m; _cols=n;
	//ignore arguments since operator is not block cyclic distributed
	_rowsBlockFactor=1; _colsBlockFactor=1; 
	if (name==NULL) {
		stringstream ss;
		ss << opBaseName << setw(opDigits) << opCounter++;
		_name = ss.str();
	}
	else _name = string(name);

	//Local operator and process grid info
	_localRows=_rows;
	_localCols=_cols;

	//Scalapack matrix descriptors (not used since operator is not distributed)
	_gridDesc = NULL;
	_rootDesc = NULL;

	//Init operator data area
	_localElements=_rows*_cols;
	_data = new T[_localElements];
}


export template <class T> Operator2D<T>::
Operator2D	(int m, int n, T* data, const char* name, int mb, int nb) : blacs(Blacs::instance()) {

	init (m, n, name, mb, nb);
	//Transposed is required since C allocates matrices by rows and Fortran by columns
	linearTranspose(data, _data);
}


/* Distribute operators on process grid and collect them on root process
   Not need to use it in Lapack since operator is not block cyclic distributed
	Defined with no action for version compatibility only
 */

export template <class T> void Operator2D<T>::distribute(T* globalData) const { }
export template <class T> void Operator2D<T>::collect(T* globalData) const { }

/* Print operator 
 */
export template <class T>
ostream& Operator2D<T>::toStream(ostream &os, const bool transpose) const {
	return linearToStream(os, _data, transpose);
}

/* Access individual elements
	(transpose on addressing elements, since F77 stores by columns and C by rows)
 */
export template <class T>
void Operator2D<T>::set(int r, int c, T alpha) { _data [c*_rows+r]=alpha; }
	
export template <class T>
T Operator2D<T>::get(int r, int c) const { return _data [c*_rows+r]; }


#endif