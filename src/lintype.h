/* LinAPI operator types for C++ interface v1.0
   @2009 hdaniel, mmadeira
 */

#ifndef _LINTYPE_H_
#define _LINTYPE_H_

#include <complex>

/* Fortran DATA TYPES equivalence to C++

   In wrappers for external Fortran functions, arguments must be declared as 
   pointers of this types:
*/
typedef long INTEGER;
typedef char CHAR;
typedef float REAL;   //single precision float (same as SINGLE)
typedef float SINGLE; 
typedef double DOUBLE;
//the following complex types report same storage has Fortran 8 and 16 bytes with sizeof()
//Also it seams that members have same order {real,imag}
//(see /usr/include/c++/4.3/complex) 
typedef complex<float> COMPLEX;
typedef complex<double> DBLCOMPLEX; //COMPLEX*16
//The following complex types have the same storage space as Fortan
//CAUTION: operator << is not overloaded as above template complex<>
//struct COMPLEX { REAL r, i; };
//struct DBLCOMPLEX { DOUBLE r, i; };  //COMPLEX*16
typedef long LOGICAL;

/* Force typeof to be same as __typeof__
   typeof is available only with GNU extension to C99 enabled (e.g.  -fasm)
   typeof is NOT available in strict C99 (-std=c99)
	__typeof__ is available without GNU extensions
 */
#define typeof __typeof__


//return true if is valid Sca/Lapack arithmetic type
#include <typeinfo>

template <class T>
bool isValidLinType() {
	if (typeid(T)==typeid(SINGLE)) return true;
	if (typeid(T)==typeid(DOUBLE)) return true;
	if (typeid(T)==typeid(COMPLEX)) return true;
	if (typeid(T)==typeid(DBLCOMPLEX)) return true;
	return false;
}

#endif