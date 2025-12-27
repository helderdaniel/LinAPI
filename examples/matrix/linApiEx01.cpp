/* LinAPI C++ demo
   @2009  hdaniel, mmadeira
*/


/* This example do not address type conversion, between Sca/Lapack data types
	
   DOUBLE, SINGLE, COMPLEX, ....

	Because tihs support is not yet available.
	One way to do implement that support can be defining specialized
	constructors for type convertion like:

	LinMatrix<DOUBLE>(SINGLE Op)

	which creates a new matrix LinMatrix<Double> object
	and copies Op._data to this._data

	Then C++ type convertion constructor can be called automatically
	in operations like:

	LinMatrix<DOUBLE> A
	LinMatrix<SINGLE> B

	LinMatrix<DOUBLE> C = A + B

*/

  
   
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <linexception.h>
using namespace std;

//LinAPI headers
#include <linapi.h>
#include <linmatrix.h>
#include <linvector.h>
#include <linsolver.h>

const int _CHECKCOMM_=1; //If >0 Check Blacs process grid
const int _CHECKTYPE_=1; //If >0 Print Scalapack types size in bytes
const int _CHECKDIST_=1; //If >0 Check Operator2D matrix distribution

#define _EXCEPTION_

#ifdef _EXCEPTION_
int mainThread (int argc, char **argv, LinApi& la) {
#else
int main (int argc, char ** argv) {
#endif
  
	//Blacs blacs& = la.blacs();	//is the same as:
	Blacs& blacs = Blacs::instance();
	if (_CHECKCOMM_>0) BlacsCheck::comm(cout); 	//Check Blacs process grid

   //define operators
   /* If Blacs singleton not defined, first call to Operator2D or its descendants
      will init singleton Blacs as a square process grid
    */
	//typedef SINGLE type;
	typedef DOUBLE type;
	//typedef COMPLEX  type; //prints elements in form (3,0) ... (5,0)
	//typedef DBLCOMPLEX  type; //prints elements in form (3,0) ... (5,0)
	//typedef bool type; //must throw Not supported class exception

	if (_CHECKTYPE_) 
		if (blacs.isRoot()) {
			cout << "Size of INTEGER (bytes)=    " << sizeof(INTEGER) << endl;
			cout << "Size of SINGLE (bytes)=     " << sizeof(SINGLE) << endl;
			cout << "Size of DOUBLE (bytes)=     " << sizeof(DOUBLE) << endl;
			cout << "Size of COMPLEX (bytes)=    " << sizeof(COMPLEX) << endl;
			cout << "Size of DBLCOMPLEX (bytes)= " << sizeof(DBLCOMPLEX) << endl;
		}

   type Adata[]={
    3, 1, 3,
    1, 5, 9,
    2, 6, 5
    };

	if (_CHECKDIST_) {
		LinMatrix<type> X(3, 3, Adata, "X");
		printf ("Proc=%d has rxc= %dx%d\n", blacs.grid().processId(),
			X.thisProcOperatorRows(), X.thisProcOperatorCols());

		type n=X.get(2,2); //C matrix have zero based indexes
		X.set(2,2, ((type)-2)*n); //-2*n complex<T> do not implement int*complex<> ?? must explicit convert -2 to complex
		cout << X; 
		X.print(cout, true); //X Transpose 

		LinMatrix<type> Z(3, 3, "Z"); //ident = false is zero operator
		Z.print(cout);

		LinMatrix<type> I0(3, 3, "I0", true); //ident = true is unity operator
		LinMatrix<type> I1; //I1=I0; uses copy constructor
		if (blacs.isRoot()) cout << I1.name() << ": " << I1.rows() << 'x' << I1.cols() << endl;
		I1=I0; //here operator = is used
		I0.set(0,0,-1.5);
		I0.print(cout);
		I1.print(cout);

		LinMatrix<type> I2(3, 4, "I2", true, 1, 2); //ident = true is unity operator
		LinMatrix<type> I3(I2);
		I3.print(cout);
	}

	LinMatrix<type> A(3, 3, Adata, "A"); 
	LinMatrix<type> I3(3, 3, "I3", true);

	type xdata[] = { -1, -1, 1 };
	LinVector<type> x(3, xdata, "x", 1, 1);

	type ydata[] = { 1, 2, 3 };
	LinVector<type> y(3, ydata, "y", 2, 2);

	LinVector<type> i(3, "i3", true);

	//print operators
	A.print(cout);
	y.print(cout);
	x.print(cout);

  /* Level 1 PBLAS:vector operations
     man (P)SDAXPY, DAXPY, ... (Double A*X+Y)
     (P)xAXPY perform:

     y = alpha * x + y
   */
	//PBLAS::Level1::add(2.5,x,y);
	//const Operator2D<type>& z=x+y; //OK but must be declared a reference since Operator2D<T> is abstract
	{	LinVector<type> z=x+y; //x+=y; // also valid but modifies x
		z.print(cout);
	} //z declared on a block to be freed at the end

	//2.0*x+y
	{	LinVector<type> z(y);
		x.axpy(2.0,z);
		z.print(cout);
	}

	/* Level 2 BLAS: matrix-vector operations
      man SGEMV, DGEMV, .... (Double GEneral Matrix Vector)
      xGEMV perform:
            
       y = alpha*A*x + beta*y or
       y = alpha*A'*x + beta*y
    */
	// y = A*x
/*  char TRANSA='N'; 
  DOUBLE BETA=0.0;
  pdgemv_(&TRANSA, &M, &N, &ALPHA, Al, &subRow, &subCol, DESCA, xl, &subRow, 
          &subCol, DESCX, &INCX, &BETA, yl, &subRow, &subCol, DESCY, &INCY);
  Cpdgemr2d(M, NRHS, yl, 1, 1, DESCY, y, 1, 1, DESCYr, DESCY[1]);
  if (thisProcRoot) printm ("y=A.x", y, 3);
*/            
#ifdef _LAPACK_
	//y=A(:,1)+...+A(:,3)+y

	/* Level 3 BLAS: matrix-matrix operations
     man SGEMM, DGEMM, .... (Double GEneral Matrix Matrix)
     xGEMM perform:

    C = alpha*op(A)*op(B) + beta*C, where op(X) can be the tranposed matrix
    specified by 2nd and 3rd parameters
  */      

	{ // C = A * A
		LinMatrix<type> B(A);
		B *= A;
		B.name("B=A; B*=A");
		B.print(cout);
		LinMatrix<type> C=A*A;
		C.name("C=B*A");
		C.print(cout);
	}
	{ // x=A*y
		LinMatrix<type> C=*((LinMatrix<type>*) &(A*y));	//When 2nd arg is LinVector, operator return
		C.name("C=A*y");											//LinVector. It can  be cast to a LinMatrix
		C.print(cout);											//only thru a pointer (why?)
		//this is the only way to convert an LinVector to an  LinMatrix (both descendants of Operator2D)
		//or define a convert operator on both classes

		LinVector<type> x =A*y;
		x.name("x=A*y");
		x.print(cout);
	}

	{ // C = A + A
		LinMatrix<type> C = A+A;
		C.name("C=A+A"); //rename it or name will be A(clone) atributed by copy
						     //constructor/assignment operator used in operation above
		C.print(cout);
	}
#endif

	/* Sca/Lapack linear system solver for general matrices:
   	man SGESV, DGESV, ... (Double GEneral SolVer)

   	y = A.x <=> x = inv(A) * y <=> x' = 0.244, 0.477, -0.070
	*/
	{  LinVector<type> z = LinSolver::gesv(A, y);
	  	z.name("y = A.x");
		z.print(cout);
	}

	//if (Blacs::instance().processId()==1)
	//throw 20;
	//throw LinException();
	//throw LinTypeNotSupportedException<type, LinMatrix<type> > ("method");
	//LinMatrix<char> xx; //not supported type exception

	return EXIT_SUCCESS;
}  

#ifdef _EXCEPTION_
int main (int argc, char ** argv) {
	LinApi& la=LinApi::square();
	int r=la.run(mainThread, argc, argv);

//Stop before glibc memory error on scalapck version (exception handling messaging has bugs
//  . run scast 7 should show error
//getchar(); 

	return r;
}
#endif