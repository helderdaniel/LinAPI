/* LinAPI 2D operators base class v1.0
   @2009 hdaniel, mmadeira
 */

#ifndef _OPERATOR2D_H_
#define _OPERATOR2D_H_

#include "blacs.h"
#include "../linexception.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>

using namespace std;

typedef INTEGER Operator2Ddesc;

template <class T>
class Operator2D {

	friend class LinSolver;

	virtual void abstract()=0; //make this base class abstract
	void init(int m, int n, const char* name, int mb, int nb);
	//transpose 2D operator stored in linear addressing fashion
	void linearTranspose(T* original, T* trans); 
	//send to stream 2D operator elements stored in linear addressing fashion
	ostream& linearToStream(ostream &os, T* data, const bool transpose) const;

public:
	static const int sourceProcRow=0; 
	static const int sourceProcCol=0;

protected:
	static const int Operator2DdescLen=9;
	static const int opDigits=6;
	static const int defaultWidth=5;
	static const int defaultPrecision=3;
	static const int fMBase=1; //fortran matrix base index
	static const int cMBase=0; //C/C++ matrix base index
	static const int interleave; //interleave for interlaced operators
	static const int subRow; //subMatrix start row
	static const int subCol; //and col
	static const string opBaseName;
	static int opCounter;


	//print field width and precision digits
	int _width;
	int _precision;

	/* Blacs singleton reference can not be static variable because if so it is
      executed before program, definning default square grid, avoiding user to
      define its topology.
		However if not defined by user before first instanciation of operator,
      it is then defined as a square process grid.
    */
	Blacs& blacs;
	int _rowsBlockFactor, _colsBlockFactor;
	int _rows, _cols;
	int _localRows, _localCols; //local stands for this process
	int _localElements;
	string _name;
	T* _data;
protected:
	Operator2Ddesc *_gridDesc, *_rootDesc;
	Operator2Ddesc* initDesc (int MBp, int NBp, int ICTXT, int LLD);

	//misc
	void setAll(const T&);
	void swap(const Operator2D<T>& op); //swap op with *this

	//Distribute and collect operators
public:
void distribute(T* globalData) const;
	void collect(T* globalData) const;


public:
   /* Default const called when: Operator2D<?> A;, must leave object in a
	   consistent state, thus call init for a scalar 1x1 matrix with
	   When assigned something: A = B, previous object will be destructed by =operator
	   Note: reference instance variable must be initialized in constructor init list
	 */
	Operator2D () : blacs(Blacs::instance()) { init(1, 1, "Default to scalar", 1, 1); }
	Operator2D (int m, int n, const char* name=NULL, int mb=2,int nb=2);
	Operator2D (int m, int n, T* data, const char* name=NULL, int mb=2,int nb=2);
	Operator2D (const Operator2D<T>& op);
	virtual ~Operator2D() { delete[] _data; delete[] _gridDesc; delete[] _rootDesc; };

	int rows() const { return _rows; }
	int cols() const { return _cols; }
	int thisProcOperatorRows() const { return _localRows; }
	int thisProcOperatorCols() const { return _localCols; }
	string name() const { return _name; }
	void name(string name) { _name=name; }

	int width() const { return _width; }
	void width(int w) { _width = w; }
	int precision() const { return _precision; }
	void precision(int p) { _precision = p; }

	//set all elements to 0
	void reset() { setAll((T) 0); };

	//Access individual elements
	T get(int r, int c) const;
	void set(int r, int c, T alpha);

	//print
	ostream& toStream(ostream &os, const bool transpose=false) const;
	void print (ostream &os, bool transpose=false) const;
	template <class X>
	friend ostream& operator << (ostream& os, const Operator2D<X>& op);

	//linear algebra
	//Defines interface for Operator2D
	//all abstract functions must be implemented on concrete class
	virtual Operator2D<T>& operator += (const Operator2D<T>& op) = 0;
	//uncomment only if implemented on Linvector (dot product ?)
	//virtual Operator2D<T>& operator *= (const Operator2D<T>& op) = 0;
};



/* Common Scalapack / Lapack version implementation
 */
template <class T> int Operator2D<T>::opCounter=0;
template <class T> const string Operator2D<T>::opBaseName="Operator";
template <class T> const int Operator2D<T>::interleave=1;
template <class T> const int Operator2D<T>::subRow=1;
template <class T> const int Operator2D<T>::subCol=1;

//Generic array clone if not NULL
//could be in another general library
template <class X> X* cloneArray(X* src, const int& size) {
X* dst=NULL;

	if (src!=NULL) {
		dst = new X[size];
		for (int i=0; i<size; ++i)  dst[i]=src[i];
	}
	return dst;
}

//swap op with *this
template <class T> void Operator2D<T>::swap(const Operator2D<T>& op) {
char tmp;
char* thisObject = (char*) this;
char* opObject = (char*) &op;
const int size = sizeof(*this);

	for (int i=0; i<size; ++i) {
		tmp=thisObject[i];
		thisObject[i]=opObject[i];
		opObject[i]=tmp;
	}
}		


//Copy constructor
template <class T> Operator2D<T>::Operator2D(const Operator2D<T>& op) 
	//References
	: blacs(op.blacs) {

	//primitive types
	_width = op._width;
	_precision = op._precision;
	_rowsBlockFactor = op._rowsBlockFactor;
	_colsBlockFactor = op._colsBlockFactor;
	_rows = op._rows;
	_cols = op._cols;
	_localRows = op._localRows;
	_localCols = op._localCols;
	_localElements = op._localElements;

	//heap (dynamic) objects
	_name=string(op._name+"(clone)"); //add "clone" to indicate its a clone (?)
	_data = cloneArray (op._data, _localElements);
	_gridDesc = cloneArray (op._gridDesc, (int) Operator2DdescLen);
	_rootDesc = cloneArray (op._rootDesc, (int) Operator2DdescLen);
}


//Constructor 2
template <class T>Operator2D<T>::
Operator2D(int m, int n, const char* name, int mb,int nb) : blacs(Blacs::instance()) {
	init (m, n, name, mb, nb);
	setAll((T) 0);
}

//Set all elements to n (For complex type it works only with std::complex)
template <class T> void Operator2D<T>::setAll(const T& n) {
const int elements = _localRows*_localCols;
	for (int i=0; i<elements; ++i) _data [i] = n;
}

template <class T>
void Operator2D<T>::linearTranspose(T* original, T* trans) {
	for (int r=0; r<_rows; ++r)
		for (int c=0; c<_cols; ++c)
			trans[c*_rows+r]=original[r*_cols+c];
}

template <class T>
ostream& Operator2D<T>::linearToStream(ostream &os, T* data, const bool transpose) const {
	for (int i=0; i<_rows; ++i) {
			for (int j=0; j<_cols; ++j) {
				os << setw(4) << fixed << setw(_width) << setprecision(_precision);
				if (transpose) os << data [i*_cols+j] << ' ';
				else os << data [j*_rows+i] << ' '; //transpose on sending since F77 stores by columns and C by rows
			}
			os << '\n';
		}
	return os;
}

template <class T> void Operator2D<T>::print(ostream &os, bool transpose) const {
	//check if root not really needed on Lapack implementation
	if (Blacs::instance().isRoot()) os << _name << " =" << endl;
	toStream(os, transpose); //already check if it is root
	
}

template <class X>
ostream& operator << (ostream& os, const Operator2D<X>& op) { return op.toStream(os); }


//Fix no export keyword for separate templates in *. h and *.cpp
//This allows compiler that implement export keyword or not to
//compile this source code.
//However it can be removed when support arrives at gcc (remove also from *.cpp)
//Makefile define 2 diferent MACROS _SCALAPACK_ and _LAPACK_
//depending on version lapack/scalapack
#ifndef USE_EXPORT_KEYWORD
	#ifdef _SCALAPACK_
   	#include "scalapack/operator2d_scalapack.cpp"
	#endif
	#ifdef _LAPACK_
		#include "lapack/operator2d_lapack.cpp"
	#endif
#endif


#endif