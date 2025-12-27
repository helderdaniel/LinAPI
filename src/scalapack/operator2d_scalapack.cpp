/* LinAPI 2D operators base class v1.0
   Scalapack version
   @2009 hdaniel mmadeira
 */


//There is a bug probably here!!!
//*** glibc detected *** ./linApiEx01Scalapackst: corrupted double-linked list: 0x0811b520 ***

#include "../operator2d.h"
#include "../linexception.h"

#include <iomanip>
/*
#include <cmath>
#include <cctype>
#include <cstdlib>
*/

//Fix no export keyword for separate templates in *. h and *.cpp
//This allows compiler that implement export keyword or not to
//compile this source code.
//However it can be removed when support arrives at gcc (remove also from *.h)
//and remove this filename from exclusion (with patsubst) in makefile 
//so that makefile compiles all *.cpp on this dir 
#ifndef _OPERATOR2D_SCALAPACK_
	#define _OPERATOR2D_SCALAPACK_

	#ifndef USE_EXPORT_KEYWORD
   	#define export 
	#endif

static const int descCTXT=1; //index of process grid context store on descriptor

/* External SCALAPACK functions prototypes
*/
extern "C" {
	/* Scalapack Tools wrappers
 	 */
   //Return number of row or columns (RC) of operator in this process
	int numroc_(const INTEGER *globalOperatorRC, const INTEGER *globalOperatorRCBlockSize,
					const INTEGER *thisProcGridRC, const INTEGER *SrcProcGridRC,
					const INTEGER *gridRC);
	//Init scalapack distributed matrix descriptor
	int descinit_ (const INTEGER *DESC, const INTEGER *Mp, const INTEGER *Np,
						const INTEGER* MBp, const INTEGER *NBp, const INTEGER *IRSRCp,
						const INTEGER *ICSRCp, const INTEGER *ICTXTp, const INTEGER *LLDp,
						const INTEGER *INFOp);
	//Get start row and col of Global matrix alocated in this process
	void infog1l_ (const INTEGER* GINDXp, const INTEGER* NBp, const INTEGER* NPROCSp,
						const INTEGER* MYROCp, const INTEGER* ISRCp, const INTEGER* LINDXp,
						const INTEGER* ROCSRCp);
	void infog2l_ (const INTEGER* GRINDXp, const INTEGER* GCINDXp, const INTEGER* DESC,
						const INTEGER* NPROWp, const INTEGER* NPCOLp, const INTEGER* MYROWp,
						const INTEGER* MYCOLp, const INTEGER* LRINDXp, const INTEGER* LCINDXp,
						const INTEGER* RSRCp, const INTEGER* CSRCp);

	//get/set operator single elements
	void pielset_(INTEGER* A, INTEGER* IAp, INTEGER* JAp, INTEGER* DESCA, INTEGER* ALPHAp);
	void pselset_(SINGLE* A, INTEGER* IAp, INTEGER* JAp, INTEGER* DESCA, SINGLE* ALPHAp);
	void pdelset_(DOUBLE* A, INTEGER* IAp, INTEGER* JAp, INTEGER* DESCA, DOUBLE* ALPHAp);
	void pcelset_(COMPLEX* A, INTEGER* IAp, INTEGER* JAp, INTEGER* DESCA, COMPLEX* ALPHAp);
	void pzelset_(DBLCOMPLEX* A, INTEGER* IAp, INTEGER* JAp, INTEGER* DESCA, DBLCOMPLEX* ALPHAp);

	void pielget_(CHAR *SCOPEp, CHAR *TOPp, INTEGER *ALPHAp, INTEGER *A,
   	           INTEGER* IAp, INTEGER* JAp, INTEGER *DESCA);
	void pselget_(CHAR *SCOPEp, CHAR *TOPp, SINGLE *ALPHAp, SINGLE *A,
   	           INTEGER* IAp, INTEGER* JAp, INTEGER *DESCA);
	void pdelget_(CHAR *SCOPEp, CHAR *TOPp, DOUBLE *ALPHAp, DOUBLE *A,
   	           INTEGER* IAp, INTEGER* JAp, INTEGER *DESCA);
	void pcelget_(CHAR *SCOPEp, CHAR *TOPp, COMPLEX *ALPHAp, COMPLEX *A,
   	           INTEGER* IAp, INTEGER* JAp, INTEGER *DESCA);
	void pzelget_(CHAR *SCOPEp, CHAR *TOPp, DBLCOMPLEX *ALPHAp, DBLCOMPLEX *A,
   	           INTEGER* IAp, INTEGER* JAp, INTEGER *DESCA);

	//Operator redistribute/copy
	void Cpigemr2d(INTEGER m, INTEGER n, INTEGER *A, INTEGER IA, INTEGER JA, INTEGER *DESCA,
						INTEGER *B, INTEGER IB, INTEGER JB, INTEGER *descB, INTEGER gcontext);
	void Cpsgemr2d(INTEGER m, INTEGER n, SINGLE *A, INTEGER IA, INTEGER JA, INTEGER *DESCA,
						SINGLE *B, INTEGER IB, INTEGER JB, INTEGER *descB, INTEGER gcontext);
	void Cpdgemr2d(INTEGER m, INTEGER n, DOUBLE *A, INTEGER IA, INTEGER JA, INTEGER *DESCA,
						DOUBLE *B, INTEGER IB, INTEGER JB, INTEGER *descB, INTEGER gcontext);
	void Cpcgemr2d(INTEGER m, INTEGER n, COMPLEX *A, INTEGER IA, INTEGER JA, INTEGER *DESCA,
						COMPLEX *B, INTEGER IB, INTEGER JB, INTEGER *descB, INTEGER gcontext);
	void Cpzgemr2d(INTEGER m, INTEGER n, DBLCOMPLEX *A, INTEGER IA, INTEGER JA, INTEGER *DESCA,
						DBLCOMPLEX *B, INTEGER IB, INTEGER JB, INTEGER *descB, INTEGER gcontext);

}


export template <class T>
Operator2Ddesc* Operator2D<T>::initDesc (int MBp, int NBp, int ICTXTp, int LLDp) {
INTEGER info;
GridCoord rootProc=blacs.processRowColId(Blacs::rootId);
const INTEGER rootProcRow=rootProc.row();
const INTEGER rootProcCol=rootProc.col();
Operator2Ddesc *desc = new Operator2Ddesc[Operator2DdescLen];

	descinit_ (desc, (INTEGER *)&_rows, (INTEGER *)&_cols, (INTEGER *)&MBp, (INTEGER *)&NBp,
					&rootProcRow, &rootProcCol, (INTEGER *)&ICTXTp, (INTEGER *)&LLDp, &info);
	if (info!=0) {
		stringstream ss;
		ss << "Error on DESCINIT (Scalapack):" << -info;
		ss << ". argument, when initializing operator: " << _name << endl;
		throw out_of_range(ss.str());
	}
	return desc;
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
	_rowsBlockFactor=mb; _colsBlockFactor=nb;
	if (name==NULL) {
		stringstream ss;
		ss << opBaseName << setw(opDigits) << opCounter++;
		_name = ss.str();
	}
	else _name = string(name);

	//Local operator and process grid info
	INTEGER gridRows = blacs.grid().gridRows();
	INTEGER thisProcRow = blacs.grid().processRow();
	/* Need to pass static var value to instance var, because g++-4.3 can not
		resolve address of static member:
				&sourceProcRow
		if value defined inside class:
				static const int sourceProcRow=0;
		But if value of static var defined outside class:
				template <class T> int Operator2D<T>::sourceProcRow=0;
		 gcc finds it (gcc BUG?):
	*/
	INTEGER srcProcRow = sourceProcRow;
	_localRows = numroc_ ((INTEGER *)&_rows, (INTEGER *)&_rowsBlockFactor, &thisProcRow, &srcProcRow, &gridRows);

	INTEGER gridCols=blacs.grid().gridCols();
	INTEGER thisProcCol=blacs.grid().processCol();
	INTEGER srcProcCol = sourceProcCol;
	_localCols = numroc_ ((INTEGER *)&_cols, (INTEGER *)&_colsBlockFactor, &thisProcCol, &srcProcCol, &gridCols);

	//Scalapack matrix descriptors
	//process grid descriptor
	int lld = max/*<int>*/(fMBase, _localRows);
	int gridContext = blacs.grid().context();
	_gridDesc = initDesc (_rowsBlockFactor, _colsBlockFactor, gridContext, lld);

	//root collector descriptor
	int rootContext = blacs.root().context();
	if (blacs.isRoot()) _rootDesc = initDesc (_rows, _cols, rootContext, _rows);
	else  { //if not root processor, grid context is invalid
		_rootDesc = new Operator2Ddesc[Operator2DdescLen];
		_rootDesc[descCTXT]=NOGRIDCONTEXT;
	}

	//Init local operator data area
	_localElements=_localRows*_localCols;
	_data = new T[_localElements];

}

export template <class T> Operator2D<T>::
Operator2D (int m, int n, T* data, const char* name, int mb, int nb) : blacs(Blacs::instance()) {

	init (m, n, name, mb, nb);
	//Distribute operator
	//making a version of CpXgemr2d where rows and cols indexes are exchanged
	//(maybe not as easy as it sounds, since matrix is distributed on a process grid)
	//will eliminate the need of a transpose step
	//Transposed is required since C allocates matrices by rows and Fortran by columns
	T* trans = new T[m*n];  //or _rows*_cols inited in int function
	linearTranspose(data, trans);
	distribute(trans);
	delete[] trans;
}


/* Distribute operators on process grid and collect them on root process
   Can not collect or distribute a section of a matrix, only the whole operator
 */
export template <class T> void Operator2D<T>::distribute(T* globalData) const {
	throw LinTypeNotSupportedException<T, typeof(*this)> ("distribute");
}
/* Not specialiing INTEGER (using Cpigemr2d) avoid support for Operator2D<INTEGER>,
	since Sca/Lapack has no math functions that support integers
*/
export template <> void Operator2D<SINGLE>::distribute(SINGLE* globalData) const {
	Cpsgemr2d(_rows, _cols, globalData, fMBase, fMBase, _rootDesc, _data, fMBase, fMBase, _gridDesc, _gridDesc[descCTXT]);
}
export template <> void Operator2D<DOUBLE>::distribute(DOUBLE* globalData) const {
	Cpdgemr2d(_rows, _cols, globalData, fMBase, fMBase, _rootDesc, _data, fMBase, fMBase, _gridDesc, _gridDesc[descCTXT]);
}
export template <> void Operator2D<COMPLEX>::distribute(COMPLEX* globalData) const {
	Cpcgemr2d(_rows, _cols, globalData, fMBase, fMBase, _rootDesc, _data, fMBase, fMBase, _gridDesc, _gridDesc[descCTXT]);
}
export template <> void Operator2D<DBLCOMPLEX>::distribute(DBLCOMPLEX* globalData) const {
	Cpzgemr2d(_rows, _cols, globalData, fMBase, fMBase, _rootDesc, _data, fMBase, fMBase, _gridDesc, _gridDesc[descCTXT]);
}

export template <class T> void Operator2D<T>::collect(T* globalData) const {
	throw LinTypeNotSupportedException<T, typeof(*this)> ("collect");
}
/* Not specialiing INTEGER (using Cpigemr2d) avoid support for Operator2D<INTEGER>,
	since Sca/Lapack has no math functions that support integers
*/
export template <> void Operator2D<SINGLE>::collect(SINGLE *globalData) const {
	Cpsgemr2d(_rows, _cols, _data, fMBase, fMBase, _gridDesc, globalData, fMBase, fMBase, _rootDesc, _gridDesc[descCTXT]);
}
export template <> void Operator2D<DOUBLE>::collect(DOUBLE *globalData) const {
	Cpdgemr2d(_rows, _cols, _data, fMBase, fMBase, _gridDesc, globalData, fMBase, fMBase, _rootDesc, _gridDesc[descCTXT]);
}
export template <> void Operator2D<COMPLEX>::collect(COMPLEX *globalData) const {
	Cpcgemr2d(_rows, _cols, _data, fMBase, fMBase, _gridDesc, globalData, fMBase, fMBase, _rootDesc, _gridDesc[descCTXT]);
}
export template <> void Operator2D<DBLCOMPLEX>::collect(DBLCOMPLEX *globalData) const {
	Cpzgemr2d(_rows, _cols, _data, fMBase, fMBase, _gridDesc, globalData, fMBase, fMBase, _rootDesc, _gridDesc[descCTXT]);
}


/* Print operator 
 */
export template <class T>
ostream& Operator2D<T>::toStream(ostream &os, const bool transpose) const {
T* tmp = new T[_rows*_cols];

	collect(tmp);
	if (blacs.isRoot()) linearToStream(os, tmp, transpose);
	delete[] tmp;
	return os;
}


/* Access individual elements

   Rows and columns are incremented, so that matrix will be indexed
	one based when calling Fortran routines
 */
export template <class T> void Operator2D<T>::set(int r, int c, T alpha) {
	//this should not happen since Operator2D constructor prevents
	//creation of unsupported types
	throw LinTypeNotSupportedException<T, typeof(*this)> ("set");
}
/* Not defining set<INTEGER> (using pielset_) since this type of matrices
	have no support in Sca/Lapack
 */
export template <> void Operator2D<SINGLE>::set(int r, int c, SINGLE alpha) {
	pselset_ (_data, (INTEGER*)&(++r), (INTEGER*)&(++c), _gridDesc, &alpha);
}
export template <> void Operator2D<DOUBLE>::set(int r, int c, DOUBLE alpha) {
	pdelset_ (_data, (INTEGER*)&(++r), (INTEGER*)&(++c), _gridDesc, &alpha);
}
export template <> void Operator2D<COMPLEX>::set(int r, int c, COMPLEX alpha) {
	pcelset_ (_data, (INTEGER*)&(++r), (INTEGER*)&(++c), _gridDesc, &alpha);
}
export template <> void Operator2D<DBLCOMPLEX>::set(int r, int c, DBLCOMPLEX alpha) {
	pzelset_ (_data, (INTEGER*)&(++r), (INTEGER*)&(++c), _gridDesc, &alpha);
}


/* p?elget_ wrapper
 */
template <class T> static
void elgetw (CHAR SCOPEp, CHAR TOPp, T *ALPHAp, T *A, INTEGER IAp,
					INTEGER JAp, INTEGER *DESCA) {
	//this should not happen since Operator2D contructor prevents
	//creation of unsupported types
	throw LinTypeNotSupportedException<T, Operator2D<T> > ("get");
}

/* Not defining set<INTEGER> (using pielget_) since this type of matrices
	have no support in Sca/Lapack
 */
template <> 
void elgetw<SINGLE> (CHAR SCOPEp, CHAR TOPp, SINGLE *ALPHAp, SINGLE *A,
					 			INTEGER IAp, INTEGER JAp, INTEGER *DESCA) {
	pselget_(&SCOPEp, &TOPp, ALPHAp, A, &IAp, &JAp, DESCA);
}
template <> 
void elgetw<DOUBLE> (CHAR SCOPEp, CHAR TOPp, DOUBLE *ALPHAp, DOUBLE *A,
					 			INTEGER IAp, INTEGER JAp, INTEGER *DESCA) {
	pdelget_(&SCOPEp, &TOPp, ALPHAp, A, &IAp, &JAp, DESCA);
}
template <> 
void elgetw<COMPLEX> (CHAR SCOPEp, CHAR TOPp, COMPLEX *ALPHAp, COMPLEX *A,
					 			INTEGER IAp, INTEGER JAp, INTEGER *DESCA) {
	pcelget_(&SCOPEp, &TOPp, ALPHAp, A, &IAp, &JAp, DESCA);
}
template <> 
void elgetw<DBLCOMPLEX> (CHAR SCOPEp, CHAR TOPp, DBLCOMPLEX *ALPHAp, DBLCOMPLEX *A,
					 			INTEGER IAp, INTEGER JAp, INTEGER *DESCA) {
	pzelget_(&SCOPEp, &TOPp, ALPHAp, A, &IAp, &JAp, DESCA);
}


export template <class T> T Operator2D<T>::get(int r, int c) const {
CHAR SCOPE='A'; //Get element considering all processes in grid
CHAR TOP=' ';   //Default system dependant broadcast
T alpha;
	elgetw (SCOPE, TOP, &alpha, _data, (INTEGER)(++r), (INTEGER)(++c), _gridDesc);
	return alpha;
}

#endif