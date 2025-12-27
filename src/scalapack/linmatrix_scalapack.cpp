/* LinAPI matrix class v1.0
	Scalapack version
   @2009 hdaniel mmadeira
 */

#include "../linmatrix.h"
#include "../linexception.h"

#include <iostream>
#include <iomanip>
#include <sstream>

//Fix no export keyword for separate templates in *. h and *.cpp
//This allows compiler that implement export keyword or not to
//compile this source code.
//However it can be removed when support arrives at gcc (remove also from *.h)
//and remove this filename from exclusion (with patsubst) in makefile 
//so that makefile compiles all *.cpp on this dir 
#ifndef _LINMATRIX_SCALAPACK_
	#define _LINMATRIX_SCALAPACK_

	#ifndef USE_EXPORT_KEYWORD
   #define export 
#endif


/* External SCALAPACK functions prototypes
*/
extern "C" {
	//Scalapack tools (Not used)
	int indxg2p_ (const int* INDXGLOB, const int* NB, const int* IPROC, const int* ISRCPROC, const int* NPROCS );
	int indxg2l_(const int* INDXGLOB, const int* NB, const int* IPROC, const int* ISRCPROC, const int* NPROCS);
 	int indxl2g_ (const int* INDXLOC, const int* NB, const int* IPROC, const int* ISRCPROC, const int* NPROCS );

	//Level 2 PBLAS

	//Level 3 PBLAS
}

/* Wrapper functions
 */

// End wrapper functions

/* Global matrix index descodification to local processor and address
   as in Scalapack tools: INDXG2P and INDXG2L

	The principle of operation:

		(prior to calling this method all elements have been assigned value zero)

		In every process:
			For evey elements in the main diagonal
				if it is alocated on this process its value is turned to 1.

	This forces checking if all main diagonal elements are allocated on every process,
	even the elements that are known not to be on local process (matrix is distributed),
	using indxg2p_ and indxg2l or similar which have linear assimptotic complexity

	Alternatively it can be checked, for all the elements allocated on a process, if they
	are on main diagonal, using indxl2g or similar which also have linear assimptotic
	complexity.

	This last technique involves checking all matrix elements (however distributed by a
	process grid),	while the implemnented technique only checks all main diagonal elements
	(on every processor)

	For matrice operators larger than the process grid, which is the common case
	the implemeted technique should be faster:

		main diagonal elements to be checked on every processor for an mxn matrix:

			min (m,n)

		average of all the elements to be checked on a processor of the grid:

			(mxn)/ nº process on the grid

	assuming a square grid:

		n		(implemented technique is linear)

		n^2/nº process   (the alternate technique is quadratic)
 */
export template <class T>
void LinMatrix<T>::identity() {
const int nprows=this->blacs.grid().gridRows(); //for readability
const int npcols=this->blacs.grid().gridCols();	
const int thisProcRow=this->blacs.grid().processRow(); 
const int thisProcCol=this->blacs.grid().processCol(); 
const int MB = this->_rowsBlockFactor; 
const int NB = this->_colsBlockFactor; 
const int MBxNPROWS = MB*nprows; //calculate only once (optimize perf)
const int NBxNPCOLS = NB*npcols;	//calculate only once (optimize perf)

	//Only diagonal needes initialization, all other elements are already 0 in base constructor
	for (int i=0; i<this->_rows; ++i)
		for (int j=0; j<this->_cols; ++j)
			if (i==j)
				if(((i / MB) % nprows)==thisProcRow) //proc that has global row i
					if(((j / NB) % npcols)==thisProcCol) { //proc that has global col j
						int r = MB * (i / MBxNPROWS) + (i % MB);  //r = local row of global row i
						int c = NB * (j / NBxNPCOLS) + (i % NB);  //c = local col of global col j
						this->_data [c*this->_localRows+r] = (T) 1;
					}
}


/* Linear algebra
 */
export template <class T>
Operator2D<T>& LinMatrix<T>::operator += (const Operator2D<T>& op) {
	return *this;
}

#endif