/* BLACS C++ interface
   Scalapack MPI version v1.0
   @2009 hdaniel mmadeira
 */

#include "../blacs.h"
#include "mpiobj.h"
#include "../linexception.h"
#include <cmath>
#include <cctype>
#include <cstdlib>

const string VERSION("Blacs v1.0 Scalapack");

/* External BLACS C interface functions prototypes
   If not defined it gives warnings with mpicc -Wall (gcc ANSI C)
   If not defined it gives errors with mpicxx (g++ C++)   
   Anyway prototypes are advisable to check variables type on compile
*/
extern "C" {
	int Cblacs_pinfo (int *thisProcNum, int *nprocs);
	//void Cblacs_setup (int *thisProcNum, int *nprocs);
	void Cblacs_get (int context, int what, int *val);
	void Cblacs_gridinit (int *context, char *rowColMajor, int prow, int pcol);
	void Cblacs_gridinfo (int context, int *npRow, int *npCol, int *thisPRow, int *thisCol);
	int Cblacs_pnum (int context, int prow, int pcol);
	int Cblacs_pcoord (int context, int pnum, int *prow, int *pcol);
	void Cblacs_exit (int doneflag);
	void Cblacs_gridexit (int context);
	void Cblacs_barrier(int context, char* scope);

	/* 4th argument of function 
   
     C?gesd2d  (send functions)
     C?gebs2d
     C?trsd2d
     C?trbs2d

	  C?gerv2d  (receive functions)
     C?gebr2d
     C?trrv2d
     C?trbr2d     
     
   depends on Precision:
   
	Single precision real (?==s) or complex (?==c) -> float *
	Double precision real (?==d) or complex (?==z) -> double *
	integer (?==i) -> int * (the one used above)
	*/
   //only C?gesd2d  (send functions) and C?gerv2d  (receive functions) currently implemented:
	void Csgesd2d(int context, int M, int N, SINGLE *A, int LDA, int RDEST, int CDEST);
	void Cdgesd2d(int context, int M, int N, DOUBLE *A, int LDA, int RDEST, int CDEST);
	void Ccgesd2d(int context, int M, int N, COMPLEX *A, int LDA, int RDEST, int CDEST);
	void Czgesd2d(int context, int M, int N, DBLCOMPLEX *A, int LDA, int RDEST, int CDEST);
	void Cigesd2d(int context, int M, int N, INTEGER *A, int LDA, int RDEST, int CDEST);

	void Csgerv2d(int context, int M, int N, SINGLE *A, int LDA, int RSRC, int CSRC);
   void Cdgerv2d(int context, int M, int N, DOUBLE *A, int LDA, int RSRC, int CSRC);
   void Ccgerv2d(int context, int M, int N, COMPLEX *A, int LDA, int RSRC, int CSRC);
	void Czgerv2d(int context, int M, int N, DBLCOMPLEX *A, int LDA, int RSRC, int CSRC);
	void Cigerv2d(int context, int M, int N, INTEGER *A, int LDA, int RSRC, int CSRC);

//distribute matrix PUT it here or on blacs
/*void Cpsgemr2d(int m, int n, SINGLE *A, int IA, int JA, int *DESCA,
   	            DOUBLE *B, int IB, int JB, int *descB, int gcontext);
  void Cpdgemr2d(int m, int n, DOUBLE *A, int IA, int JA, int *DESCA,
   	            DOUBLE *B, int IB, int JB, int *descB, int gcontext);

   //Collect matrix
*/

}



/* ProcessGrid class implementation
 */
void ProcessGrid::initGrid (int rows, int cols) {
   _context=NOGRIDCONTEXT;
   //Get default system context (1st arg==0)
	Cblacs_get(0, 0, &_context);

   //Create process grid row major ordering
	Cblacs_gridinit(&_context, (char * const) "Row", _gridRows=rows, _gridCols=cols);
   Cblacs_gridinfo(_context, &_gridRows, &_gridCols, &thisPRow, &thisPCol);

	//Get this process grid ID (an int) from this process grid coordinates
	//only if _context is valid
	if (_context>NOGRIDCONTEXT) thisProc = Cblacs_pnum (_context, thisPRow, thisPCol);
}

/* Default constructor, init _context with invalid value
 */
ProcessGrid::ProcessGrid () { _context=NOGRIDCONTEXT; }

//Create grid with rows and cols suplied
ProcessGrid::ProcessGrid (int rows, int cols) { initGrid (rows, cols); }

/* Since in C++ destructors executes body first and then call destructors for
   member variables, destructing Blacs object will call ProcessGrid destructors
   after CBlacs_exit(0); which is an errors since MPI is already finalized.

   Got to find a way to fix this!
   Temporary solution uses function clear(), which terminates GridProcess
   with Blacs_gridexit() and set _context variable to an invalid value, indicating no grid context defined.
   This solves Blacs object destruction since clear() is called before Blacs destructor, destructing ProcessGrid contexts and indicating it on _context==-1,
	so that after call to CBlacs_exit(0), ProcessGrid destructor do not try to free
	context again.
 */
void ProcessGrid::clear () {
	if (_context>NOGRIDCONTEXT) {
		Cblacs_gridexit(_context);
		_context=NOGRIDCONTEXT;
	}
}

ProcessGrid::~ProcessGrid () { clear(); }

//Return process processId coordinates in process grid 
GridCoord ProcessGrid::processRowColId(int processId) {
int r, c;
	Cblacs_pcoord(_context, processId, &r, &c);
	return GridCoord(r, c);
}

int ProcessGrid::processAtRowCol(const GridCoord& procCoord) {
	return Cblacs_pnum (_context, procCoord.row(), procCoord.col());
}

void ProcessGrid::barrier() {
char scope='A';
	Cblacs_barrier(_context, &scope);
}



/* Blacs singleton class implementation
 */
string Blacs::_version=VERSION;
const int Blacs::MPIprocs = Mpi::instance().procCount();
const int Blacs::rootId = Mpi::rootId;

Blacs& Blacs::custom(int rows, int cols) {
	static Blacs b(rows, cols);
	return b;
}

Blacs& Blacs::linear() { return custom(MPIprocs, 1); }

Blacs& Blacs::square() {
  	int rows = (int) sqrt((double) MPIprocs );
	int cols = MPIprocs / rows;
	Blacs &b =custom(rows, cols);

	//check if this process is out of grid and end process
   ProcessGrid &g = b.grid();
	if ((g.processRow() >= g.gridRows()) || (g.processCol() >= g.gridCols())) {
		printf ("MPI proc %d out of GRID\n", Mpi::instance().procId());
   	exit(EXIT_SUCCESS); //terminating program, calls all destructors including Blacs
	}

	return b;
}

Blacs::Blacs (int rows, int cols) : _root(1,1), _grid (rows, cols) { _procCount = rows*cols; }

Blacs::~Blacs () {
   root().clear();
	grid().clear();
	Cblacs_exit(0);
}

/* Blacs general matrix point to point send and receive functions

	All Sca/Lapack types supported have a specialization,	so that if user
	tries to instantiate an	unsupported type, compiler issues an error
	However this should not happen since Operator2D contructor prevents
	creation of unsupported types

	Not needed a general template ges2d<class T> since it is already
	declared as a member function
 */
//Blacs general matrix point to point send functions
template <>
void Blacs::gesd2d<INTEGER>(int M, int N, INTEGER *A, int LDA, int RDEST, int CDEST) {
	Cigesd2d(grid().context(), M, N, A, LDA, RDEST, CDEST);
}
template <>
void Blacs::gesd2d<SINGLE>(int M, int N, SINGLE *A, int LDA, int RDEST, int CDEST) {
	Csgesd2d(grid().context(), M, N, A, LDA, RDEST, CDEST);
}
template <>
void Blacs::gesd2d<DOUBLE>(int M, int N, DOUBLE *A, int LDA, int RDEST, int CDEST) {
	Cdgesd2d(grid().context(), M, N, A, LDA, RDEST, CDEST);
}
template <>
void Blacs::gesd2d<COMPLEX>(int M, int N, COMPLEX *A, int LDA, int RDEST, int CDEST) {
	Ccgesd2d(grid().context(), M, N, A, LDA, RDEST, CDEST);
}
template <>
void Blacs::gesd2d<DBLCOMPLEX>(int M, int N, DBLCOMPLEX *A, int LDA, int RDEST, int CDEST) {
	Czgesd2d(grid().context(), M, N, A, LDA, RDEST, CDEST);
}


//Blacs general matrix point to point receive functions
template <>
void Blacs::gerv2d<INTEGER>(int M, int N, INTEGER *A, int LDA, int RSRC, int CSRC) {
	Cigerv2d(grid().context(), M, N, A, LDA, RSRC, CSRC);
}
template <>
void Blacs::gerv2d<SINGLE>(int M, int N, SINGLE *A, int LDA, int RSRC, int CSRC) {
	Csgerv2d(grid().context(), M, N, A, LDA, RSRC, CSRC);
}
template <>
void Blacs::gerv2d<DOUBLE>(int M, int N, DOUBLE *A, int LDA, int RSRC, int CSRC) {
	Cdgerv2d(grid().context(), M, N, A, LDA, RSRC, CSRC);
}
template <>
void Blacs::gerv2d<COMPLEX>(int M, int N, COMPLEX *A, int LDA, int RSRC, int CSRC) {
	Ccgerv2d(grid().context(), M, N, A, LDA, RSRC, CSRC);
}
template <>
void Blacs::gerv2d<DBLCOMPLEX>(int M, int N, DBLCOMPLEX *A, int LDA, int RSRC, int CSRC) {
	Czgerv2d(grid().context(), M, N, A, LDA, RSRC, CSRC);
}


//Low-level process signaling (using MPI)
void Blacs::signalProc(int p) { Mpi::instance().signalProc(p); };
void Blacs::waitProc(int p) { Mpi::instance().waitProc(p); };
int Blacs::waitAny() { return Mpi::instance().waitAny(); };

//Blacs grid may have less procs then specified in MPI if it is a
//forced close to square grid, so _procCount indicates number of Blacs procs
void Blacs::signalAll() { Mpi::instance().signalAll(_procCount); };
void Blacs::waitAll() { Mpi::instance().waitAll(_procCount); };
void Blacs::signalAllOthers() { Mpi::instance().signalAllOthers(_procCount); };
void Blacs::waitAllOthers() { Mpi::instance().waitAllOthers(_procCount); };

