/* BLACS C++ interface
   Lapack version (only one process on grid is considered)
   @2009 hdaniel mmadeira
 */

#include "../blacs.h"
#include "linexception.h"
#include <cstdlib>
#include <sstream>
#include <stdexcept>
using namespace std;

const string VERSION("Blacs v1.0 lapack");

/* Since there will be only on process this buffer will be used only when
   process 0 (0,0) send messages to itself.
   The general case of MEssage passing simulation requires a multiple
   buffer solution, one for each process. In this case a class Buffer
   is more appropriate
 */
void *comm_buffer=NULL;

/* Single rocess grid constants
 */
static const int GRIDROWS = 1;
static const int GRIDCOLS = 1;
static const int MPIPROCS = 1;
static const int ROOTID = 0;
static const int SINGLEPROCID = 0;
static const int SINGLEPROCROW = 0;
static const int SINGLEPROCCOL = 0;
static const int FIRSTCONTEXT = 0;
static int contextCount = FIRSTCONTEXT;


/* ProcessGrid class implementation
 */
void ProcessGrid::initGrid (int rows, int cols) {
   _context = contextCount ++; 	//atribute context ID and inc counter
	_gridRows = GRIDROWS;
	_gridCols = GRIDCOLS; 			//Only one process in grid (ignores rows, cols arguments)
	thisProc = SINGLEPROCID;		//Single process ID is 0
	thisPRow = SINGLEPROCCOL;
   thisPCol = SINGLEPROCCOL; 		//Single process is at (0,0) grid position
}

/* Default constructor, init _context with invalid value
 */
ProcessGrid::ProcessGrid () { _context=NOGRIDCONTEXT; }

/* Create grid with rows and cols suplied  (however initGrid ignore rows and cols
   values, assuming each equal to 1
 */
ProcessGrid::ProcessGrid (int rows, int cols) { initGrid (rows, cols); }

//Set context to invalid
void ProcessGrid::clear () { _context=NOGRIDCONTEXT; }
ProcessGrid::~ProcessGrid () { clear(); }

//Return process processId coordinates in process grid 
GridCoord ProcessGrid::processRowColId(int processId) {
	if (processId == SINGLEPROCID) return GridCoord(SINGLEPROCROW, SINGLEPROCCOL);
	else return GridCoord(NOTINGRID, NOTINGRID);
}

int ProcessGrid::processAtRowCol(const GridCoord& procCoord) {
	if (procCoord.row()==SINGLEPROCROW && procCoord.col()==SINGLEPROCCOL) return SINGLEPROCID;
	else return NOTINGRID;
}

//Since only one process exists barrier is always reached
void ProcessGrid::barrier() { }


/* Blacs singleton class implementation
 */
string Blacs::_version=VERSION;
const int Blacs::MPIprocs = MPIPROCS;
const int Blacs::rootId = ROOTID;

Blacs& Blacs::custom(int rows, int cols) {
	static Blacs b(rows, cols);
	return b;
}
Blacs& Blacs::linear() { return custom(MPIprocs, 1); } //MPIprocs=1 so equal to square
Blacs& Blacs::square() { return custom(1, 1); } 		 //only one process
Blacs::Blacs (int rows, int cols) : _root(1,1), _grid (rows, cols) { _procCount = MPIPROCS; }
Blacs::~Blacs () {
   root().clear();
	grid().clear();
}


void procInGrid(int procRow, int procCol) {
stringstream ss;
	if (procRow!=SINGLEPROCROW || procCol!=SINGLEPROCROW) {
		ss << "Process at grid position (" << procRow << "x" << procCol;
		ss << ") do not exist. Only one process at (" << SINGLEPROCROW;
		ss << "x" << SINGLEPROCCOL << ") exists.";
		throw out_of_range(ss.str());
	}
}


/* Blacs general matrix point to point send functions
   use single comm_buffer since with only one process
   only proc 0 can send data to itself

	All Sca/Lapack types supported have a specialization,	so that if user
	tries to instantiate an	unsupported type, compiler issues an error
	However this should not happen since Operator2D contructor prevents
	creation of unsupported types

	Not needed a general template ges2d<class T> since it is already
	declared as a member function
 */
//Lapack send / receive simulation
template <class T>
void send (int M, int N, T *A, int procRow, int procCol) {
	procInGrid (procRow, procCol);
	T *ptr = (T *) comm_buffer;
	ptr = new T[M*N];
	for (int elements=M*N, i=0; i<elements; ++i) ptr[i]=A[i];
}

template <class T>
void receive (int M, int N, T *A, int procRow, int procCol) {
	procInGrid (procRow, procCol);
	T *ptr = (T *) comm_buffer;
	for (int elements=M*N, i=0; i<elements; ++i) A[i]=ptr[i];
	delete [] ptr;
}

//Blacs general matrix point to point send functions
template<>
void Blacs::gesd2d<INTEGER>(int M, int N, INTEGER *A, int LDA, int RDEST, int CDEST) {
	send (M, N, A, RDEST, CDEST);
}
template<>
void Blacs::gesd2d<SINGLE>(int M, int N, SINGLE *A, int LDA, int RDEST, int CDEST) {
	send (M, N, A, RDEST, CDEST);
}
template<>
void Blacs::gesd2d<DOUBLE>(int M, int N, DOUBLE *A, int LDA, int RDEST, int CDEST) {
	send (M, N, A, RDEST, CDEST);
}
template<>
void Blacs::gesd2d<COMPLEX>(int M, int N, COMPLEX *A, int LDA, int RDEST, int CDEST) {
	send (M, N, A, RDEST, CDEST);
}
template<>
void Blacs::gesd2d<DBLCOMPLEX>(int M, int N, DBLCOMPLEX *A, int LDA, int RDEST, int CDEST) {
	send (M, N, A, RDEST, CDEST);
}

//Blacs general matrix point to point receive functions
template<>
void Blacs::gerv2d<INTEGER>(int M, int N, INTEGER *A, int LDA, int RSRC, int CSRC) {
	receive (M, N, A, RSRC, CSRC);
}
template<>
void Blacs::gerv2d<SINGLE>(int M, int N, SINGLE *A, int LDA, int RSRC, int CSRC) {
	receive (M, N, A, RSRC, CSRC);
}
template<>
void Blacs::gerv2d<DOUBLE>(int M, int N, DOUBLE *A, int LDA, int RSRC, int CSRC) {
	receive (M, N, A, RSRC, CSRC);
}
template<>
void Blacs::gerv2d<COMPLEX>(int M, int N, COMPLEX *A, int LDA, int RSRC, int CSRC) {
	receive (M, N, A, RSRC, CSRC);
}
template<>
void Blacs::gerv2d<DBLCOMPLEX>(int M, int N, DBLCOMPLEX *A, int LDA, int RSRC, int CSRC) {
	receive (M, N, A, RSRC, CSRC);
}


//Low-level process signaling (dummy on Lapack since only one process exists)
void Blacs::signalProc(int p) { };
void Blacs::waitProc(int p) { };
void Blacs::signalAll() { };
void Blacs::waitAll() { };
void Blacs::signalAllOthers() { };
void Blacs::waitAllOthers() { };
int Blacs::waitAny() { return NOTINGRID; };
