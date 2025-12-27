/* BLACS C++ interface v1.0
   @2009 hdaniel, mmadeira
 */

#ifndef _BLACS_H_
#define _BLACS_H_

#include <iostream>
using namespace std;

#include "lintype.h"


const int NOGRIDCONTEXT = -1;
const int NOTINGRID     = -1;


//Process Grid Coordinate class
class GridCoord {
	int _row, _col;
public:
	GridCoord (int r=0, int c=0) { _row=r; _col=c; }
   inline int row() const { return _row; }
   inline int col() const { return _col; }
   friend ostream& operator << (ostream& os, const GridCoord& gc) {
		return os << '(' << gc.row() << 'x' << gc.col() << ')';
	}
};


//ProcessGrid class
class ProcessGrid {
	int thisProc, thisPRow, thisPCol;
	int _gridCols, _gridRows;
   int _context;

   void initGrid (int rows, int cols);

public:
   /* Default constructor, called when variables are declared without parameters,
      for example Blacs class instance variables, only allocate space for variables.
      Initialization of grid must be done by other constructor, not at declaration time, but when grid row and cols or topology are known.
    */
   ProcessGrid ();
	/* create Process grid (rowsxcols)
      Arguments can not have default values, to avoid this being a Default constructor
    */
   ProcessGrid (int rows, int cols);
   ~ProcessGrid ();
   void clear ();
   inline int  context() const { return _context; }
   inline bool processInGrid() const { return !(thisProc==NOTINGRID); }
   inline int  processRow() const { return thisPRow; }
   inline int  processCol() const { return thisPCol; }
   inline int  processId() const  { return thisProc; }
   inline int  gridRows() const   { return _gridRows; }
   inline int  gridCols() const   { return _gridCols; }

   //Return proc processId coordinates in process grid
   GridCoord processRowColId(int processId);

   //Return processID of process at grid coordinates procCoord
	int processAtRowCol(const GridCoord& procCoord);

	//Force all process wait at this point until all of them call this function
	void barrier();

private:
   /* Avoid copy of instance or assignment
      ProcessGrid a(b) or ProcessGrid a=b (copy ctor)
      ProcessGrid a; a=b; (= operator)

		just allows references or pointers
		ProcessGrid& a=b;
		ProcessGrid* a=&b;

		Blacs can not have 2 grids with same context number
		(it must be constructed another one overlaping)
    */
	ProcessGrid(ProcessGrid const& g) {};
	ProcessGrid& operator=(ProcessGrid const& g) { return *this; /*suppress warning*/ }; 
};


//Blacs singleton class
class Blacs { //Singleton
	static string _version;
	int _procCount;

protected:
	ProcessGrid _root, _grid;

public:
	static const int MPIprocs;
	static const int rootId;

   inline ProcessGrid& root() { return _root; }
   inline ProcessGrid& grid() { return _grid; }
   //Return process processId coordinates in underlying process grid
   inline GridCoord processRowColId(int processId){ return _grid.processRowColId(processId);}
   //Return processID of process at underlying process grid coordinates procCoord
	inline int processAtRowCol(const GridCoord& procCoord){ return _grid.processAtRowCol(procCoord);}
   inline string version() const { return _version; }
   inline bool isRoot() const { return _grid.processId()==rootId; }
	inline int processId() const { return _grid.processId(); }
	inline int procCount () { return _procCount; }

	//Force all process on underlying grid wait at this point until
   //all of them call this function
	void barrier() { grid().barrier(); }

	//General matrix point to point send functions
	template <class T> void gesd2d(int M, int N, T *A, int LDA, int RDEST, int CDEST);
	template <class T> void gerv2d(int M, int N, T *A, int LDA, int RSRC, int CSRC);

	//Low level process signaling
	void signalProc(int p);
	void waitProc(int p);
	void signalAll();
	void waitAll();
   void signalAllOthers();
	void waitAllOthers();
   int waitAny();


   /* Create/Get sole (static) instance
    */
   //Create custom grid 
	static Blacs& custom(int rows, int cols);
	//Create Linear grid with all MPI processes.
	static Blacs& linear();
	//Create grid as close to square as possible some MPI processes may be left out of grid
	static Blacs& square();
	//get sole instance. If not created yet, create square instance and return it (on Lapack is (1,1))
   static Blacs& instance() { return square(); };


private: //or protected if derivation allowed. Hide ctors and dctor
	Blacs (int rows, int cols);
   ~Blacs (); 
  	//To avoid instantiation like: Blacs BLACS=Blacs::instance();
  	//Instanciate it and access it like: Blacs::instance().procId()
	Blacs(Blacs const&) {};    // copy ctor
  	Blacs& operator=(Blacs const& m) { return *this; /*suppress warning*/ };  // assign op

};


#endif