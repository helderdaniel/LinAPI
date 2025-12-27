/* MPIch2 C++ interface v1.0
   @2009 hdaniel, mmadeira
 */

#ifndef _MPIOBJ_H_
#define _MPIOBJ_H_

#include <cstdio>


//MPI singleton class
class Mpi {
	int _procId, _procCount;

private: //or protected if derivation allowed. Hide ctors and dctor
  	Mpi();
	~Mpi() {};
  	//To avoid instantiation like: Mpi MPI=Mpi::instance();
  	//Instanciate it and access it like: Mpi::instance().procId()
  	Mpi(Mpi const&) {};    // copy ctor
  	Mpi& operator=(Mpi const& m) { return *this; /*suppress warning*/ };  // assign op

public:
	static const int rootId = 0;

	//Get sole (static) instance
	static Mpi& instance() { static Mpi s; return s; }
	int procId () { return _procId; }
	int procCount () { return _procCount; }

	void signalProc(int p);
	void waitProc(int p);
	int waitAny(); //return sender proc id

	//l can be used as an upper limit to processe range
	//if l < 0 all MPI procs are used to broadcast
	void signalAll(int l=-1);
	void waitAll(int l=-1);
	void signalAllOthers(int l=-1);
	void waitAllOthers(int l=-1);
	
};

#endif