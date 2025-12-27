/* MPIch2 C++ interface
   @2009 hdaniel mmadeira
 */

#include <mpi.h>
#include "mpiobj.h"
#include <cmath>
#include <cctype>
#include <cstdlib>

/* External BLACS C interface functions prototypes
   If not defined it gives warnings with mpicc -Wall (gcc ANSI C)
   If not defined it gives errors with mpicxx (g++ C++)   
   Anyway prototypes are advisable to check variables type on compile
*/
extern "C" {
	int MPI_Initialized(int *flag);
	int MPI_Init(int *argc, char ***argv);
	int MPI_Comm_size(MPI_Comm comm, int *size);
	int MPI_Comm_rank(MPI_Comm comm, int *rank);

	int MPI_Send(void* buf, int count, MPI_Datatype datatype, int dest,
	int tag, MPI_Comm comm);	
	int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source,
	int tag, MPI_Comm comm, MPI_Status *status);
}

/* Mpi singleton class implementation
 */

/* Determine this process number and the number of processes
   passed on "mpiexec -np <no. processes>" and this Process number (zero based)

   use MPI primitive function called on Cblacs_pinfo(int*, int*);
   Functions can be called only after MPI initialized with MPI_init(...)
	which must be called only once. Cblacs_pinfo calls MPI_init only if not
	called already.
   This check is done with MPI_Initialized(int *flag), where if flag != 0
	MPI_init already called
   
   This constructor body is equivalent to: Cblacs_pinfo(&_procId, &_procCount);
   however there is no point in call a Blacs function since MPI does not
   need Blacs
 */
Mpi::Mpi() {
int flag;

	MPI_Initialized(&flag);
	if (!flag) MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &_procCount);
	MPI_Comm_rank(MPI_COMM_WORLD, &_procId);
}


void Mpi::signalProc(int p) {
int token=0;
	MPI_Send(&token, 1, MPI_INT, p, 99, MPI_COMM_WORLD);
}

void Mpi::signalAll(int l) {
int token=0;
const int limit = (l<0) ? _procCount : l;
	for (int dest=0; dest<limit; ++dest)
		MPI_Send(&token, 1, MPI_INT, dest, 99, MPI_COMM_WORLD);
}

void Mpi::signalAllOthers(int l) {
int token=0;
const int limit = (l<0) ? _procCount : l;
	for (int dest=0; dest<limit; ++dest)
		if (dest!=_procId)
			MPI_Send(&token, 1, MPI_INT, dest, 99, MPI_COMM_WORLD);
}

void Mpi::waitProc(int p) {
int token;
	MPI_Recv(&token, 1, MPI_INT, p, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void Mpi::waitAll(int l) {
//int token;
const int limit = (l<0) ? _procCount : l;
	for (int src=0; src<limit; ++src)
		//MPI_Recv(&token, 1, MPI_INT, src, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		waitAny(); //cause it is imprevisible which one send first
		//however if same proc send X times it is accounted as if X procs have send signal
		//A more secure function should check if all processors sends signal
		//using for instance an array for marking
		//if one processor sends more than one signal it should print an error and abort
}

void Mpi::waitAllOthers(int l) {
//int token;
const int limit = (l<0) ? _procCount : l;
	for (int src=0; src<limit; ++src)
		if (src!=_procId)
			//MPI_Recv(&token, 1, MPI_INT, src, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			waitAny();  //(read comment on waitAll)
}

int Mpi::waitAny() {
int token;
MPI_Status status;
	MPI_Recv(&token, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	return status.MPI_SOURCE;
}