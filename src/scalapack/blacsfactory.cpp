/* BlacsFactory v0.3
   Not needed on current implementation
   Test version available only for Scalapack (not for LAPACK)
   2009 hdaniel
*/

#include "blacsfactory.h"
#include "mpiobj.h"
#include <cstdlib>

/* BlacsFactory class implementation (static factory)
 */
int BlacsFactory::MPIprocs = Mpi::instance().procCount();

Blacs& BlacsFactory::custom(int rows, int cols) {
  	Blacs& b=Blacs::custom(rows, cols);
	return b;
}
Blacs& BlacsFactory::linear() { return custom(MPIprocs, 1); }
Blacs& BlacsFactory::square() {
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
