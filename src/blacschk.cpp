/* Check BLACS C++ interface v1.0
   Scalapack version
   @2009 hdaniel mmadeira
 */

#include "blacschk.h"
#include <iomanip>
#include <sstream>
#include <cstdlib>


/* BlacsCheck class implementation
 */


void BlacsCheck::comm (ostream& os) {
/* Blacs singleton reference can not be static variable because if so it is
   executed before program, defining default square grid, avoiding user to
   define its own topology
*/
Blacs &blacs=Blacs::instance();
stringstream outs;
ProcessGrid& g = blacs.grid();
GridCoord gc;

//Get grid Rows and Colums and this process grid coordinates
int thisPRow=g.processRow(), thisPCol=g.processCol(),
    npRow=g.gridRows(), npCol=g.gridCols();

INTEGER thisProc=g.processId();

INTEGER sendProcRow, sendProcCol, senderProc;

	//If thisProc is located at grid (0,0), receive check-in messages from all other nodes
	//if ((thisPRow == 0) && (thisPCol == 0)) {
	if (blacs.isRoot()) {
		outs << endl << blacs.version() << endl;
		outs << "Proc=" << blacs.grid().processId() << endl;
		outs << "Grid topology = " << thisPRow << "x" << thisPCol << endl;
		for (int i=0; i <= npRow-1; ++i) {
			for (int j=0; j <= npCol-1; ++j) {
				//Proc (0,0) receive check-in messages from all process but itself
				if ((i != 0) || (j != 0)) {
					//Receive message from process
					blacs.gerv2d(1, 1, &senderProc, 1, i, j);

					//Check sender process (procGridId) is located on right process grid coordinates
					gc=g.processRowColId(senderProc);
					sendProcRow=gc.row();
					sendProcCol=gc.col();
					if ((sendProcRow != i) || (sendProcCol != j)) {
						outs << "Process " << senderProc << " grid position should be (";
						outs << i << "x" << j << ") but is reported to be (";
						outs <<  sendProcRow << "," << sendProcCol << ")" << endl;
						outs << "Grid error!  Halting . . ." << endl;
						os << outs.str();
						exit (EXIT_FAILURE);
					}
					outs << "SENDER process   {" << setw(2) << i << "," << setw(2) << j;
					outs << "} (node no.=" << setw(4) << senderProc;
					outs << ") has checked in at RECEIVER process = " << thisProc;
					outs << "." << endl;
				} // If not proc (0,0) 
				else {
					outs << "RECEIVER process {" << setw(2) << i << ",";
					outs << setw(2) << j << "} (node no.=" << setw(4) << thisProc;
					outs << ") STARTED!" << endl;
				}
			} // for j
		} // for i
		outs << endl << "All processes checked in.  Run finished." << endl;
	} // outer if

	//if thisProc != proc (0,0) send check-in message with its ID to proc (0,0)
	else {
		//This code runs in all processors but (0,0) so concurrent writes need
		//be mounted on a string before send it to cout to avoid mixing
		outs << "Process " << thisProc << " sent check in message" << endl;
		blacs.gesd2d(1, 1, &thisProc, 1, 0, 0);
	}
	os << outs.str();
        blacs.barrier(); //wait to not mix processes output        
}
