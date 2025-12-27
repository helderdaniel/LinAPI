/* Blacs class for LinAPI C++ demo
   @2009  hdaniel, mmadeira
*/

#define _EXCEPTION_HANDLER_

//LinAPI headers includes all required such as blacs.h
#include <linapi.h>

#include <cstdlib>

/* On Lapack must show only:

argc=1, argv=./blacstestLapackst
Blacs v1.0 lapack
Proc=0
Grid topology = 0x0
RECEIVER process { 0, 0} (node no.=   0) STARTED!
All processes checked in.  Run finished.
*/

#ifdef _EXCEPTION_HANDLER_
int mainThread (int argc, char ** argv, LinApi& la) {

	Blacs& blacs=la.blacs();
	if (blacs.isRoot()) printf("argc=%d, argv=%s\n", argc, argv[0]);
	blacs.barrier(); //wait for all proceses to output info to avoid mixing with next code segment

	//Check Blacs process grid
	BlacsCheck::comm(cout);
	blacs.barrier();

	return EXIT_SUCCESS;
}

int main (int argc, char ** argv) {
	LinApi& la=LinApi::square();
	return la.run(mainThread, argc, argv);
}  
#else

int main (int argc, char ** argv) {

	// Singleton LinApi instanciates singleton Blacs
	LinApi& la1=LinApi::square();
	LinApi& la2=LinApi::linear();
	LinApi& la3=LinApi::custom(1,2);
	
	//Like LinApi, Blacs is singleton only first object is created
	//Next attempts only return reference for already created object
	Blacs& b1=Blacs::square();
	Blacs& b2=Blacs::linear();
	Blacs& b3=Blacs::custom(1,2);
	Blacs& b4=la3.blacs();
	if (b4.isRoot()) printf("argc=%d, argv=%s\n", argc, argv[0]);
	b4.barrier(); //wait for all proceses to output info to avoid mixing with next code segment

	//Check Blacs process grid
	BlacsCheck::comm(cout);
	b4.barrier();
   
	return EXIT_SUCCESS;
}  

#endif