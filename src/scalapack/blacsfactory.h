/* BlacsFactory v0.3
   Not needed on current implementation
   Test version available only for Scalapack (not for LAPACK)
   2009 hdaniel
*/

#ifndef _LINAPI_H_
#define _LINAPI_H_

#include "../blacs.h"

/* Blacs Factory class
   that need only Blacs::custom(int r, int c) get/creation instance function

   Blacs::linear()
	Blacs::square()
	Blacs::instance()  are all defined here calling Blacs::custom(int r, int c)
 */

//Blacs factory class (static factory pattern)
class BlacsFactory {
static int MPIprocs;

public:
	static Blacs& custom(int rows, int cols);
	// Create Linear grid with all MPI processes.
	static Blacs& linear();
	/* Create grid as close to square as possible
      some MPI processes may be left out of grid
    */
	static Blacs& square();
   static Blacs& instance() { return square(); };

//Avoid instanciation of full static class
private: //Hide ctors, dctor and = operator
	BlacsFactory () { }
   ~BlacsFactory ();
  	BlacsFactory (BlacsFactory const&) {};    // copy ctor
  	BlacsFactory& operator=(BlacsFactory const& bF) { return *this; /*suppress warning*/ };
};

#endif