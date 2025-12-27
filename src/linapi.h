/* LinAPI main C++ interface v1.0
   @2009 hdaniel mmadeira
 */

#ifndef _LINAPI_H_
#define _LINAPI_H_

/*	Include all modules so that client program only include linapi.h
	(and template definitions, see comment below)
*/
#include "blacs.h"
#include "blacschk.h"
/* Can not include in library file template classes which
   will be included in client file.
	Al least while gcc do not supports export keyword
	if so linker finds multiple references.

	So following template classes must be included in clien (if needed)
	along with linapi.h

	#include "operator2d.h"
	#include "linmatrix.h"
	#include "linvector.h"
*/

#include <iostream>
#include <sstream>
#include <pthread.h>
#include <semaphore.h>
using namespace std;


//LinApi singleton class
class LinApi {
	int thisProc;
	Blacs* _blacs;

public:
	static LinApi& custom(int rows, int cols);
   static LinApi& linear();
	static LinApi& square();
	static LinApi& instance() { return square(); };

	Blacs& blacs() const { return *_blacs; };
	string blacsVersion () const;
	int run(int (*mainptr) (int argc, char** argv, LinApi&), int &argc, char** &argv);

private: //Hide ctors and dctor
	LinApi (int gridRows, int gridcols);
	LinApi (bool linear);
   ~LinApi () {};  
	LinApi(LinApi const&) {};    // copy ctor
  	LinApi& operator=(LinApi const&) { return *this; /*suppress warning*/ };
	//Note: Also ctors, and = op can be declared without body
	void init();
};


#endif
