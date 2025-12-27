/* Check BLACS C++ interface v1.0
   @2009 hdaniel, mmadeira
 */

#ifndef _BLACSCHK_H_
#define _BLACSCHK_H_

#include <iostream>
using namespace std;

#include "blacs.h"


/* BlacsCheck class
 */
class BlacsCheck {

public:
        //check basic procs comm
	static void comm (ostream& s);

//Avoid instanciation
private: //Hide ctors, dctor and = operator
	BlacsCheck () { }
	~BlacsCheck  (); 
  	BlacsCheck (BlacsCheck  const&) {};    // copy ctor
  	BlacsCheck & operator=(BlacsCheck const& bc) { return *this; /*suppress warning*/ };
};

#endif