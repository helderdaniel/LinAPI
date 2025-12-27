/* LinAPI main C++ interface v1.0
   @2009 hdaniel mmadeira
 */

#include "linapi.h"
#include "linexception.h"
#include <string>
#include <cstdlib>


void LinApi::init() {
	thisProc = blacs().processId();
}


LinApi& LinApi::custom(int rows, int cols) {
	static LinApi la = LinApi(rows, cols); 
	return la;
}
LinApi& LinApi::linear() { static LinApi la = LinApi(true); return la; }
LinApi& LinApi::square() { static LinApi la = LinApi(false); return la;  }

LinApi::LinApi (int rows, int cols) { _blacs = &Blacs::custom(rows, cols);	init(); }
LinApi::LinApi (bool linear) {
	if (linear) _blacs = &Blacs::linear();
	else _blacs = &Blacs::square();
	init();
}


/* Althougth Blacs and MPI are Singleton, LinApi is Not since 
   three creation methods with local static LinApi instance exist.
   So it is possible to create one square, one linear and one 
   custom instance. However _blacs pointer point the first _blacs instance
   created, since is singleton.
   
   One way to solve this is to have a static pointer that points the instance
	 after created and is null before creation.
	 custom(), square() and linear() check if pointer null before creating instance
	 If NOT null value of pointer (the sole instance) is returned, 
	 like this (NOT tested in this file yet but tested on demo program 
	 evaluation/singleton class with more than a creation method/singletonmc.cpp):

	#1 Add static class variable in LinApi.h:

		class LinApi {
			static LinApi* pInstance;
	 		(...)
		}

	#2 In LinApi.cpp
	#2.1 Init pointer as NULL		

		LinApi* LinApi::pInstance=0; //or NULL
		
	#2.2 modify methods that create instance to actualize and check pointer

		LinApi& LinApi::custom(int rows, int cols) {
			if (!pInstance) {
				static LinApi la = LinApi(rows, cols); 
				pInstance = &la;
			}
			return *pInstance;
		}
		
		LinApi& LinApi::linear() { 
			if (!pInstance) {			
				static LinApi la = LinApi(true); 
				pInstance = &la;
			}
			return *pInstance;
		}
		
		LinApi& LinApi::square() { 
			if (!pInstance) {			
				static LinApi la = LinApi(false); 
				pInstance = &la;
			}
			return *pInstance;
		}
	
*/



string LinApi::blacsVersion () const {
	return blacs().version();
}


/* Execute LinApi program with Distributed exception handling
  (Stack trace and exception handling could be in a different class)
 */
//stack trace streamming
static void stacktrace(stringstream &ss) {
int thisProc = Blacs::instance().grid().processId();
void * array[25];
int nSize = backtrace(array, 25);
char ** symbols = backtrace_symbols(array, nSize);

	ss << "Process: " << thisProc << " stack trace" << endl;
	for (int i = 0; i < nSize; i++) ss << symbols[i] << endl;
}

//common exception handling
static void common_exp(const char* msg) {
stringstream ss;

	stacktrace(ss);
	ss << msg << endl;
	cerr << ss.str();
	Blacs::instance().signalAll(); //easiest but buggie
	//Blacs::instance().signalProc(Blacs::rootId); //more secure but also buggie
	
	}

//Distributed exception handling
static sem_t endSem;  	//wait for termination semaphore
static int ret=EXIT_SUCCESS;  //return status by default OK

/* Thread that ends this process if exception found in any other process
	Threads can not be member function(it can be static class function or static
	linkage functions like this
 */
void* exceptionGuard (void *arg) {
	Blacs::instance().waitAny(); //easiest but buggie
	//more secure but also buggie
	/*if (Blacs::instance().isRoot()) {
		Blacs::instance().waitAny();
		Blacs::instance().signalAllOthers();
		Blacs::instance().waitAllOthers();
	}
	else {
		Blacs::instance().waitProc(Blacs::rootId);
		Blacs::instance().signalProc(Blacs::rootId);
		exit(0);
	}*/
	ret = EXIT_FAILURE;
	sem_post(&endSem);
	return NULL;
}

//main function for running with exception handling
int LinApi::run(int (*mainptr)(int argc, char **argv, LinApi& la), int &argc, char** &argv) {
stringstream ss;
pthread_t expGuard;

	int s = sem_init (&endSem, 0, 0);
	int t = pthread_create(&expGuard, NULL, exceptionGuard, NULL);
	 //Can not create Guard thread or semaphore, end program
	if (s || t) {
		ss << "Could not start LinApi exception Guard thread or Semaphore" << endl;
		ss << "Press <ctrl>-C to abort" << endl;
		cerr << ss.str();
		exit(EXIT_FAILURE);
	}

	try {
		ret=mainptr(argc, argv, *this);
		sem_post(&endSem);
	}
	catch (const exception &exp) { common_exp(exp.what()); }
	//Not needed since LinException hiearchy descends from exception class
   //catch (const LinException &exp) { common_exp(exp.what()); }
	catch (...) { common_exp("unhandled exception"); }

	//wait endSem
	sem_wait(&endSem);

	return ret;
}



	
