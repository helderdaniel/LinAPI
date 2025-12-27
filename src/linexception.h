/* LinAPI exception handling v1.0
   @2009 hdaniel mmadeira
 */

#ifndef _LINEXCEPTION_H_
#define _LINEXCEPTION_H_

#include <execinfo.h>
#include <exception>
#include <stdexcept>
#include <typeinfo>
#include <sstream>
using namespace std;


/* Default LinApi exceptions behaviour
	Do not descent from runtime_error since member strin _M_msg is private
 */
class LinException : public exception {
protected:
	string _msg;

public:
	explicit LinException (const string &s="undefined LinException", int c=0) {
		stringstream ss;
			ss << s << "code: " << c << endl;
			_msg=ss.str();
	}
	~LinException() throw() {}
	const char* what() const throw() { return _msg.c_str(); }
};

//Arithmetic type not supported
template <class T, class C>
class LinTypeNotSupportedException : public LinException {

public:
	explicit	LinTypeNotSupportedException (const string &method="") {
		stringstream ss;
			ss << "LinApi not supported type <";
			ss << typeid(T).name() << ">";
			ss << " in " << typeid(C).name() << "::" << method << endl;
			_msg=ss.str();
	}
};

#endif
