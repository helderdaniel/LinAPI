#Makefile to build Example for static and shared LinAPI libraries
#Lapack and Scalapack versions

all: lapst lapsh scast scash

lapst:
	$(MAKE) -f MakeLapack static
lapsh:
	$(MAKE) -f MakeLapack shared
scast:
	$(MAKE) -f MakeScalapack static
scash:
	$(MAKE) -f MakeScalapack shared

clean: 
	#-(cmd) supress error reporting from (cmd) back to make, so that
	#make continues without needed for -i (ignore error report) switch:
	# supress error if file do not exist (ex: -rm *.~)
	# supress error if file: link do already exist (ex: -ln somefile link)
	$(MAKE) -f MakeLapack clean
	$(MAKE) -f MakeScalapack clean
	
cleanall: 
	$(MAKE) -f MakeLapack cleanall
	$(MAKE) -f MakeScalapack cleanall
	
