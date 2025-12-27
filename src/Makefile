#Makefile to build static and shared LinAPI libraries
#Lapack and Scalapack versions

all:
	#-C change to dir
	$(MAKE) -C lapack all
	$(MAKE) -C scalapack all

clean: 
	#-(cmd) supress error reporting from (cmd) back to make, so that
	#make continues without needed for -i (ignore error report) switch:
	# supress error if file do not exist (ex: -rm *.~)
	# supress error if file: link do already exist (ex: -ln somefile link)
	#$(MAKE) -C lapack clean
	$(MAKE) -C scalapack clean	
	@-rm *~	
	
cleanall: 
	$(MAKE) -C lapack cleanall
	$(MAKE) -C scalapack cleanall
	@-rm *~
	
