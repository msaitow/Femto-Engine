
include make.inc

all: femtolib

clean: femtoclean

femtolib:
	( cd core;     $(MAKE) )
#	( cd orz;      $(MAKE) )
	( cd reaktor;  $(MAKE) )
	( cd reaktor2; $(MAKE) )
	( cd portal;   $(MAKE) )
	( cd obj;      $(MAKE) )

femtoclean:
	( cd core;     $(MAKE) clean )
#	( cd orz;      $(MAKE) clean )
	( cd reaktor;  $(MAKE) clean )
	( cd reaktor2; $(MAKE) clean )
	( cd portal;   $(MAKE) clean )
	( cd obj;      $(MAKE) clean )

