#------------------------------------------------------------------------------
# CSparse Makefile
#------------------------------------------------------------------------------

C:
	( cd Lib ; $(MAKE) )
	( cd Demo ; $(MAKE) )
	( cd Test ; $(MAKE) )

all: C 

library:
	( cd Lib ; $(MAKE) )

test:
	( cd Test ; $(MAKE) )

clean:
	( cd Lib ; $(MAKE) clean )
	( cd Demo ; $(MAKE) clean )

purge:
	( cd Lib ; $(MAKE) purge )
	( cd Demo ; $(MAKE) purge )


