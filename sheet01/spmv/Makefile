#------------------------------------------------------------------------------
# CSparse Makefile
#------------------------------------------------------------------------------

all: C test

C: library
	( cd Demo ; $(MAKE) )

test: library
	( cd Test; $(MAKE) )

bench: library
	( cd Demo; $(MAKE) )

library:
	( cd Lib ; $(MAKE) )

clean:
	( cd Lib ; $(MAKE) clean )
	( cd Demo ; $(MAKE) clean )
	( cd Test ; $(MAKE) clean )

purge:
	( cd Lib ; $(MAKE) purge )
	( cd Demo ; $(MAKE) purge )
	( cd Test ; $(MAKE) purge )


