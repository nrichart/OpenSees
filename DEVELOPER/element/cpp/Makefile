include ../../Makefile.def

OBJS = Truss2D.o
LIBRARIES = $(OBJS:%.o=%$(SUFFIX))

all:   $(LIBRARIES)
	@$(CD) Macroelement3d/Macroelement3d/; $(MAKE);

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core
	@$(CD) Macroelement3d/Macroelement3d/; $(MAKE) tidy;

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so *.dylib
	@$(CD) Macroelement3d/Macroelement3d/; $(MAKE) clean;

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core
	@$(CD) Macroelement3d/Macroelement3d/; $(MAKE) spotless;

wipe: spotless
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core $(LIBRARIES)
	@$(CD) Macroelement3d/Macroelement3d/; $(MAKE) wipe;


# DO NOT DELETE THIS LINE -- make depend depends on it.


