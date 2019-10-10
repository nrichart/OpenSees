include ../../Makefile.def

OBJS = ElasticPPcpp.o BeamFricitonSupport.o CompressionDamage1d.o OrthotropicMembraneSection.o TensionDamage1d.o
LIBRARIES = $(OBJS:%.o=%$(SUFFIX))

all:   $(LIBRARIES)

BeamFricitonSupport.o: BeamFrictionSupport/BeamFrictionSupport/BeamFrictionSupport.cpp
	$(CC++) -g -Wall $(INCLUDES) -fPIC -c $< -o $@

CompressionDamage1d.o: CompressionDamage1d/CompressionDamage1d/CompressionDamage1d.cpp
	$(CC++) -g -Wall $(INCLUDES) -fPIC -c $< -o $@

OrthotropicMembraneSection.o: OrthotropicMembraneSection/OrthotropicMembraneSection.cpp
	$(CC++) -g -Wall $(INCLUDES) -fPIC -c $< -o $@

TensionDamage1d.o: TensionDamage1d/TensionDamage1d/TensionDamage1d.cpp
	$(CC++) -g -Wall $(INCLUDES) -fPIC -c $< -o $@

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so *.dylib

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core $(LIBRARY)

# DO NOT DELETE THIS LINE -- make depend depends on it.
