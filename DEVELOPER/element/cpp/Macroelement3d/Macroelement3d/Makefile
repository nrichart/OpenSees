include ../../../../Makefile.def

INCLUDES = -I../../../../core
OBJS = Macroelement3d.o CohesiveSurface.o DamageShearInterface.o GambarottaLagomarsinoModel.o Macroelement3d.o NoTensionSection3d.o WrappedMaterial.o
LIBRARIES = Macroelement3d$(SUFFIX)

all:  $(LIBRARIES)

Macroelement3d.so: $(OBJS)
	$(CC++) -shared -Wl,-soname,$@ -o $@ $^ -lc
	$(CP) $@ $(prefix)/lib

Macroelement3d.dylib: $(OBJS)
	$(CC++) $(OS_FLAG) -dynamiclib -Wl,-undefined,suppress,-flat_namespace $^  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $@
	$(CP) $@ $(prefix)/lib

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so *.dylib

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core $(LIBRARY)
