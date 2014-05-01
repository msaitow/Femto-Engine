
include ../make.inc

OBJ = $(notdir $(wildcard *.o))

../$(FEMTOLIB): $(OBJ)
	rm -f $@
	$(ARCH) -r $(ARCHFLAGS) $@ $(OBJ)
	$(RANLIB) $@

clean:
	rm -f *.o

%.o: %.cc
	$(CPP) $(CPPFLAG) $(CPPINCLUDE) -c $< -o $@
