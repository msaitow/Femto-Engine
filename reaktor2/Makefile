
include ../make.inc

SRC = $(notdir $(wildcard *.cc))
OBJ = $(SRC:%.cc=%.o)

# ../$(FEMTOLIB): $(OBJ)
# 	rm -f $@
# 	$(ARCH) -r $(ARCHFLAGS) $@ $(OBJ)
# 	$(RANLIB) $@

../$(FEMTOLIB): $(OBJ)
	cp $(OBJ) ../obj

clean:
	rm -f *.o

%.o: %.cc
	$(CPP) $(CPPFLAG) $(CPPINCLUDE) -c $< -o $@
