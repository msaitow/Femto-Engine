
CXX = g++
CXXFLAG = -g

LIB_FEMTO = ~/Drops/femto
CXXINCLUDE = -I $(LIB_FEMTO)/include

.PHONY: default clean

default: mrci_sigma

SRC = $(notdir $(wildcard *.cc))
OBJ = $(SRC:%.cc=%.o)

mrci_sigma: $(OBJ)
	$(CXX) -o $@ $(OBJ) -L $(LIB_FEMTO) -lfemto
%.o: %.cc
	$(CXX) $(CXXFLAG) $(CXXINCLUDE) -c $<

clean: 
	-rm -f *.o
