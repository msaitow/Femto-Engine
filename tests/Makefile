
CXX = g++
CXXFLAG = -O -g

LIB_FEMTO = ~/Drops/femto
CXXINCLUDE = -I $(LIB_FEMTO)/include

.PHONY: default clean

default: femto_test

SRC = $(notdir $(wildcard *.cc))
OBJ = $(SRC:%.cc=%.o)

femto_test: $(OBJ)
	$(CXX) -o $@ $(OBJ) -L $(LIB_FEMTO) -lfemto
%.o: %.cc
	$(CXX) $(CXXFLAG) $(CXXINCLUDE) -c $<

clean: 
	-rm -f *.o
