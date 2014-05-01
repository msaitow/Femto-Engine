
CXX = g++
CXXFLAG = -O3 -g

LIB_FEMTO = ~/Dropbox/femto
CXXINCLUDE = -I $(LIB_FEMTO)/include
C_TARGET = ~/orz/src/sci/ctnew2
F_TARGET = ${C_TARGET}/fsrc/mrci_tensors

.PHONY: default clean

default: eclipse

SRC = $(notdir $(wildcard *.cc))
OBJ = $(SRC:%.cc=%.o)

eclipse: $(OBJ)
	$(CXX) -o $@ $(OBJ) -L $(LIB_FEMTO) -lfemto
%.o: %.cc
	$(CXX) $(CXXFLAG) $(CXXINCLUDE) -c $<

clean: 
	-rm -f *.o
	-rm -f c_*
	-rm -f f_*

install:
	mv c_* ${C_TARGET}
	mv f_* ${F_TARGET}
