
CXX = g++
CXXFLAG = -O3 -std=c++0x -g
#CXXFLAG = -g

#*NOT YET* CXX = clang++
#*NOT YET* CXXFLAGS = -O3 -g -I${HOME}/include 

LIB_FEMTO  = ~/Forge/femto
CXXINCLUDE = -I $(LIB_FEMTO)/include -I/Users/saitow/include
CXXLIB     = -L /Users/saitow/lib
C_TARGET   = ~/orz/src/sci/ctnew2
F_TARGET   = ${C_TARGET}/fsrc/mrci_tensors

.PHONY: default clean

default: eclipse

SRC = $(notdir $(wildcard *.cc))
OBJ = $(SRC:%.cc=%.o)

eclipse: $(OBJ)
	$(CXX) -o $@ $(OBJ) $(CXXLIB) -L $(LIB_FEMTO) -lfemto -lboost_serialization -lboost_wserialization

%.o: %.cc
	$(CXX) $(CXXFLAG) $(CXXINCLUDE) -c $<

clean: 
	-rm -f *.o
	-rm -f c_*
	-rm -f f_*

install_c:
	-mv c_* ${C_TARGET}

install_f:
	-mv f_* ${C_TARGET}

install:
	-mv c_* ${C_TARGET}
	-mv f_* ${F_TARGET}
