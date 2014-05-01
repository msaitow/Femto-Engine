
CXX = g++
CXXFLAG = -O -std=c++0x -g

LIB_FEMTO = ~/Forge/femto
CXXINCLUDE = -I $(LIB_FEMTO)/include -I /Users/saitow/include
CXXLIB = -L/Users/saitow/lib

.PHONY: default clean

default: worldEngine

SRC = $(notdir $(wildcard *.cc))
OBJ = $(SRC:%.cc=%.o)

worldEngine: $(OBJ)
	$(CXX) -o $@ $(OBJ) -L $(LIB_FEMTO) ${CXXLIB} -lfemto -lboost_filesystem -lboost_system
%.o: %.cc
	$(CXX) $(CXXFLAG) $(CXXINCLUDE) -c $<

clean: 
	-rm -f *.o

install:
	-rm -f ./pack/*
	-mv -f c_* ./pack
	-mv -f f_* ./pack
	tar -zcvf ./pack.tar.gz ./pack
	scp ./pack.tar.gz saitow@qcl.ims.ac.jp:
