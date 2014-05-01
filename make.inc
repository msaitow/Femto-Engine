
SHELL = /bin/sh

# Home directory for FEMTO
FEMTO_HOME=${HOME}/Forge/Femto

# Include sirectories
CPPINCLUDE = -I$(FEMTO_HOME)/include -L${HOME}/Forge/lib

# C++ compiler specification
CPP      =  /usr/local/bin/g++
CPPFLAG  = -O -std=c++0x -g -I/Users/saitow/Forge/include
#CPPFLAG = -std=c++11 -stdlib=libc++ -I/Users/saitow/include
#CPPFLAG = -g

#*NOT YET* CPP     = clang++
#*NOT YET* CPPFLAG = -O3 -g -I${HOME}/include  

# The archiver specification
ARCH      = ar
APCHFLAGS = cr
RANLIB    = ranlib

# Name of the resultant FEMTO library
FEMTOLIB = libfemto.a
