# Compiler
CC=gcc
CXX=g++
AR=ar
LD=ld

DYN_SUFFIX=.so
DYN_OPT=-shared -Wl,-soname,$(shell basename $(DYN_PRODUCT))

VERSION=1.10.1
PREFIX=/usr/local

PATH_nuSQUIDS=/home/dcianci/Physics/GlobalFits/IcecubePack/nuSQuIDS
PATH_SQUIDS=/home/dcianci/Physics/GlobalFits/IcecubePack/SQuIDS

SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(patsubst src/%.cpp,build/%.o,$(SOURCES))

CXXFLAGS= -std=c++17

# Directories

GSL_CFLAGS=-I/usr/local/include
GSL_LDFLAGS=-L/usr/local/lib -lgsl -lgslcblas -lm
HDF5_CFLAGS=-I/usr/include/hdf5/serial
HDF5_LDFLAGS=-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_hl -lhdf5 -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lsz -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial
SQUIDS_CFLAGS=-Wno-abi -I/usr/local/include
SQUIDS_LDFLAGS=-L/usr/local/lib -lSQuIDS -lgsl -lgslcblas -lm

INCnuSQUIDS=$(PATH_nuSQUIDS)/include
LIBnuSQUIDS=$(PATH_nuSQUIDS)/lib

# FLAGS
CFLAGS= -O3 -fPIC -I$(INCnuSQUIDS) $(SQUIDS_CFLAGS) $(GSL_CFLAGS) $(HDF5_CFLAGS)
LDFLAGS= -Wl,-rpath -Wl,$(LIBnuSQUIDS) -L$(LIBnuSQUIDS)
LDFLAGS+= $(SQUIDS_LDFLAGS) $(GSL_LDFLAGS) $(HDF5_LDFLAGS)
EXMAPLES_FLAGS=-I$(INCnuSQUIDS) -I$(PATH_SQUIDS)/include $(CXXFLAGS) $(CFLAGS)

# Project files
NAME=nuSQuIDS
STAT_PRODUCT=lib/lib$(NAME).a
DYN_PRODUCT=lib/lib$(NAME)$(DYN_SUFFIX)

d1:	src/icecube.cxx
	@$(CXX) $(EXMAPLES_FLAGS) src/icecube.cxx -lnuSQuIDS $(LDFLAGS) -o icecube
