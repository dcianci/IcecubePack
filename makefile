# Compiler
CC=gcc
CXX=g++
AR=ar
LD=ld

DYN_SUFFIX=.so
DYN_OPT=-shared -Wl,-soname,$(shell basename $(DYN_PRODUCT))

VERSION=1.10.1
PREFIX=/usr/local

PATH_nuSQUIDS=/a/share/amsterdam/dcianci/nuSQuIDS
PATH_SQUIDS=/a/share/amsterdam/dcianci/SQuIDS
PATH_GSL=/a/share/amsterdam/dcianci/gsl

SOURCES = $(wildcard src/*.cpp)
OBJECTS = $(patsubst src/%.cpp,build/%.o,$(SOURCES))

CXXFLAGS=-std=c++17

# Directories
GSL_CFLAGS=-I$(PATH_GSL)/include
GSL_LDFLAGS=-Wl,-rpath -Wl,$(PATH_GSL)/lib -L$(PATH_GSL)/lib -lgsl -lgslcblas -lm
HDF5_CFLAGS=-I/usr/include
#HDF5_CFLAGS=-I/usr/include/hdf5/serial
HDF5_LDFLAGS=-L/usr/lib64/ -lhdf5_hl -lhdf5# -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib64/
#HDF5_LDFLAGS=-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_hl -lhdf5 -Wl,-Bsymbolic-functions -Wl,-z,relro -lpthread -lz -ldl -lm -Wl,-rpath -Wl,/usr/lib/x86_64-linux-gnu/hdf5/serial
SQUIDS_CFLAGS=-I$(PATH_SQUIDS)/include
SQUIDS_LDFLAGS=-Wl,-rpath -Wl,$(PATH_SQUIDS)/lib -L$(PATH_SQUIDS)/lib -lSQuIDS
NUSQUIDS_CFLAGS=-I$(PATH_nuSQUIDS)/include
NUSQUIDS_LDFLAGS=-Wl,-rpath -Wl,$(PATH_nuSQUIDS)/lib -L$(PATH_nuSQUIDS)/lib -lnuSQuIDS
ROOT_CFLAGS=$(shell root-config --cflags)
ROOT_LDFLAGS=$(shell root-config --cflags  --glibs) -lMinuit

# FLAGS
CFLAGS= $(CXXFLAGS) -O3 -fPIC $(SQUIDS_CFLAGS) $(NUSQUIDS_CFLAGS) $(GSL_CFLAGS) $(HDF5_CFLAGS) $(ROOT_CFLAGS)
LDFLAGS= $(SQUIDS_LDFLAGS) $(NUSQUIDS_LDFLAGS) $(GSL_LDFLAGS) $(HDF5_LDFLAGS) $(ROOT_LDFLAGS)
$(info $(CFLAGS))
$(info $(LDFLAGS))


# Project files
NAME=nuSQuIDS
STAT_PRODUCT=lib/lib$(NAME).a
DYN_PRODUCT=lib/lib$(NAME)$(DYN_SUFFIX)

d1:	
	@$(CXX) $(CFLAGS) src/icecubeScan_v0.cxx $(LDFLAGS) -o build/icecubescan_v0
	@$(CXX) $(CFLAGS) src/icecubeScan_v1.cxx $(LDFLAGS) -o build/icecubescan_v1
