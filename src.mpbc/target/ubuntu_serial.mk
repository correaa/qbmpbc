#knonw to work with ubuntu 10.10 and qbox-1.52.3 use as 'export TARGET=ubuntu', 'make ubuntu', 'make qb lib'
all: lib
	mv libqb.a libqb_serial.a
ubuntu: 
	sudo apt-get install \
          fftw-dev libxerces-c-dev libscalapack-mpi-dev g++ \
          libatlas-dev liblapack-dev libblas-dev libblacs-mpi-dev \
          uuid-dev libopenmpi-dev 
#last 2 lines above are redundant

PLT=$(shell uname -m)
USRDIR=

MPIDIR=$(USRDIR)

XERCESCQ=USE_XERCES
XERCESCDIR=$(USRDIR)
XERCESCFLAG=XERCESC_3

FFTWDIR=$(USRDIR)
BLASDIR=$(USRDIR)
LAPACKDIR=$(USRDIR)

FFTW3DIR=$(USRDIR)

CXXDIR=
CXX=$(CXXDIR)c++
LD=$(CXX) 
PLTFLAGS += \
    -DUSE_FFTW \
    -D$(XERCESCQ) -D$(XERCESCFLAG) \
    -DADD_ -DAPP_NO_THREADS -DXML_USE_NO_THREADS
INCLUDE =    

LIBPATH = 

CXXFLAGS = -O4 -D$(PLT) \
           $(INCLUDE) \
           $(PLTFLAGS) \
           $(DFLAGS) \
           -Wno-write-strings -Wno-unused-result -Wfatal-errors

LIBS =  $(PLIBS) \
        -lfftw \
        -lxerces-c \
        -llapack \
        -luuid \
	-boost_system

LDFLAGS = $(LIBPATH) $(LIBS) -Wl,-rpath=$(HOME)/usr/lib 

