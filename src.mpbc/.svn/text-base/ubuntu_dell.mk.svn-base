#-------------------------------------------------------------------------------
#  ubuntu_hp.mk
#  makefile-include for an Ubuntu distribution,
#  do 'make ubuntu' to install all the system packages
#  tested on Ubuntu 9.10
#-------------------------------------------------------------------------------
 PLT=Linux_x8664
#-------------------------------------------------------------------------------

USRDIR=/usr
MPIIMPL=
#.openmpi
#.mpich
MPIDIR=$(USRDIR)/lib/$(MPIIMPL)
XERCESCDIR=$(USRDIR)
FFTWDIR=$(USRDIR)
BLASDIR=$(HOME)/usr
ATLASDIR=$(USR)
LAPACKDIR=$(USR)
HDF5DIR=$(HOME)/usr
BOOSTDIR=$(HOME)/usr

FFTW3DIR=$(USRDIR)
MPBCDIR=./include

CXXDIR=$(HOME)/usr/bin/
CXX=$(CXXDIR)mpicxx$(MPIIMPL)
LD=$(CXX) 

PLTFLAGS += -DIA32 -DUSE_FFTW -DUSE_XERCES -DUSE_MPI -DADD_ -DSCALAPACK -DAPP_NO_THREADS -DXML_USE_NO_THREADS

INCLUDE_RAW = \
    -I$(FFTWDIR)/include \
    -I$(XERCESCDIR)/include \
    -I$(FFTW3DIR)/include \
    -I$(HDF5DIR)/include \
    -I$(MPIDIR)/include \
    -I$(MPBCDIR) \
    -I$(BOOSTDIR)/include
    
INCLUDE = $(INCLUDE_RAW)
    

CXXFLAGS = -O4 -D$(PLT) \
           $(INCLUDE) \
           $(PLTFLAGS) \
           $(DFLAGS) \
           -Wno-write-strings -Wfatal-errors

LIBPATH = -L$(FFTWDIR)/lib \
          -L$(MPIDIR)/lib \
          -L$(LAPACKDIR)/lib \
          -L$(BLASDIR)/lib \
          -L$(ATLASDIR)/lib \
          -L$(XERCESCDIR)/lib \
          -L$(HDF5DIR)/lib

#for system-wide boost libraries, leave empty, don't mix versions of boost
BOOST_VER = 
#-gcc44-mt


LIBS =  $(PLIBS) \
        -lfftw \
        -lfftw3 \
        -lxerces-c \
        -lboost_filesystem$(BOOST_VER) -lboost_system$(BOOST_VER) \
        -llapack -lptf77blas -latlas \
        -lgfortran
#         -lmpich
#         -Xlinker -Bstatic \
#         -lgfortran -static-libgcc -lmpich  -lrt -pthread \
#         -Xlinker -Bdynamic -lsicuuc -lsicudata
#         -lm -lc

LDFLAGS = $(LIBPATH) $(LIBS) -Wl,-rpath=$(HOME)/usr/lib 

PLAT=Linux_x8664

# Blacs libraries
# BLACSDIR = $(USR)
# for some reason the BLACS included in ubuntu 9.10 returns an incorrect number of processors
# so, in the following I use a manually installed version

#user compiled version
 BLACSDIR   = $(HOME)/usr/lib
#distribution compiled version
#BLACSDIR   = $(USRDIR)/lib

#user compiled version only
 BLACSPOSTFIX  = -mpi
#or -mpich distribution compiled version only
#BLACSPOSTFIX = -openmpi

BLACSFINIT  = $(BLACSDIR)/libblacsF77init$(BLACSPOSTFIX).a
BLACSCINIT  = $(BLACSDIR)/libblacsCinit$(BLACSPOSTFIX).a
BLACSLIB    = $(BLACSDIR)/libblacs$(BLACSPOSTFIX).a
CBLACSLIB   = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
FBLACSLIB   = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

# Scalapack libraries
#distribution compiled version
#SCALAPACKDIR = $(USRDIR)
#user compiled version
 SCALAPACKDIR = $(HOME)/usr

#distribution compiled version
#SCLPCKPOSTFIX = -openmpi
#user compiled version
 SCLPCKPOSTFIX =
 
SCALAPACKLIB = $(SCALAPACKDIR)/lib/libscalapack$(SCLPCKPOSTFIX).a

# Parallel libraries
PLIBS = $(SCALAPACKLIB) $(CBLACSLIB) -lmpich

#eof

