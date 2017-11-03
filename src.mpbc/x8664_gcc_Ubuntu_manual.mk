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
MPIDIR=$(HOME)/usr
XERCESCDIR=$(HOME)/usr
FFTWDIR=$(HOME)/usr
BLASDIR=$(HOME)/usr
ATLASDIR=$(HOME)/usr
LAPACKDIR=$(HOME)/usr
HDF5DIR=$(HOME)/usr
BOOSTDIR=$(HOME)/usr
BOOST_VER = 
FFTW3DIR=$(HOME)/usr
BLACSDIR=$(HOME)/usr/lib
SCALAPACKDIR = $(HOME)/usr

CXXDIR=
CXX=c++
LD=$(CXX) 

PLTFLAGS += -DIA32 -DUSE_FFTW -DUSE_XERCES -DUSE_MPI -DADD_ -DSCALAPACK -DAPP_NO_THREADS -DXML_USE_NO_THREADS

INCLUDE_RAW = \
	-I$(BOOSTDIR)/include \
	-I$(FFTWDIR)/include \
	-I$(FFTW3DIR)/include \
	-I$(HDF5DIR)/include \
	-I$(MPIDIR)/include \
	-I$(XERCESCDIR)/include
INCLUDE=$(shell awk 'BEGIN{split("${INCLUDE_RAW}",a);for(i in a)b[a[i]]=1;r="";for(i in b)r=r" "i;print r}')

CXXFLAGS = \
	-O4 -D$(PLT) \
	$(INCLUDE) \
	$(PLTFLAGS) \
	$(DFLAGS) \
	-Wno-write-strings -Wfatal-errors

LIBPATH_RAW = \
	-L$(BOOSTDIR)/lib \
	-L$(FFTWDIR)/lib \
	-L$(MPIDIR)/lib \
	-L$(LAPACKDIR)/lib \
	-L$(BLASDIR)/lib \
	-L$(ATLASDIR)/lib \
	-L$(XERCESCDIR)/lib \
	-L$(HDF5DIR)/lib
LIBPATH = $(shell awk 'BEGIN{split(${LIBPATH_RAW},a);for(i in a)b[a[i]]=1;r="";for(i in b)r=r" "i;print r}')

LIBS = \
	$(PLIBS) \
	-lboost_filesystem$(BOOST_VER) -lboost_system$(BOOST_VER) \
	-lfftw \
	-lfftw3 \
	-lhdf5 \
	-lmpich -lmpl \
	-llapack -lblas \
	-lgfortran \
	-lxerces-c

LDFLAGS = $(LIBPATH) $(LIBS) 
#-Wl,-rpath=$(HOME)/usr/lib 

PLAT=$(MACHTYPE)

BLACSPOSTFIX  = -mpi

BLACSFINIT  = $(BLACSDIR)/libblacsF77init$(BLACSPOSTFIX).a
BLACSCINIT  = $(BLACSDIR)/libblacsCinit$(BLACSPOSTFIX).a
BLACSLIB    = $(BLACSDIR)/libblacs$(BLACSPOSTFIX).a
CBLACSLIB   = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
FBLACSLIB   = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

SCLPCKPOSTFIX = 
SCALAPACKLIB = $(SCALAPACKDIR)/lib/libscalapack$(SCLPCKPOSTFIX).a

PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)
#eof
