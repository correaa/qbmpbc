#knonw to work with ubuntu 10.10 and qbox-1.52.3 use as 'export TARGET=ubuntu', 'make ubuntu', 'make qb lib'
all: qb lib lib_depend

ubuntu: 
	sudo apt-get install \
          fftw-dev libxerces-c-dev libscalapack-mpi-dev g++ \
          libatlas-dev liblapack-dev libblas-dev libblacs-mpi-dev \
          uuid-dev libopenmpi-dev
#last 2 lines above are redundant

PLT=$(shell uname -m)
USRDIR=/usr/lib

BOOSTDIR=$(USRDIR)
BOOST_FILESYSTEM_NAME=boost_filesystem
MPIDIR=$(USRDIR)

XERCESCQ=USE_XERCES
XERCESCDIR=$(USRDIR)
XERCESCFLAG=XERCESC_3_0_1

FFTWDIR=$(USRDIR)
FFTWNAME=fftw
BLASDIR=$(USRDIR)
LAPACKDIR=$(USRDIR)

FFTW3DIR=$(USRDIR)

HDF5DIR=$(USRDIR)
HDF5NAME=hdf5

CXXDIR=
CXX=$(CXXDIR)mpicxx$(MPIIMPL)
LD=$(CXX) 
PLTFLAGS += \
    -DUSE_FFTW \
    -D$(XERCESCQ) -D$(XERCESCFLAG) \
    -DUSE_MPI -DADD_ -DSCALAPACK -DAPP_NO_THREADS -DXML_USE_NO_THREADS
INCLUDE = -Iinclude/

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
	-lfftw3 -lboost_system -lboost_filesystem -l$(HDF5NAME)

LDFLAGS = $(LIBPATH) $(LIBS) -Wl,-rpath=$(HOME)/usr/lib 

ARCHIVE_LIBS = 

BLACSDIR = /usr/lib
BLACSFINIT_NAME = blacsF77init-openmpi
BLACSCINIT_NAME  = blacsCinit-openmpi
BLACSLIB_NAME = blacs-openmpi
BLACSFINIT  = -l$(BLACSFINIT_NAME)
BLACSCINIT  = -l$(BLACSCINIT_NAME)
BLACSLIB    = -l$(BLACSLIB_NAME)
CBLACSLIB   = -L$(BLACSDIR) $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
FBLACSLIB   = -L$(BLACSDIR) $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

SCALAPACKDIR = /usr/lib
SCALAPACKVER = -openmpi
SCALAPACKNAME = scalapack$(SCALAPACKVER)
SCALAPACKLIB = -l$(SCALAPACKNAME)
PLIBS = $(SCALAPACKLIB) $(CBLACSLIB) 

ARCHIVE_LIBS = \
 $(BLACSDIR)/lib$(BLACSFINIT_NAME).a  $(BLACSDIR)/lib$(BLACSCINIT_NAME).a $(BLACSDIR)/lib$(BLACSLIB_NAME).a \
# $(SCALAPACKDIR)/lib$(SCALAPACKNAME).a
# $(BOOSTDIR)/lib$(BOOST_FILESYSTEM_NAME).a
# $(FFTWDIR)/lib$(FFTWNAME).a
# $(HDF5DIR)/lib$(HDF5NAME).a \
# $(FFTW_LIB) $(MPI_F77_LIB) \
# $(PLIBS) \
# $(XERCESC_LIB) \
# $(BOOST_MPI_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_SERIALIZATION_LIB) $(BOOST_SYSTEM_LIB)


#warning: both boost_mpi and fftw3 have a timer.o object file inside
#warning: this assumes that /tmp/ is writable
lib_depend: 
	rm -rf /tmp/qb_depend && \
	mkdir -p /tmp/qb_depend && cd /tmp/qb_depend && \
	for a in $(ARCHIVE_LIBS); do echo $$a; ar -x $$a; done &&\
	rm -f libqb_depend.a libqb_depend.so && \
	pwd && \
	ar -cr libqb_depend.a *.o *.o && \
	gcc -fPIC -shared *.o -o libqb_depend.so && \
	cd - && mv /tmp/qb_depend/libqb_depend.a . 

