
PLT = Linux_x8664

USR=/usr
ATLAS_DIR=$(USR)
BLACS_LIBDIR = $(HOME)$(USR)/lib
BLAS_DIR=$(USR)
BOOST_DIR = $(HOME)$(USR)
# .. because cygwin comes with boost 1.33
BOOST_VER = 
FFTW_DIR = $(USR)
FFTW3_DIR = $(USR)
HDF5_DIR = $(USR)
MPI_IMPL = 
MPI_DIR = $(USR)/lib/$(MPI_IMPL)
LAPACK_DIR = $(USR)
SCALAPACK_LIBDIR = $(HOME)$(USR)/lib
XERCESC_DIR = $(USR)

FFTW_LIB = $(FFTW_DIR)/lib/libfftw.a
FFTW3_LIB = $(FFTW3_DIR)/lib/libfftw3.a
BLAS_LIB = $(BLAS_DIR)/lib/libblas.a
BOOST_DIR = $(HOME)/usr
BOOST_SYSTEM_LIB = $(BOOST_DIR)/lib/libboost_system.a
BOOST_FILESYSTEM_LIB = $(BOOST_DIR)/lib/libboost_filesystem.a
BOOST_MPI_LIB = $(BOOST_DIR)/lib/libboost_mpi.a
LAPACK_LIB = $(BOOST_DIR)/lib/liblapack.a
HDF5_LIB = $(HDF5_DIR)/lib/libhdf5.a
Z_LIB = $(HDF5_DIR)/lib/libz.a
MPI_F77_LIB = $(MPI_DIR)/lib/libmpi_f77.a
GFORTRAN_LIB = /usr/lib/gcc/i686-pv-cygwin/4.3.4/libgfortran.a

CXXDIR =
CXX = $(HOME)/usr/bin/mpicxx
LD = $(CXX) 

PLTFLAGS += -DIA32 -DUSE_FFTW -DUSE_XERCES -DUSE_MPI -DADD_ -DSCALAPACK -DAPP_NO_THREADS -DXML_USE_NO_THREADS

INCLUDE_RAW = \
	-I$(BOOST_DIR)/include \
	-I$(FFTW_DIR)/include \
	-I$(FFTW3_DIR)/include \
	-I$(HDF5_DIR)/include \
	-I$(MPI_DIR)/include \
	-I$(XERCESC_DIR)/include \
	-I./include

#INCLUDE=$(INCLUDE_RAW)
INCLUDE=$(shell awk 'BEGIN{split("${INCLUDE_RAW}",a);for(i in a)b[a[i]]=1;r="";for(i in b)r=r" "i;print r}')

CXXFLAGS = \
	-O4 -D$(PLT) \
	$(INCLUDE) \
	$(PLTFLAGS) \
	$(DFLAGS) \
	-Wno-write-strings -Wfatal-errors

LIBPATH_RAW = \
	-L$(BOOST_DIR)/lib \
	-L$(FFTW_DIR)/lib \
	-L$(MPI_DIR)/lib \
	-L$(LAPACK_DIR)/lib \
	-L$(BLAS_DIR)/lib \
	-L$(ATLAS_DIR)/lib \
	-L$(XERCESC_DIR)/lib \
	-L$(HDF5_DIR)/lib

#LIBPATH = $(LIBPATH_RAW)

LIBPATH=$(shell awk 'BEGIN{split("${LIBPATH_RAW}",a);for(i in a)b[a[i]]=1;r="";for(i in b)r=r" "i;print r}')

LIBS = \
	$(PLIBS) \
	-lboost_filesystem$(BOOST_VER) -lboost_system$(BOOST_VER) \
	-lfftw -lfftw3 \
	-lhdf5 \
	-llapack -lblas \
	-lxerces-c \
	-lgfortran
#	-lmpi_f77 -lmpi_cxx 

LDFLAGS = $(LIBPATH) $(LIBS)

PLAT=$(MACHTYPE)

BLACS_VER = -mpi
BLACSFINIT  = $(BLACS_LIBDIR)/libblacsF77init$(BLACS_VER).a
BLACSCINIT  = $(BLACS_LIBDIR)/libblacsCinit$(BLACS_VER).a
BLACSLIB    = $(BLACS_LIBDIR)/libblacs$(BLACS_VER).a
CBLACSLIB   = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
FBLACSLIB   = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

SCALAPACK_VER = $(MPI_IMPL)
SCALAPACK_LIB = $(SCALAPACK_LIBDIR)/libscalapack$(SCALAPACK_VER).a

PLIBS = $(SCALAPACK_LIB) $(CBLACSLIB)
#eof

