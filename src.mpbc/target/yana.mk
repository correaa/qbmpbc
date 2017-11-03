all: qb lib lib_depend

PLT=Linux_x8664

USR=$(HOME)/usr

ATLAS_DIR=$(USR)
BLACS_LIBDIR = $(USR)/lib$(CLUSTER)
BLAS_LIBDIR= /usr/lib64
BOOST_DIR = $(USR)
BOOST_LIBDIR = $(USR)/lib$(CLUSTER)
BOOST_VER = 
FFTW_DIR = $(USR)
FFTW3_LIBDIR = /usr/lib64
HDF5_DIR = $(USR)
#MPI_IMPL = openmpi
#MPI_DIR = $(USR)/lib$(CLUSTER)/$(MPI_IMPL)
#MPI_LIBDIR = $(MPI_DIR)/lib
LAPACK_LIBDIR = $(BLAS_LIBDIR)
SCALAPACK_LIBDIR = $(USR)/lib$(CLUSTER)
XERCESC_DIR = $(USR)

FFTW_LIB = $(FFTW_DIR)/lib$(CLUSTER)/libfftw.a
FFTW3_LIB = $(FFTW3_LIBDIR)/libfftw3.a
BLAS_LIB = $(BLAS_LIBDIR)/libblas.a
BOOST_SYSTEM_LIB = $(BOOST_DIR)/lib$(CLUSTER)/libboost_system.a
BOOST_FILESYSTEM_LIB = $(BOOST_DIR)/lib$(CLUSTER)/libboost_filesystem.a
BOOST_MPI_LIB = $(BOOST_DIR)/lib$(CLUSTER)/libboost_mpi.a
BOOST_SERIALIZATION_LIB = $(BOOST_DIR)/lib$(CLUSTER)/libboost_serialization.a
LAPACK_LIB = $(LAPACK_LIBDIR)/liblapack.a
HDF5_LIB = $(HDF5_DIR)/lib$(CLUSTER)/libhdf5.a
Z_LIB = /usr/lib64/libz.a
#MPI_F77_LIB = $(MPI_DIR)/lib/libmpi_f77.a
XERCESC_LIB = $(XERCESC_DIR)/lib$(CLUSTER)/libxerces-c.a
#GFORTRAN_LIB = /usr/lib/gcc/x86_64-redhat-linux/4.1.1/libgfortran.a

CXX_DIR=
CXX=$(CXX_DIR)mpiCC
LD=$(CXX) 

PLTFLAGS += -DIA32 -DUSE_FFTW -DUSE_XERCES -DUSE_MPI -DADD_ -DSCALAPACK -DAPP_NO_THREADS -DXML_USE_NO_THREADS
#other PLTFLAGS -DXERCES_3_0_1

INCLUDE_RAW = \
	-I$(BOOST_DIR)/include \
	-I$(FFTW_DIR)/include \
	-I$(FFTW3_DIR)/include \
	-I$(HDF5_DIR)/include \
	-I$(MPI_DIR)/include \
	-I$(XERCESC_DIR)/include \
	-I./include
INCLUDE=$(INCLUDE_RAW)
#INCLUDE=$(shell awk 'BEGIN{split("${INCLUDE_RAW}",a);for(i in a)b[a[i]]=1;r="";for(i in b)r=r" "i;print r}')

CXXFLAGS = \
	-O4 -D$(PLT) \
	$(INCLUDE) \
	$(PLTFLAGS) \
	$(DFLAGS) \
	-Wno-write-strings -Wfatal-errors

LIBPATH_RAW = \
	-L$(BOOST_LIBDIR) \
	-L$(FFTW_LIBDIR) \
	-L$(MPI_LIBDIR) \
	-L$(BLAS_LIBDIR) \
	-L$(LAPACK_LIBDIR) \
	-L$(XERCESC_DIR)/lib \
	-L$(HDF5_DIR)/lib
LIBPATH=$(LIBPATH_RAW)
#LIBPATH=$(shell awk 'BEGIN{split(${LIBPATH_RAW},a);for(i in a)b[a[i]]=1;r="";for(i in b)r=r" "i;print r}')

LIBS = \
	$(PLIBS) \
	-lboost_filesystem$(BOOST_VER) -lboost_system$(BOOST_VER) \
	-lfftw \
	-lfftw3 \
	-lhdf5 -lz \
	-llapack -lblas \
	-lgfortran \
	-lxerces-c

#-Wl,-rpath=$(HOME)/usr/lib 

PLAT=$(MACHTYPE)

BLACSPOSTFIX  = -mpi

BLACSFINIT  = $(BLACS_LIBDIR)/libblacsF77init$(BLACSPOSTFIX).a
BLACSCINIT  = $(BLACS_LIBDIR)/libblacsCinit$(BLACSPOSTFIX).a
BLACS_LIB    = $(BLACS_LIBDIR)/libblacs$(BLACSPOSTFIX).a
CBLACS_LIB   = $(BLACSCINIT) $(BLACS_LIB) $(BLACSCINIT)
FBLACS_LIB   = $(BLACSFINIT) $(BLACS_LIB) $(BLACSFINIT)

SCLPCKPOSTFIX = 
SCALAPACK_LIB = $(SCALAPACK_LIBDIR)/libscalapack$(SCLPCKPOSTFIX).a

PLIBS = $(SCALAPACK_LIB) $(CBLACS_LIB)

ARCHIVE_LIBS = \
 $(BLACSCINIT) $(BLACS_LIB) $(SCALAPACK_LIB) \
 $(HDF5_LIB) \
 $(FFTW_LIB) $(MPI_F77_LIB) \
 $(PLIBS) \
 $(XERCESC_LIB) \
 $(BOOST_MPI_LIB) $(BOOST_FILESYSTEM_LIB) $(BOOST_SERIALIZATION_LIB) $(BOOST_SYSTEM_LIB)
#$(Z_LIB) \
#$(LAPACK_LIB) \
#$(BLAS_LIB)
#$(GFORTRAN_LIB) 
LDFLAGS = $(LIBPATH) $(LIBS)

lib_depend: 
	rm -rf /tmp/qb_depend && \
	mkdir -p /tmp/qb_depend && cd /tmp/qb_depend && \
	for a in $(ARCHIVE_LIBS); do echo $$a; ar -x $$a; done &&\
	rm -f libqb_depend.a libqb_depend.so && \
	pwd && \
	ar -cr libqb_depend.a *.o *.o && \
	gcc -fPIC -shared *.o -o libqb_depend.so && \
	cd - && mv /tmp/qb_depend/libqb_depend.a . 


#eof
