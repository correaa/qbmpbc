#-------------------------------------------------------------------------------
#
# Copyright (c) 2008 The Regents of the University of California
#
# This file is part of Qbox
#
# Qbox is distributed under the terms of the GNU General Public License 
# as published by the Free Software Foundation, either version 2 of 
# the License, or (at your option) any later version.
# See the file COPYING in the root directory of this distribution
# or <http://www.gnu.org/licenses/>.
#
#-------------------------------------------------------------------------------
#
#  wcr.mk
#
#-------------------------------------------------------------------------------
# $Id: wcr.mk,v 1.11 2008/06/18 03:39:53 fgygi Exp $
#
 PLT=wcr
#-------------------------------------------------------------------------------
 MPIDIR=/share/apps/mvapich/intel

 BOOSTDIR=$(HOME)/usr
 BOOST_VER=-gcc34-mt
 XERCESCDIR=$(HOME)/usr
 FFTWDIR=$(HOME)/usr
 BLASDIR=/usr/lib64
 LAPACKDIR=/usr/lib64
 FFTW3DIR=$(HOME)/usr
 HDF5DIR=$(HOME)/usr

 PLTOBJECTS = readTSC.o

 CXX=/share/apps/mvapich/intel/bin/mpicxx
 LD=$(CXX)

 PLTFLAGS += -DIA32 -DUSE_FFTW -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES


 INCLUDE = -I./include -I$(XERCESCDIR)/include \
           -I$(MPIDIR)/include -I$(FFTWDIR)/include -I$(FFTW3DIR)/include -I$(HDF5DIR)/include


 CXXFLAGS= -O  -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)
#CXXFLAGS= -g  -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR) -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L$(LAPACKDIR) -L$(BLASDIR) \
           -L$(XERCESCDIR)/lib \
           -L$(BOOSTDIR)/lib -L$(HDF5DIR)/lib -L/opt/intel/fce/9.1.052/lib

 LIBS =  $(PLIBS) -lfftw -lfftw3 \
          -llapack -lblas -lm \
          -lxerces-c -lc \
          -lboost_system$(BOOST_VER) -lboost_filesystem$(BOOST_VER) -lhdf5 -lifcore \

 LDFLAGS = $(LIBPATH) $(LIBS)

 PLAT=LINUX_x8664
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = $(HOME)/usr
 BLACSFINIT    = $(BLACSdir)/lib/libblacsF77init-mpi.a
 BLACSCINIT    = $(BLACSdir)/lib/libblacsCinit-mpi.a
 BLACSLIB      = $(BLACSdir)/lib/libblacs-mpi.a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = $(HOME)/usr
 SCALAPACKLIB  = $(SCALAPACK_DIR)/lib/libscalapack.a

 #BLACSdir      = $(HOME)/Codes/BLACS/LIB
 #BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-$(PLAT)-$(BLACSDBGLVL).a
 #BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLAT)-$(BLACSDBGLVL).a
 #BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLAT)-$(BLACSDBGLVL).a

 # Scalapack libraries
 #SCALAPACK_DIR = $(HOME)/Codes/scalapack-1.8.0
 #SCALAPACKLIB  = -L$(SCALAPACK_DIR) -lscalapack

# Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
