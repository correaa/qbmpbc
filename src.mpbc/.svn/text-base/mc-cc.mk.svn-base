#-------------------------------------------------------------------------------
# Notes: gcc-34 version too low to compile driver program
#        /opt/mpich/intel/bin/mpicxx compiler lacks appropriate license
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
#  mc-cc.mk
#
#-------------------------------------------------------------------------------
# $Id: mc-cc.mk,v 1.11 2008/06/18 03:39:53 fgygi Exp $
#
 PLT=mc-cc
#-------------------------------------------------------------------------------
 MPIDIR=/opt/mpich/intel
#MPIDIR=/opt/mpich/gnu

 BOOSTDIR=$(HOME)/usr
 BOOST_VER=-gcc34-mt
 XERCESCDIR=$(HOME)/usr
 FFTWDIR=/opt/fftw-2.1.5
 BLASDIR=/usr/lib
 LAPACKDIR=/usr/lib
 FFTW3DIR=$(HOME)/usr
 HDF5DIR=$(HOME)/usr

 PLTOBJECTS = readTSC.o

#CXX=/usr/bin/g77
#CXX=/usr/bin/g++
#CXX=/opt/mpich/gnu/bin/mpiCC
#CXX=/opt/mpich/intel/bin/mpicc
 CXX=icc
 LD=$(CXX)

# PLTFLAGS += -DIA32 -DUSE_FFTW -D_LARGEFILE_SOURCE \
#             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DADD_ \
#             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES

 PLTFLAGS += -DIA32 -DUSE_FFTW -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=32 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES

 INCLUDE = -I./include -I$(XERCESCDIR)/include \
           -I$(MPIDIR)/include -I$(FFTWDIR)/include -I$(FFTW3DIR)/include -I$(HDF5DIR)/include

 CXXFLAGS= -O $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR)/lib -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L$(LAPACKDIR) -L$(BLASDIR) \
           -L$(XERCESCDIR)/lib \
           -L$(BOOSTDIR)/lib -L$(HDF5DIR)/lib -L/export/apps/intel/fc/9.1.043/lib

 LIBS =  $(PLIBS) -lfftw -lfftw3 \
          -llapack -lblas -lm \
          -lxerces-c -lc \
          -lboost_system$(BOOST_VER) -lboost_filesystem$(BOOST_VER) -lhdf5 -lifcore -lmpich \

# LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR)/lib -L/usr/X11R6/lib \
#           -L$(MPIDIR)/lib -L$(LAPACKDIR) -L$(BLASDIR) \
#           -L$(XERCESCDIR)/lib -L/export/apps/matlab/bin/glnx86 -L$(BOOSTDIR)/lib

 LDFLAGS = $(LIBPATH) $(LIBS)

 PLAT=LINUX
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = $(HOME)/usr
 BLACSFINIT    = $(BLACSdir)/lib/libblacsF77init-mpi.a
 BLACSCINIT    = $(BLACSdir)/lib/libblacsCinit-mpi.a
 BLACSLIB      = $(BLACSdir)/lib/libblacs-mpi.a

# BLACSdir      = $(HOME)/Codes/oldfiles/BLACS/LIB
# BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-$(PLAT)-$(BLACSDBGLVL).a
# BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLAT)-$(BLACSDBGLVL).a
# BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLAT)-$(BLACSDBGLVL).a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = $(HOME)/usr
 SCALAPACKLIB  = $(SCALAPACK_DIR)/lib/libscalapack.a

# SCALAPACK_DIR = $(HOME)/Codes/oldfiles/scalapack-1.8.0
# SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

# Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
