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
#  su-ahpcrc.mk
#
#-------------------------------------------------------------------------------
# $Id: wcr.mk,v 1.11 2008/06/18 03:39:53 fgygi Exp $
#
 PLT=suahpcrc
#-------------------------------------------------------------------------------
 MPIDIR=/usr/mpi/intel/mvapich2-1.2

 BOOSTDIR=$(HOME)/usr-intel
 BOOST_VER=
 XERCESCDIR=$(HOME)/usr-intel
 FFTWDIR=$(HOME)/usr-intel
 BLASDIR=$(HOME)/usr-intel/local/blas
 LAPACKDIR=$(HOME)/usr-intel/local/lapack
 FFTW3DIR=$(HOME)/usr-intel
 HDF5DIR=$(HOME)/usr-intel

 PLTOBJECTS = readTSC.o

 CXX=mpicxx
 LD=$(CXX) 

 PLTFLAGS += -DMPICH_IGNORE_CXX_SEEK -DIA32 -DUSE_FFTW -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK  -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES


 INCLUDE = -I./include -I$(XERCESCDIR)/include \
           -I$(MPIDIR)/include -I$(FFTWDIR)/include -I$(FFTW3DIR)/include -I$(HDF5DIR)/include

 CXXFLAGS= -O  -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)
#CXXFLAGS= -g  -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR)/lib -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L$(LAPACKDIR)/lib -L$(BLASDIR)/lib \
           -L$(XERCESCDIR)/lib \
           -L$(BOOSTDIR)/lib -L$(HDF5DIR)/lib

 LIBS =  $(PLIBS) -lfftw -lfftw3 \
          -llapack -lblas -lm \
          -lxerces-c -L/opt/intel/fce/10.1.015/lib/ -lifcore \
          -lboost_system$(BOOST_VER) -lboost_filesystem$(BOOST_VER) -lhdf5\

 LDFLAGS = $(LIBPATH) $(LIBS)

 PLAT=Linux_x8664
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = $(HOME)/usr-intel/local/blacs
 BLACSFINIT    = $(BLACSdir)/lib/libmpiblacsF77init.a
 BLACSCINIT    = $(BLACSdir)/lib/libmpiblacsCinit.a
 BLACSLIB      = $(BLACSdir)/lib/libmpiblacs.a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = $(HOME)/usr-intel/local/scalapack
 SCALAPACKLIB  = $(SCALAPACK_DIR)/lib/libscalapack.a

# Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
