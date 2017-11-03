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
#  x8664_gcc.mk
#
#-------------------------------------------------------------------------------
# $Id: x8664_gcc.mk,v 1.11 2008/06/18 03:39:53 fgygi Exp $
#
 PLT=Linux_x8664
#-------------------------------------------------------------------------------
 MPIDIR=$(HOME)/usr
 XERCESCDIR=$(HOME)/usr
 FFTWDIR=$(HOME)/usr
 BLASDIR=$(HOME)/usr/local/atlas
 LAPACKDIR=$(HOME)/usr/local/lapack

 FFTW3DIR=$(HOME)/usr

 PLTOBJECTS = readTSC.o

 CXX=/usr/bin/g++
 LD=$(CXX)

 PLTFLAGS += -DIA32 -DUSE_FFTW -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES


 INCLUDE = -I./include -I$(MPIDIR)/include -I$(FFTWDIR)/include -I$(XERCESCDIR)/include -I$(FFTW3DIR)/include

 CXXFLAGS= -g -O4 -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR)/lib -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L$(LAPACKDIR)/lib -L$(BLASDIR)/lib \
           -L$(XERCESCDIR)/lib

 LIBS =  $(PLIBS) -lfftw -lfftw3\
         -llapack -lf77blas -latlas -lm -lboost_filesystem-gcc41-mt -lboost_system-gcc41-mt \
         -Xlinker -Bstatic \
          -lc -lgfortran -static-libgcc -lmpich -lxerces-c -lrt -pthread \
         -Xlinker -Bdynamic -lcurl

 LDFLAGS = $(LIBPATH) $(LIBS)

 PLAT=Linux_x8664
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = $(HOME)/usr/local/blacs/lib
 BLACSFINIT    = $(BLACSdir)/libmpiblacsF77init.a
 BLACSCINIT    = $(BLACSdir)/libmpiblacsCinit.a
 BLACSLIB      = $(BLACSdir)/libmpiblacs.a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = $(HOME)/usr/local/scalapack/lib
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

# Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
