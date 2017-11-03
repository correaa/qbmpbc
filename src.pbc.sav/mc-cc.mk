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
#  mc-cc.mk
#
#-------------------------------------------------------------------------------
# $Id: mc-cc.mk,v 1.11 2008/06/18 03:39:53 fgygi Exp $
#
 PLT=mc-cc
#-------------------------------------------------------------------------------
 MPIDIR=/opt/mpich/gnu
 BOOSTDIR=$(HOME)/usr
 XERCESCDIR=$(HOME)/usr
 FFTWDIR=/opt/fftw-2.1.5
 BLASDIR=/usr/lib
 LAPACKDIR=/usr/lib

 PLTOBJECTS = readTSC.o

 #CXX=/usr/bin/g77
 #CXX=/usr/bin/g++
 CXX=/opt/mpich/gnu/bin/mpiCC
 LD=$(CXX)

# PLTFLAGS += -DIA32 -DUSE_FFTW -D_LARGEFILE_SOURCE \
#             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DADD_ \
#             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES

 PLTFLAGS += -DIA32 -DUSE_FFTW -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=32 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES


 INCLUDE = -I./include -I$(XERCESCDIR)/include \
           -I$(MPIDIR)/include -I$(BOOSTDIR)/include -I$(FFTWDIR)/include

 CXXFLAGS= -g -O4 $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR)/lib -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L$(LAPACKDIR) -L$(BLASDIR) \
           -L$(XERCESCDIR)/lib -L/export/apps/matlab/bin/glnx86 -L$(BOOSTDIR)/lib

 LIBS =  $(PLIBS) -lfftw -lfftw3 \
          -llapack -lblas -lm \
          -lxerces-c -L/export/apps/intel/fc/9.1.043/lib -lifcore \
          -lboost_system-gcc34-mt -lboost_filesystem-gcc34-mt

 LDFLAGS = $(LIBPATH) $(LIBS)

 PLAT=LINUX
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = $(HOME)/Codes/oldfiles/BLACS/LIB
 BLACSFINIT    = $(BLACSdir)/blacsF77init_MPI-$(PLAT)-$(BLACSDBGLVL).a
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLAT)-$(BLACSDBGLVL).a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLAT)-$(BLACSDBGLVL).a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = $(HOME)/Codes/oldfiles/scalapack-1.8.0
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

# Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
