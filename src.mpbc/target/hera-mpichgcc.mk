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
#  hera_gcc.mk
#
#-------------------------------------------------------------------------------
# $Id: x8664_gcc.mk,v 1.15 2009/12/04 03:16:03 fgygi Exp $
#
 PLT=Linux_x8664
#-------------------------------------------------------------------------------
 MPIDIR=/opt/mvapich-gnu-gen2-0.9.9/
 XERCESCDIR=$(HOME)/usr
 FFTWDIR=/usr/lib64
 BLASDIR=/usr/lib64
 LAPACKDIR=$(HOME)/usr/lib-hera

 PLTOBJECTS = readTSC.o
 CXX_DIR= $(MPIDIR)bin/
 CXX=$(CXX_DIR)mpiCC
 LD=$(CXX)

 PLTFLAGS += -DIA32 -DUSE_FFTW -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES 
#            -DXERCESC_3_0_1


 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include

 CXXFLAGS= -g -Wunused -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR)/lib -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L$(LAPACKDIR) -L$(BLASDIR) \
           -L$(XERCESCDIR)/lib
 LIBS = $(PLIBS) -lblas -lfftw -lxerces-c -llapack
# LIBS =  $(PLIBS) -lpthread -lfftw \
#         -Xlinker -Bstatic \
#          -lc -lgfortran -static-libgcc -lmpich -lxerces-c -luuid \
#         -Xlinker -Bdynamic

 LDFLAGS = $(LIBPATH) $(LIBS) -Wl,-rpath=/opt/mvapich-gnu-gen2-0.9.9/lib/shared

 PLAT=Linux_x8664
 # Blacs libraries
  
 BLACSDBGLVL   = 0
 BLACSdir      = $(HOME)/usr/lib-hera
 BLACSFINIT    = $(BLACSdir)/libblacsF77init-mpi.a
 BLACSCINIT    = $(BLACSdir)/libblacsCinit-mpi.a
 BLACSLIB      = $(BLACSdir)/libblacs-mpi.a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = $(HOME)/usr/lib-hera
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

# Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
