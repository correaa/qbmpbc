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
#  cygwin.mk
#
#-------------------------------------------------------------------------------
# $Id: x8664_gcc.mk,v 1.11 2008/06/18 03:39:53 fgygi Exp $
#
 PLT=cywgin
#-------------------------------------------------------------------------------
 HOME=/home/caiwei
 MPIDIR=$(HOME)/soft/MPICH2
 XERCESCDIR=$(HOME)/usr
 FFTWDIR=$(HOME)/usr
 BLASDIR=/usr/lib
 LAPACKDIR=/usr/lib
 HDF5DIR=$(HOME)/usr
 FFTW3DIR=$(HOME)/usr

 PLTOBJECTS = readTSC.o

 CXX=/usr/bin/g++
 LD=$(CXX) 

# cygwin does not use SCALAPACK
 PLTFLAGS += -DIA32 -DUSE_FFTW -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DADD_ -DAdd_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES


 INCLUDE = -I./include -I$(MPIDIR)/include -I$(FFTWDIR)/include -I$(XERCESCDIR)/include 
 CXXFLAGS= -g -O4 -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS) -Wno-write-strings 

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR)/lib -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L$(LAPACKDIR)/lib -L$(BLASDIR)/lib \
           -L$(XERCESCDIR)/lib -L$(HDF5DIR)/lib

 LIBS =  $(PLIBS) -lfftw \
         -llapack -lblas -lm \
         /home/caiwei/usr/lib/libxerces-c.a \
         -Xlinker -Bdynamic -lsicuuc -lsicudata -lmpi

 LDFLAGS = $(LIBPATH) $(LIBS)

 PLAT=Linux_x8664
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = /usr/lib
 BLACSFINIT    = $(BLACSdir)/libblacsF77init-mpich.a
 BLACSCINIT    = $(BLACSdir)/libblacsCinit-mpich.a
 BLACSLIB      = $(BLACSdir)/libblacs-mpich.a

 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)
 FBLACSLIB     = $(BLACSFINIT) $(BLACSLIB) $(BLACSFINIT)

 # Scalapack libraries
 SCALAPACK_DIR = /usr/lib
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack-mpich.a

# Parallel libraries
# PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)
 PLIBS =

#-------------------------------------------------------------------------------
