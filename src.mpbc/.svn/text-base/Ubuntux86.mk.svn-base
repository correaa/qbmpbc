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
#  Ubuntux86.mk
#
#-------------------------------------------------------------------------------
# $Id: Ubuntux86.mk,v 1.11 2008/10/18 03:39:53 correaa Exp $
#
 PLT=Linux_x8664
#-------------------------------------------------------------------------------
 MPIDIR=/usr/lib/mpich
 XERCESCDIR=/usr/lib
 FFTWDIR=/usr/lib
 BLASDIR=/usr/lib
 LAPACKDIR=/usr/lib

 PLTOBJECTS = readTSC.o

 CXX=g++-4.2
 LD=$(CXX)

 PLTFLAGS += -DIA32 -DUSE_FFTW -D_LARGEFILE_SOURCE \
             -D_FILE_OFFSET_BITS=64 -DUSE_MPI -DSCALAPACK -DADD_ \
             -DAPP_NO_THREADS -DXML_USE_NO_THREADS -DUSE_XERCES


 INCLUDE = -I$(MPIDIR)/include -I$(FFTWDIR) -I$(XERCESCDIR)/include

 CXXFLAGS= -g -O4 -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(GCCDIR)/lib -L$(FFTWDIR)/.libs -L/usr/X11R6/lib \
           -L$(MPIDIR)/lib -L$(LAPACKDIR) -L$(BLASDIR) \
           -L$(XERCESCDIR)/lib

 LIBS =  $(PLIBS) -lfftw -lfftw3 -lboost_filesystem \
         -llapack -lf77blas -latlas -lm \
         -Xlinker -Bstatic \
          -lc -lgfortran -static-libgcc -lmpich \
         -Xlinker -Bdynamic -lxerces-c -lrt

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
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)

#-------------------------------------------------------------------------------
