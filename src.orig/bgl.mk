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
#  bgl-anl.mk
#
#-------------------------------------------------------------------------------
# $Id: $
#
 PLT=BGL
#-------------------------------------------------------------------------------
 BGL_ROOT=/bgl/BlueLight/ppcfloor
 BGL_SYS=$(BGL_ROOT)/bglsys

 LIBS_MPI     += -L $(BGL_ROOT)/bglsys/lib -lmpich.rts \
                 -lmsglayer.rts -lrts.rts -ldevices.rts \
                 -lc -lnss_files -lnss_dns -lresolv

 GNU_ROOT=/BlueLight/ppcfloor
 BLRTS_GNU_ROOT=$(GNU_ROOT)/blrts-gnu
 CXX=blrts_xlC

 LD=$(CXX)

 PLTFLAGS += -DUSE_FFTW
 PLTFLAGS += -DUSE_MPI -DSCALAPACK
 PLTFLAGS +=  -D__linux__ -DPLT_BIG_ENDIAN
 PLTFLAGS += -DUSE_XERCES
 PLTFLAGS += -DUSE_CSTDIO_LFS -D_LARGEFILE64_SOURCE  -D_FILE_OFFSET_BITS=64
 PLTFLAGS += -DMPICH_IGNORE_CXX_SEEK

 FFTWDIR=$(HOME)/software/fftw/bgl/fftw-2.1.5
 FFTWINCLUDEDIR=$(FFTWDIR)/fftw
 FFTWLIBDIR=$(FFTWDIR)/fftw

 XERCESCDIR=$(HOME)/software/xml/xerces-c-src_2_6_0
 XERCESCLIBDIR=$(XERCESCDIR)/lib

 BLASLIB=-L/soft/tools/GotoBLAS -lgoto

 INCLUDE =  -I$(XERCESCDIR)/include \
            -I$(FFTWINCLUDEDIR) -I$(BGL_ROOT)/bglsys/include

 CXXFLAGS= -g -O3 -qarch=440d -D$(PLT) $(INCLUDE) $(PLTFLAGS) $(DFLAGS)

 LIBPATH = -L$(FFTWLIBDIR) -L$(XERCESCLIBDIR) \
           -L/opt/ibmcmp/xlf/bg/10.1/blrts_lib

 LIBS =  $(PLIBS) -lfftw $(DGEMMLIB) $(BLASLIB) -lg2c \
         -lxlf90 -lxlopt -lxlomp_ser -lxl -lxlfmath -lmassv -lxerces-c

 LDFLAGS = $(LIBPATH) $(LIBS) $(LIBS_MPI)

 PLAT=BGL
 # Blacs libraries
 BLACSDBGLVL   = 0
 BLACSdir      = $(HOME)/software/blacs/bgl/BLACS/LIB
 BLACSCINIT    = $(BLACSdir)/blacsCinit_MPI-$(PLAT)-$(BLACSDBGLVL).a
 BLACSLIB      = $(BLACSdir)/blacs_MPI-$(PLAT)-$(BLACSDBGLVL).a
 CBLACSLIB     = $(BLACSCINIT) $(BLACSLIB) $(BLACSCINIT)

 # Scalapack libraries
 SCALAPACK_DIR = $(HOME)/software/scalapack/bgl/SCALAPACK
 PBLASLIB      = $(SCALAPACK_DIR)/pblas_$(PLAT).a
 SCALAPACKLIB  = $(SCALAPACK_DIR)/libscalapack.a

 # Parallel libraries
 PLIBS = $(SCALAPACKLIB) $(CBLACSLIB)
#-------------------------------------------------------------------------------

