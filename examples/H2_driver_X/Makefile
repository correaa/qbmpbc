# example: make DFLAGS+="-D_MPBC -D_USE_HDF5 -D_SHORTFFT"

include ../../src.mpbc/$(TARGET).mk

TARGET = hydrogen current
INCLUDE += -Wfatal-errors -I../../src.mpbc/driver -I../../src.mpbc -I../../src.mpbc/include -I. 
LIBPATH += -L${HOME}/usr/lib
LIBS    += -lboost_serialization$(BOOST_VER) -lboost_mpi$(BOOST_VER) -lhdf5

all: current 

gedit:
	gedit current.cpp hydrogen.cpp Makefile &


clean:
	rm -f *.o ${TARGET} *~

current: current.cpp ../../src.mpbc/driver/qbox.hpp
	$(LD) $(INCLUDE) $(DFLAGS) $(PLTFLAGS) -o $@ `ls ../../src.mpbc/*.o | grep -v qb.o` $< $(LDFLAGS) 

