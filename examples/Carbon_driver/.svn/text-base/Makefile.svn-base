include ../../src.wei/$(TARGET).mk

TARGET = current
INCLUDE += -I../../src.wei/driver -I../../src.wei -I../../src.wei/include -I. 
LIBPATH += -L${HOME}/usr/lib
LIBS    += -lboost_mpi$(BOOST_VER) -lhdf5

all: current

gedit:
	gedit current.cpp hydrogen.cpp Makefile &

clean:
	rm -f *.o ${TARGET} *~

current: current.cpp ../../src.wei/driver/qbox.hpp
	$(LD) $(LDFLAGS) $(INCLUDE) $(DFLAGS) $(PLTFLAGS) -o $@ `ls ../../src.wei/*.o | grep -v qb.o` $< $(LDFLAGS) -Wfatal-errors 

