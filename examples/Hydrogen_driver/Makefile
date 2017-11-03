include ../../src.wei/$(TARGET).mk

TARGET = hydrogen current
INCLUDE += -I../../src.wei/driver -I../../src.wei -I../../src.wei/include -I. 
LIBPATH += -L${HOME}/usr/lib
#LIBS    += -lboost_mpi$(BOOST_VER) -lboost_serialization$(BOOST_VER) -lhdf5
LIBS    += -lboost_serialization$(BOOST_VER) -lhdf5

all: current 

gedit:
	gedit current.cpp hydrogen.cpp Makefile &


#hydrogen.o: hydrogen.cpp
#	$(CXX) $(DFLAGS) $(CXXFLAGS) -c $^ 
#	
#current.o: current.cpp 
#	$(CXX) $(DFLAGS) $(CXXFLAGS) -c $^

clean:
	rm -f *.o ${TARGET} *~

hydrogen: hydrogen.cpp ../../src.wei/driver/qbox.hpp
	$(LD) -Wl,-rpath=${HOME}/usr/lib $(INCLUDE) $(DFLAGS) -o $@ `ls ../../src.wei/*.o | grep -v qb.o` $< $(LDFLAGS) -Wfatal-errors

current: current.cpp ../../src.wei/driver/qbox.hpp
	$(LD) -Wl,-rpath=${HOME}/usr/lib $(INCLUDE) $(DFLAGS) -o $@ `ls ../../src.wei/*.o | grep -v qb.o` $< $(LDFLAGS) -Wfatal-errors 

#	g++ -Wl,-rpath=/home/correaa/usr/lib `ls ../../src.wei/*.o | grep -v qb.o` $^ -L/lib -L/home/correaa/usr/lib -L/usr/X11R6/lib -L/usr/lib/mpich/lib -L/home/correaa/usr/local/lapack/lib -L/home/correaa/usr/local/atlas/lib -L/home/correaa/usr/lib -L/home/correaa/usr/lib /usr/lib/libscalapack-mpich.a /usr/lib/libblacsCinit-mpich.a /usr/lib/libblacs-mpich.a /usr/lib/libblacsCinit-mpich.a -lfftw -lfftw3 -llapack -lf77blas -latlas -lm -lboost_filesystem-mt -lboost_system-mt -lboost_mpi -lhdf5 -lxerces-c -Xlinker -Bstatic -lc -lgfortran -static-libgcc -lmpich  -lrt -pthread -Xlinker -Bdynamic -lsicuuc -lsicudata  -o $@
