#ifndef _YOGRT_HPP_
#define _YOGRT_HPP_
#include<yogrt.h>
#include<boost/mpi.hpp>
namespace yogrt {
int /*seconds*/ remaining() {
	return ::yogrt_remaining();
}
namespace mpi {
int /*seconds*/ remaining(boost::mpi::communicator const& comm) {
	int ret = -1;
	if(comm.rank() == 0) {
		ret =::yogrt_remaining();
	}
	boost::mpi::broadcast(comm, ret, 0);
	assert(ret != -1);
	return ret;
}
}
}
// use as:
// in mpi process (synchronized value): 
// if(yogrt::mpi::remaining(world) break;
#endif
