#if COMPILATION_INSTRUCTIONS
	ln -sf $0 $0.cpp &&
	/opt/mvapich-gnu-gen2-0.9.9/bin/mpicxx $0.cpp -o $0.x -Wfatal-errors -Wall \
		-L$HOME/usr/lib \
		-lboost_mpi -lboost_serialization -lboost_system \
		-Wl,-rpath=/opt/mvapich-gnu-gen2-0.9.9/lib/shared \
		-D_MPI_OSTREAM_HPP_TEST \
	&& mpirun -np 2 ./$0.x $1 $2 $3 $4 $5 $6 $7 $8 $9;
	exit
#endif
#ifndef MPI_OSTREAM_HPP
#define MPI_OSTREAM_HPP
#include<boost/mpi.hpp>
#include<set>
namespace boost{
namespace mpi{
class ostream /*: public std::ostringstream*/{
protected:
	std::ostringstream /*&*/ local_buffer_;
	communicator const& comm_;
	std::ostream& os_;
	ostream(ostream const&);
	public:
	ostream(communicator const& c=communicator(), std::ostream& o=std::clog) : /*local_buffer_((std::ostringstream&)(*this)),*/ comm_(c), os_(o){}
	template<typename T> ostream& operator<<(T const& t){
		local_buffer_<<t;
		return *this;
	}
       typedef ostream& (*ostream_manipulator)(ostream&);
       ostream& operator<<(ostream_manipulator manip){
         return manip(*this);
       }
       ~ostream(){flush(*this);}
	static ostream& flush(ostream& os){
		os.local_buffer_<<std::flush;
		os.comm_.barrier();
		os.comm_.send(0, os.comm_.rank(), os.local_buffer_.str());
		os.local_buffer_.str("");
		if(os.comm_.rank()==0){
			std::set<std::string> messages;
			for(int i=0; i!=os.comm_.size(); ++i){
				std::string remote_buffer_;
				os.comm_.recv(i, i, remote_buffer_);
				messages.insert(remote_buffer_); 
			}
			for(std::set<std::string>::const_iterator i=messages.begin(); i!=messages.end(); ++i){
				os.os_<<(*i);//<<std::endl;
			}
			os.os_<<std::flush;
		}
		return os;

	}
	static ostream& endl(ostream& os){
		os.local_buffer_<<'\n'; 
		return flush(os);
	}
		typedef std::ostream& (*std_manipulator)(std::ostream&);
		ostream& operator<<(std_manipulator manip){
			if(manip==(std::ostream& (*)(std::ostream&))std::endl){//need cast due to some overload
				return ostream::endl(*this); //this can avoid lots of confusion
			}
			if(manip==(std::ostream& (*)(std::ostream&))std::flush){
				return ostream::flush(*this);
			}
			manip(local_buffer_);
			return *this;
		}
	};
	ostream& endl(ostream& os){
		return ostream::endl(os);
	}
	ostream& flush(ostream& os){
		return ostream::flush(os);
	}
}}
#include<boost/filesystem/fstream.hpp>
namespace boost{
namespace mpi{
	class ofstream : public boost::mpi::ostream{
		boost::filesystem::ofstream ofs_;
		public:
		ofstream(communicator const& c, boost::filesystem::path const& p) : boost::mpi::ostream(c, ofs_){
			if(c.rank()==0) ofs_.open(p);
		}
	};
	namespace{
		ostream clog(communicator const& c, std::ostream& o);
	}
}}
#ifdef _MPI_OSTREAM_HPP_TEST
int main(int argc, char* argv[]){
   boost::mpi::environment env(argc, argv);
   boost::mpi::communicator world;	

   using boost::mpi::endl; //or using std::endl; //*
	//std::endl is also parallel aware to avoid confusion, good idea?
   boost::mpi::ostream clog(world); // or using std::clog; //*

   std::clog<<"someone says: are we the world? "<<world.rank()<<std::endl; //not parallel aware
   clog<<"world says: we are the world"<<endl;          //this messages is unique
   clog<<"rank "<<world.rank()<<" says: we are the children"<<endl; //this message is different for each process.
 
   boost::mpi::ofstream mpi_ofs(world, "ofstream.dat");
	mpi_ofs<<"hola "<<world.rank()<<std::endl;
   return 0;
}
#endif
#endif

