#ifdef compilation_instructions
ln -sf $0 $0.cpp && c++ $0.cpp -o $0.cpp.x -lboost_mpi -lboost_filesystem -lboost_system -Wfatal-errors -D_UNIQUE_INDEXED_PATH_TEST && ./$0.cpp.x; exit
#endif
#include<boost/lexical_cast.hpp>
#include<boost/filesystem.hpp>
#include"../format_cast.hpp"
#include<boost/filesystem/fstream.hpp>
//#include<boost/wave/util/cpp_include_paths.hpp>
namespace boost{
namespace filesystem{
	path unique_indexed_path(std::string const& s){
		path ret;
		for(unsigned idx=0; ; ++idx){
			std::ostringstream oss;
			oss<<boost::format(s)%idx;
			if(not exists(ret = path(oss.str()))) break;
		}
		return ret;
	}
}}

#include<boost/mpi.hpp>
namespace boost{
namespace filesystem{
	namespace mpi{
		path unique_indexed_path(std::string const& s, boost::mpi::communicator& comm){
			std::string ret;
			if(comm.rank()==0) ret=boost::filesystem::unique_indexed_path(s).string();
			boost::mpi::broadcast(comm, ret, 0);
			assert(not ret.empty());
			return ret;
		}
	}
}}


#ifdef _UNIQUE_INDEXED_PATH_TEST
int main(){
	boost::filesystem::path p = boost::filesystem::unique_indexed_path("unique-%03u.dat");
	boost::filesystem::ofstream ofs(p);
	ofs<<p.string()<<std::endl;
	return 0;
}
#endif
