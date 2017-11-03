#ifndef BOOST_FILESYSTEM_HPP
#define BOOST_FILESYSTEM_HPP
#include<boost/filesystem.hpp>
namespace boost{
namespace filesystem{
	path home(){
		return std::getenv("HOME");
	}
}
}
#endif
