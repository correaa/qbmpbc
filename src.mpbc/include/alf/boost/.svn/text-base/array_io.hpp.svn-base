#ifndef ARRAY_IO_HPP_
#define ARRAY_IO_HPP_
#include<iostream>
#include<boost/array.hpp>
//#include<boost/lexical_cast.hpp>

/*
template<typename T, unsigned I>
std::ostream& operator<<(std::ostream& os, boost::array<T,I> const& a){
	os<<"{";
	std::copy(a.begin(), a.end()-1, std::ostream_iterator<T>(os, ","));
	return os<<a.back()<<"}";
}*/

//emplate<typename T, unsigned I>
//std::ostream& operator<<(std::ostream& os, boost::array<unsigned,3u> const& a){
//	os<<"{";
//	std::copy(a.begin(), a.end()-1, std::ostream_iterator<unsigned>(os, ","));
//	return os<<a.back()<<"}";
//}

namespace boost{
	template<typename T, unsigned N>
	std::string
	to_string(array<T,N> const& a){  // I had to write this function because I couldn't make boost::lexical_cast<std::string>(boost::array<T,N> const&) to work
		std::ostringstream os;
		os<<"[";
		std::copy(a.begin(), a.end()-1, std::ostream_iterator<unsigned>(os, ","));
		os<<a.back()<<"]";
		return os.str();
	}
/*
template<>
std::string lexical_cast<std::string, boost::array<unsigned, 2u> >(boost::array<unsigned, 2u> const& a){
	std::ostringstream ss;
	ss<<"{";
	std::copy(a.begin(), a.end()-1, std::ostream_iterator<unsigned>(ss, ","));
	ss<<a.back()<<"}";
	return ss.str();
}
//template<>
std::string lexical_cast<std::string, boost::array<unsigned, 3u> >(boost::array<unsigned, 3u> const& a){
	std::ostringstream ss;
	ss<<"{";
	std::copy(a.begin(), a.end()-1, std::ostream_iterator<unsigned>(ss, ","));
	ss<<a.back()<<"}";
	return ss.str();
}
*/
}
#endif

