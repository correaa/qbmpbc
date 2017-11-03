#include<iostream>
#include<boost/array.hpp>
#include<boost/lexical_cast.hpp>

template<typename T, unsigned I>
std::ostream& operator<<(std::ostream& os, boost::array<T,I> const& a){
	os<<"{";
	std::copy(a.begin(), a.end()-1, std::ostream_iterator<T>(os, ","));
	return os<<a.back()<<"}";
}

namespace boost{
template<>
std::string lexical_cast<std::string>(boost::array<long long unsigned, 2u> const& a){
	std::ostringstream ss;
	ss<<"{";
	std::copy(a.begin(), a.end()-1, std::ostream_iterator<int>(ss, ","));
	ss<<a.back()<<"}";
	return ss.str();
}
template<>
std::string lexical_cast<std::string>(boost::array<long long unsigned, 3u> const& a){
	std::ostringstream ss;
	ss<<"{";
	std::copy(a.begin(), a.end()-1, std::ostream_iterator<int>(ss, ","));
	ss<<a.back()<<"}";
	return ss.str();
}
}


