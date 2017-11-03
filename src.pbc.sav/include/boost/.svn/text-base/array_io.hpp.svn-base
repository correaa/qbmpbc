#include<iostream>
#include<boost/array.hpp>
template<typename T, unsigned I>
std::ostream& operator<<(std::ostream& os, boost::array<T,I> const& a){
	os<<"{";
	std::copy(a.begin(), a.end()-1, std::ostream_iterator<int>(os, ","));
	return os<<a.back()<<"}";
}
