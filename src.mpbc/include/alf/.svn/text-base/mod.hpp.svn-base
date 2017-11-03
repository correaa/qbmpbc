#if 0
        cp $0 $0.cpp
        c++ $0.cpp -o $0.cpp.x -I$HOME/usr/include -D_TEST_MOD_HPP && ./$0.cpp.x
        rm -f $0.cpp
        exit
#endif
#ifndef MOD_HPP_
#define MOD_HPP_
#include <boost/numeric/conversion/cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <cmath>
namespace alf{
	unsigned mod(int i, unsigned d){
		return ((i%(int)d)+(int)d)%(int)d;
	}
	int floor(double const& d){
		return boost::numeric_cast<int>(::floor(d));
	}
	double frac(double const& d){
		return d-floor(d);
	}
	boost::tuple<int, double> modf(double const& d){
		boost::tuple<int, double> ret;
		{double integer_part_;
			ret.get<1>() = std::modf(d, &integer_part_);
			ret.get<0>() = boost::numeric_cast<int>(integer_part_);
		}
		if(d<0){
		  --ret.get<0>();
		  ret.get<1>()+=1.;
		}
		return ret;
	}
}
#endif
#ifdef _TEST_MOD_HPP
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
int main(){
  std::cout<< "std::floor(5.4) " << std::floor(5.4) << " " << typeid(std::floor(5.4)).name()<< std::endl;
  std::cout<< "alf::floor(5.4) " << alf::floor(5.4) << " " << typeid(alf::floor(5.4)).name()<< std::endl;
  std::cout<< "std::floor(-5.4) " << std::floor(-5.4) << " " << typeid(std::floor(-5.4)).name()<< std::endl;
  std::cout<< "alf::floor(-5.4) " << alf::floor(-5.4) << " " << typeid(alf::floor(-5.4)).name()<< std::endl;
  std::cout<< "alf::modf(-5.4) " << alf::modf(-5.4)<<" " <<typeid(alf::modf(-5.4)).name()<<", " << std::endl;
  return 0;
}
#endif

