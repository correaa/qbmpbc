#ifdef compile_instructions 
	 ln -sf $0 .$0.cpp; c++ -D_GSL_EPSILON_HPP_TEST .$0.cpp -o .$0.x  -Wall && ./.$0.x ; exit
#endif
#ifndef GSL_EPSILON_HPP
#define GSL_EPSILON_HPP
#include<utility> //std::pair
#include<limits>
namespace gsl{ //or gsl/detail
class epsilon : std::pair<double, double>{
	public:
	epsilon(
		double absolute_ = std::numeric_limits<double>::epsilon(), 
		double relative_ = std::numeric_limits<double>::epsilon()
	) : std::pair<double,double>(absolute_, relative_){}
	double const& absolute() const{return first;}
	double const& relative() const{return second;}
	static epsilon relative(double relative_){return epsilon(0., relative_);}
	static epsilon absolute(double absolute_){return epsilon(absolute_, 0.);}
};
}
#endif

#ifdef _GSL_EPSILON_HPP_TEST
#include<iostream>
int main(){
	std::cout << gsl::epsilon().absolute() <<" "<< gsl::epsilon().relative() << std::endl;
	return 0;
}
#endif

