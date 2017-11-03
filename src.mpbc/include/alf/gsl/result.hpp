#ifdef compile_instructions 
	 ln -sf $0 .$0.cpp; c++ -D_GSL_RESULT_HPP_TEST .$0.cpp -o .$0.x  -Wall && ./.$0.x ; exit
#endif
#ifndef GSL_RESULT_HPP_
#define GSL_RESULT_HPP_
#include<map> //for pair
#include<boost/numeric/interval.hpp>
#include<iostream>
#include<boost/array.hpp>
namespace gsl{
struct result : std::pair<double,double>{
	result(double est, double err=0) : std::pair<double, double>(est, err){}
	double estimate() const{return first;}
	double error() const{return second;}
	operator double() const{return estimate();}
	//operator boost::numeric::interval<double>() const{return boost::numeric::interval<double>(this->estimate()-error(),estimate()+error());} //make it operator?
	result(std::pair<double, double> const& p) : std::pair<double, double>(p){}
	friend std::ostream& operator<<(std::ostream& os, result const& r){return os<<r.estimate()<<" \u00B1 "<<r.error();} // \ub00B1 is unicode +/-
	private:
	//result(){} //had to add it for boost::array initialization
	template<typename T, size_t s> friend class boost::array;
};
}
#endif

#ifdef _GSL_RESULT_HPP_TEST
#include<iostream>
int main(){
	std::cout << gsl::result(2., 0.34) <<" " << std::endl;
	return 0;
}
#endif

