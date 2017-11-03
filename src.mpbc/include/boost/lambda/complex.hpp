#if 0
ln -sf $0 $0.cpp && c++ $0.cpp -o .$0.cpp.x -Wall -D_TEST_BOOST_LAMBDA_COMPLEX_HPP && ./.$0.cpp.x; exit
#endif
#ifndef BOOST_LAMBDA_COMPLEX_HPP
#define BOOST_LAMBDA_COMPLEX_HPP

#include "boost/lambda/core.hpp"
#include<boost/lambda/bind.hpp>
#include<complex>
//#include <numeric>

namespace boost{
namespace lambda{

#define BOOST_LAMBDA_COMPLEX_IMPL(function_name) \
template<typename Arg> \
lambda_functor<lambda_functor_base< \
	action<2, function_action<2> >, \
	boost::tuples::tuple< \
		const double& (* const)(const std::complex<double>&), \
		const lambda_functor<Arg> \
	> \
> > function_name(const lambda_functor<Arg>& f){ \
	return boost::lambda::bind(static_cast<double const&(*)(std::complex<double> const&)>(&std::function_name<double>), f); \
}

BOOST_LAMBDA_COMPLEX_IMPL(real)
BOOST_LAMBDA_COMPLEX_IMPL(imag)
BOOST_LAMBDA_COMPLEX_IMPL(abs)
BOOST_LAMBDA_COMPLEX_IMPL(norm)

}
}
#endif

#ifdef _TEST_BOOST_LAMBDA_COMPLEX_HPP
#include<boost/function.hpp>

#include<boost/lambda/lambda.hpp>
#include<iostream>
using std::cout; using std::endl;
int main(){
	boost::function<double(double)> f = boost::lambda::bind(static_cast<double(*)(double)>(&::std::sin), boost::lambda::_1+1.);
	cout << f(1.)<<endl;
	boost::function<double(std::complex<double>)> fi = boost::lambda::bind(static_cast<double const&(*)(std::complex<double> const&)>(&std::imag<double>), boost::lambda::_1+ std::complex<double>(1.,1.));
	cout << fi(std::complex<double>(1.,2.))<<endl;
	
	boost::function<double(std::complex<double>)> f2 = imag(boost::lambda::_1)+real(boost::lambda::_1);
	cout << f2(std::complex<double>(1.,2.))<<endl;
	return 0;
}
#endif

