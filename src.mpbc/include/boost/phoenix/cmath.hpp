#if 0
ln -sf $0 .$0.cpp && c++ .$0.cpp -o .$0.cpp.x -Wall -D_TEST_BOOST_PHOENIX_CMATH_HPP && ./.$0.cpp.x; exit
#endif
#ifndef BOOST_PHOENIX_CMATH_HPP
#define BOOST_PHOENIX_CMATH_HPP
#include<boost/spirit/home/phoenix.hpp>
#include <boost/utility/enable_if.hpp>


/// Extends any endomorphic function as a phoenix function, endomorphic in the sense that the input and output are of the same type
#define BOOST_PHOENIX_REGISTER_ENDOMORPHISM(FN) \
namespace boost{ \
namespace phoenix{ \
	struct FN##_eval{ \
		template <typename Env, typename Arg_> struct result{typedef typename boost::mpl::at_c<typename Env::args_type, 0>::type type;}; \
		template <typename RT, typename Env, typename Arg_> static RT eval(Env const& env, Arg_ const& arg_){ \
			return FN(arg_.eval(env)); \
		} \
	}; \
	template <typename Arg_> \
	typename boost::enable_if_c<is_actor<Arg_>::value, \
		actor<typename as_composite<FN##_eval, Arg_>::type> \
	>::type \
	FN(Arg_ const& arg_){ \
		return compose<FN##_eval>(arg_); \
	} \
}}

#include<cmath>
BOOST_PHOENIX_REGISTER_ENDOMORPHISM(cos);
BOOST_PHOENIX_REGISTER_ENDOMORPHISM(cosh);
BOOST_PHOENIX_REGISTER_ENDOMORPHISM(exp);
BOOST_PHOENIX_REGISTER_ENDOMORPHISM(log);
BOOST_PHOENIX_REGISTER_ENDOMORPHISM(log10);
BOOST_PHOENIX_REGISTER_ENDOMORPHISM(sin);
BOOST_PHOENIX_REGISTER_ENDOMORPHISM(sinh);
BOOST_PHOENIX_REGISTER_ENDOMORPHISM(sqrt);
BOOST_PHOENIX_REGISTER_ENDOMORPHISM(tan);
BOOST_PHOENIX_REGISTER_ENDOMORPHISM(tanh);

double sech(double x){return 1./cosh(x);}
BOOST_PHOENIX_REGISTER_ENDOMORPHISM(sech);

#endif

#ifdef _TEST_BOOST_PHOENIX_CMATH_HPP
#include<boost/function.hpp>
#include<boost/lambda/lambda.hpp>
#include<iostream>

using std::cout; using std::endl;
int main(){
	using namespace boost::phoenix::arg_names;
	//boost::function<double(std::complex<double>)> f = real(arg1);
	double dos = 2.;
	cout << sech(arg1)(dos) << " "<< sech(2.) <<endl;
//	boost::function<double(std::complex<double>)> f = abs(arg1 + 1. + imag(arg1) + sin(arg1) );
//	std::complex<double> arg(1.,2.);
//	cout << f(std::complex<double>(1.,2.))<<endl;
	return 0;
}
#endif

