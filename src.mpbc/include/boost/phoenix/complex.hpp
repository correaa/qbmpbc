#if 0
ln -sf $0 .$0.cpp && c++ .$0.cpp -o .$0.cpp.x -Wall -D_TEST_BOOST_PHOENIX_COMPLEX_HPP && ./.$0.cpp.x; exit
#endif
#ifndef BOOST_PHOENIX_COMPLEX_HPP
#define BOOST_PHOENIX_COMPLEX_HPP
#include<boost/spirit/home/phoenix.hpp>
#include<complex>
namespace boost{
namespace phoenix{
#if 0
//first try, but is either real_ or std::real
struct real_impl{
    template <typename Arg> struct result{
        typedef double type;
    };
    template <typename Arg> double
    operator()(Arg z) const{
        return std::real(z);
    }
};
boost::phoenix::function<real_impl> real;
#endif

#define BOOST_PHOENIX_UNARY_COMPLEX_TO_REAL( FunctionnamE ) \
	struct FunctionnamE##_eval{ \
		template <typename Env, typename Arg_> struct result{ typedef typename boost::mpl::at_c<typename Env::args_type, 0>::type::value_type type;}; \
		template <typename RT, typename Env, typename Arg_> static RT eval(Env const& env, Arg_ const& arg_){ \
			return FunctionnamE(arg_.eval(env)); \
		} \
	}; \
	template <typename Arg_> actor<typename as_composite<FunctionnamE##_eval, Arg_>::type> \
	FunctionnamE(Arg_ const& arg_){ \
		return compose<FunctionnamE##_eval>(arg_); \
	}
//as listed in http://www.cplusplus.com/reference/std/complex/
BOOST_PHOENIX_UNARY_COMPLEX_TO_REAL( real )
BOOST_PHOENIX_UNARY_COMPLEX_TO_REAL( imag )
BOOST_PHOENIX_UNARY_COMPLEX_TO_REAL( abs )
BOOST_PHOENIX_UNARY_COMPLEX_TO_REAL( norm )
BOOST_PHOENIX_UNARY_COMPLEX_TO_REAL( arg ) //argument complex function, do not confuse with arg_
#undef BOOST_PHOENIX_UNARY_COMPLEX_TO_REAL

#define BOOST_PHOENIX_UNARY_COMPLEX_TO_COMPLEX( FunctionnamE ) \
	struct FunctionnamE##_eval{ \
		template <typename Env, typename Arg_> struct result{typedef typename boost::mpl::at_c<typename Env::args_type, 0>::type type;}; \
		template <typename RT, typename Env, typename Arg_> static RT eval(Env const& env, Arg_ const& arg_){ \
			return FunctionnamE(arg_.eval(env)); \
		} \
	}; \
	template <typename Arg_> actor<typename as_composite<FunctionnamE##_eval, Arg_>::type> \
	FunctionnamE(Arg_ const& arg_){ \
		return compose<FunctionnamE##_eval>(arg_); \
	}
	BOOST_PHOENIX_UNARY_COMPLEX_TO_COMPLEX( conj  )
	//todo: polar 

	BOOST_PHOENIX_UNARY_COMPLEX_TO_COMPLEX( cos   )
	BOOST_PHOENIX_UNARY_COMPLEX_TO_COMPLEX( cosh  )
	BOOST_PHOENIX_UNARY_COMPLEX_TO_COMPLEX( exp   )
	BOOST_PHOENIX_UNARY_COMPLEX_TO_COMPLEX( log   )
	BOOST_PHOENIX_UNARY_COMPLEX_TO_COMPLEX( log10 )
//	BOOST_PHOENIX_UNARY_COMPLEX_TO_COMPLEX( pow   ) will be replaced by pow<>?, it is more elegant!
	BOOST_PHOENIX_UNARY_COMPLEX_TO_COMPLEX( sin   )
	BOOST_PHOENIX_UNARY_COMPLEX_TO_COMPLEX( sinh  )
	BOOST_PHOENIX_UNARY_COMPLEX_TO_COMPLEX( sqrt  )
	BOOST_PHOENIX_UNARY_COMPLEX_TO_COMPLEX( tan   )
	BOOST_PHOENIX_UNARY_COMPLEX_TO_COMPLEX( tanh  )
#undef BOOST_PHOENIX_UNARY_COMPLEX_TO_COMPLEX
}}
#endif

#ifdef _TEST_BOOST_PHOENIX_COMPLEX_HPP
#include<boost/function.hpp>
#include<boost/lambda/lambda.hpp>
#include<iostream>
using std::cout; using std::endl;
int main(){
	using namespace boost::phoenix::arg_names;
	//boost::function<double(std::complex<double>)> f = real(arg1);
	std::complex<double> z(1.,2.);
	cout << real(arg1)(z) <<endl;
//	boost::function<double(std::complex<double>)> f = abs(arg1 + 1. + imag(arg1) + sin(arg1) );
//	std::complex<double> arg(1.,2.);
//	cout << f(std::complex<double>(1.,2.))<<endl;
	boost::function<double(std::complex<double>)> f = imag(arg1);
	cout << f(std::complex<double>(3.,5.)) << endl;
	return 0;
}
#endif

