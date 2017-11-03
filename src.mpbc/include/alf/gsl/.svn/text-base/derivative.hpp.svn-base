#ifdef compile_instructions
ln -sf $0 .$0.cpp && c++ -std=c++0x .$0.cpp -Wfatal-errors -L$HOME/lib `pkg-config --libs gsl` -D_TEST_DERIVATIVE_HPP -o ./.$0.x && ./.$0.x $@
rm -f .$0.x .$0.cpp
exit
#endif

#ifndef GSL_DERIVATIVE_HPP
#define GSL_DERIVATIVE_HPP

#include <gsl/gsl_deriv.h>
#include "../gsl/function.hpp"
#include "../gsl/result.hpp"
#include<boost/spirit/home/phoenix.hpp>

namespace boost { 
	namespace phoenix{
    template <typename F>
    struct function_crtp{
        actor<typename as_composite<detail::function_eval<0>, F>::type>
        operator()() const
        {
            return compose<detail::function_eval<0> >((F const&)*this);
        }
        template <typename A0>
        actor<typename as_composite<detail::function_eval<1>, F, A0>::type>
        operator()(A0 const& _0) const
        {
            return compose<detail::function_eval<1> >((F const&)*this, _0);
        }

        template <typename A0, typename A1>
        actor<typename as_composite<detail::function_eval<2>, F, A0, A1>::type>
        operator()(A0 const& _0, A1 const& _1) const
        {
            return compose<detail::function_eval<2> >((F const&)*this, _0, _1);
        }
        //  Bring in the rest of the function call operators
        //#include <boost/spirit/home/phoenix/function/detail/function_call.hpp>
		};
	}
}



namespace gsl{
namespace derivative{ //based on 
gsl::result central(gsl::function const& f, double const& x, double const& suggested_step_h = 1e-8){
	double res;
	double abserr;
	gsl_deriv_central(&f, x, suggested_step_h , &res, &abserr);
	return gsl::result(res, abserr);
}
gsl::result forward(gsl::function const& f, double const& x, double const& suggested_step_h = 1e-8){
	double res, abserr;
	gsl_deriv_forward(&f, x, suggested_step_h , &res, &abserr);
	return gsl::result(res, abserr);
}


template<class Functor>
struct central_impl{
	mutable Functor f_;
	central_impl(Functor const& f): f_(f){}
	template<class Args> struct result{typedef double type;};
	double operator()(double const& x) const{
		return central(f_, x);
	}
	template <typename A0> //adapt this functor as phoenix function
	boost::phoenix::actor<typename boost::phoenix::as_composite<boost::phoenix::detail::function_eval<1>, central_impl<Functor>, A0>::type>
	operator()(A0 const& _0) const{
		return boost::phoenix::compose<boost::phoenix::detail::function_eval<1> >(*this, _0);
	}
};

template<class Functor>
central_impl<Functor> central(Functor const& f){
	return central_impl<Functor>(f);
}

template<double(f)(double const&)>
gsl::result central(double const& x){
	return central(f, x);
}
template<double(f)(double)>
gsl::result central(double const& x){
	return central(f, x);
}
template<gsl::result(f)(double const&)>
gsl::result central(double const& x){
	return central(f, x);
}

template<unsigned Degree>
struct central_helper{
	template<double(f)(double const&)>
	static double eval(double const& x);
};
template<>
template<double(f)(double const&)>
double central_helper<0>::eval(double const& x){
	return f(x);
}
template<>
template<double(f)(double const&)>
double central_helper<1>::eval(double const& x){
	return central(f, x);
}
template<>
template<double(f)(double const&)>
double central_helper<2>::eval(double const& x){
	return central(central<f>, x);
}

}}
#endif

#ifdef _TEST_DERIVATIVE_HPP
#include<iostream>
#include<boost/units/systems/si.hpp>


using std::cout;
using std::clog;
using std::endl;
double f (double const& x){
	return pow (x, 1.5);
}
double sine(double const& x){
	return sin(x);
}

using namespace boost::units;
quantity<si::area> sq(quantity<si::length> const& l){
	return l*l;
}

struct f_impl{
	double operator()(double const& x) const{
		return x*x-x;
	}
};


int main (void){
	gsl::function F(&f);
	cout << gsl::derivative::central(F, 2.0) << endl;
	cout << gsl::derivative::central<f>(2.0) << endl;
	cout << "exact = "<< 1.5 * sqrt(2.0) <<endl;
	cout << gsl::derivative::central_helper<0>::eval<sine>(2.0) << endl;
	cout << gsl::derivative::central_helper<1>::eval<sine>(2.0) << endl;
	cout << gsl::derivative::central_helper<2>::eval<sine>(2.0) << endl;
	cout << gsl::derivative::central<gsl::derivative::central<sine> >(2.0) << endl;
	f_impl f;
	auto df = gsl::derivative::central(f);
	//gsl::derivative::central_<f_impl> df(f);
	//cout << df(4.) << endl;
	using namespace boost::phoenix::arg_names;
	//df(arg1+1.);
	double const three = 3.;
	cout << df(arg1+1.)(three) << endl;
	return 0;
}

#endif
			
