#ifndef GSL_UNITS_ADIMENSIONAL_HPP
#define GSL_UNITS_ADIMENSIONAL_HPP
#include<boost/units/quantity.hpp>
namespace boost{namespace units{

template<class UnitDomain, class UnitResult>
struct adimensionalize_function{
	boost::function<quantity<UnitResult>(quantity<UnitDomain> const&)> f_;
	adimensionalize_function(boost::function<quantity<UnitResult>(quantity<UnitDomain> const&)> f) : f_(f){}
	double operator()(double const& x){return f_(quantity<UnitDomain>::from_value(x)).value();}
};

template<class UnitDomain>
struct adimensionalized_{
	template <typename T, typename F> struct result{ typedef T type; };
	template <typename T,typename F> T operator()(T x, F f) const{
		quantity<UnitDomain,T> forced_temporary = quantity<UnitDomain,T>::from_value(x);
		return f(forced_temporary).value();
	}
//use as:
//	boost::phoenix::function<adimensionalized_<si::length> > const adimensionalized = adimensionalized_<si::length>();
//	boost::function<double(double)> fadd = adimensionalized(arg1, lambda[arg1*arg1]);
};
}}
#endif //GSL_UNITS_ADIMENSIONAL_HPP

