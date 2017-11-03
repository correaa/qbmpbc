#if 0
ln -sf $0 .$0.cpp && c++ -std=c++0x -I$HOME/prj -I$HOME/usr/include .$0.cpp -Wall -Wfatal-errors -L$HOME/lib `pkg-config --libs gsl` -lboost_system -lboost_regex -lhunspell -D_TEST_INTERPOLATION_UNITS_HPP -o ./.$0.x && ./.$0.x $@; exit
#endif
#ifndef INTERPOLATION_UNITS_HPP
#define INTERPOLATION_UNITS_HPP
#include "../../gsl/interpolation.hpp"
#include<boost/units/quantity.hpp>
#include"boost/units/interval.hpp" //compat boost::units/boost::interval
namespace gsl{
namespace interpolation{
namespace units{
	using namespace boost::units;
	template<class U1, class U2>
	struct spline{
		gsl::interpolation::spline impl;
		spline(std::map<quantity<U1>, quantity<U2> > const& m)
		: impl(reinterpret_cast<std::map<double, double> const&>(m)){}
		template <typename Arg> struct result{
			typedef quantity<U2> type;
		};
		quantity<U2> operator()(quantity<U1> const& x) const{
			return quantity<U2>::from_value(impl(x.value()));
		}
	};
}}}
#endif //INTERPOLATION_UNITS_HPP

#ifdef _TEST_INTERPOLATION_UNITS_HPP
#include<iostream>
#include<boost/units/systems/si.hpp>
using namespace std;
using namespace boost::units;
int main(){
	std::map<quantity<si::length>, quantity<si::energy> > m;
	m[1.*si::meter] = 1.0*si::joules;
	m[2.*si::meter] = 2.1*si::joules;
	m[3.*si::meter] = 3.2*si::joules;
	m[4.*si::meter] = 4.2*si::joules;
	m[5.*si::meter] = 5.2*si::joules;

	//std::map<double, double> const& mdd = reinterpret_cast<std::map<double, double> const&>(m);
	//std::clog << mdd[2] << std::endl;
	//std::clog << m[2.*si::meter] << std::endl;
	gsl::interpolation::units::spline<si::length, si::energy> s(m);
	std::clog << s(4.1*si::meter) << std::endl;
	return 0;
}
#endif
/* vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent: */

