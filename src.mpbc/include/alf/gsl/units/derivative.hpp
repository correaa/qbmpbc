#ifdef compile_instructions
ln -sf $0 .$0.cpp && c++ -std=c++0x -I$HOME/usr/include -I$HOME/prj .$0.cpp -Wall `#-Wfatal-errors` -L$HOME/lib `pkg-config --libs gsl` -D_TEST_GSL_UNITS_DERIVATIVE_HPP -o ./.$0.x && ./.$0.x $@
rm -f .$0.x .$0.cpp; exit
#endif
#ifndef UNITS_DERIVATIVE_HPP
#define UNITS_DERIVATIVE_HPP
#include<boost/units/quantity.hpp>
#include "../../gsl/derivative.hpp"
#include "../units/adimensional.hpp"

namespace gsl{
namespace derivative{
namespace units{
using namespace boost::units;
template<class UnitDomain, class UnitFunction>
typename divide_typeof_helper<quantity<UnitFunction>, quantity<UnitDomain> >::type aa(
	boost::function<quantity<UnitFunction>(quantity<UnitDomain>)> const& f,
	quantity<UnitDomain> x
){
	boost::units::adimensionalize_function<UnitDomain, UnitFunction> ad(f);
	return divide_typeof_helper<quantity<UnitFunction>, quantity<UnitDomain> >::type::from_value(gsl::derivative::central(ad, x.value()));
}
}}}
#endif

#ifdef _TEST_GSL_UNITS_DERIVATIVE_HPP
#include<boost/units/systems/si.hpp>
#include<iostream>

using std::cout; using std::endl;
using namespace boost::units;

quantity<si::area> square_area(quantity<si::length> const& l){
	return l*l;
}

int main(){
	cout << square_area(2.*si::meter) << endl;
	boost::function<quantity<si::area>(quantity<si::length>)> f = &square_area;
	cout << gsl::derivative::units::aa< si::length, si::area >(f, 1.*si::meter) << endl;

	return 0;
}
#endif

// Editor modelines  -  http://www.wireshark.org/tools/modelines.html
// Local variables:
// c-basic-offset: 4
// tab-width: 4
// indent-tabs-mode: t
// truncate-lines: 1
// End:
// vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent:

