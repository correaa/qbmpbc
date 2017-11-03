#ifdef compile_instructions
ln -sf $0 .$0.cpp && c++ -std=c++0x -I$HOME/usr/include -I$HOME/prj .$0.cpp -Wall -Wfatal-errors -L$HOME/lib `pkg-config --libs gsl` -D_TEST_GSL_FIT_LINEAR_BOOST_UNITS_HPP -o ./.$0.x && ./.$0.x $@
rm -f .$0.x .$0.cpp
exit
#endif
#ifndef GSL_FIT_LINEAR_BOOST_UNITS_HPP
#define GSL_FIT_LINEAR_BOOST_UNITS_HPP
#include "../../gsl/fit_linear.hpp"
#include<boost/units/quantity.hpp>
namespace gsl{
namespace linear{
namespace units{
	using namespace boost::units;
	template<class UX, class UY>
	struct model : protected gsl::linear::model{
		model(
			quantity<UX> const& c0,
			quantity<typename divide_typeof_helper<UY, UX>::type> const& c1
		) : 
			gsl::linear::model(c0.value(), c1.value()){
		}
		quantity<UY> operator()(quantity<UX> const& x) const{
			return quantity<UY>::from_value(gsl::linear::model::operator()(x.value()));
		}
		quantity<UY> const& c0(){return reinterpret_cast<quantity<UY> const&>(gsl::linear::model::c0());}
		quantity<typename divide_typeof_helper<UY, UX>::type> const& c1() const{return reinterpret_cast<quantity<typename divide_typeof_helper<UY, UX>::type> const&>(gsl::linear::model::c1());}
	};
	template<class UX, class UY>
	model<UX, UY> fit(std::vector<std::pair<quantity<UX>, quantity<UY> > > const& v){
		gsl::linear::model m = gsl::linear::fit(reinterpret_cast<std::vector<std::pair<double, double> > const&>(v));
		return model<UX, UY>(
			quantity<UX>::from_value(m.c0()),
			quantity<typename divide_typeof_helper<UY, UX>::type>::from_value(m.c1())
		);
	}
	linear::model linear(std::vector<std::pair<double,double> > const& vp){
		double c0, c1;{
			double cov00, cov01, cov11, sumsq;
			gsl_fit_linear(
				&vp[0].first, /*const size_t xstride*/ 2, &vp[0].second, /*const size_t xstride*/ 2, vp.size(), 
					&c0   , &c1, 
					&cov00, /**/
					&cov01, &cov11, 
					&sumsq
			);
		}
		return linear::model(c0, c1);
	}
}
}
}

#endif
#ifdef _TEST_GSL_FIT_LINEAR_BOOST_UNITS_HPP
#include<boost/units/systems/si.hpp>
#include<boost/units/systems/si/codata_constants.hpp>
#include<boost/units/systems/si/io.hpp>
using namespace boost::units;
int main(){
	std::vector<std::pair<quantity<si::temperature>, quantity<si::energy> > > data(30); //do not use push back
	for(unsigned i = 0; i < data.size(); i++) {
		quantity<si::temperature> t = (i*2000.)*si::kelvin;
		data[i] = std::make_pair(
			t,
			(3.1 + (i%7-5)/100.) * si::constants::codata::k_B * t
		);
	};
	clog << data[3].first << " -- " << data[3].second;
	gsl::linear::units::model<si::temperature, si::energy> l(
		10.*si::kelvin,
		2.9*si::constants::codata::k_B
	);
	//gsl::linear::units::model<si::temperature, si::energy> m;
	auto m = gsl::linear::units::fit(data);
	clog << endl <<  m(1000.*si::kelvin) << endl;
	return 0;
}
#endif
