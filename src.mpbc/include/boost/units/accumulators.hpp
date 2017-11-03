#ifdef COMPILE_INSTRUCTIONS
ln -sf $0 .$0.cpp && c++ -Wfatal-errors `#-std=c++0x` .$0.cpp -o .$0.x -D_BOOST_UNITS_ACCUMULATORS_TEST -I$HOME/prj && ./.$0.x $@; exit;
#endif
#ifndef BOOST_UNITS_ACCUMULATORS
#define BOOST_UNITS_ACCUMULATORS
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/variance.hpp>
#include <boost/units/unit.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/cmath.hpp> //for power_typeof_helper
namespace boost {namespace numeric{namespace functional{
	using namespace boost::units;
	template<class Unit, typename Y> struct quantity_tag{}; // struct mydouble_tag{};
	template<class Unit, typename Y> struct tag<quantity<Unit,Y> >{ //template<> struct tag<mydouble>{
		typedef quantity_tag<Unit,Y> type;
	};
	// Specify how to divide an object by an integral count (right)
	template<typename Left, typename Right>
	struct average<Left, Right, quantity_tag<typename Left::unit_type, typename Left::value_type>, void>{
		typedef Left result_type;
		result_type operator()(Left & left, Right & right) const{
		    return left/(typename Left::value_type)right;
		}
	};
	template<typename Left, typename Right>
	struct multiplies<Left, Right, quantity_tag<typename Left::unit_type, typename Left::value_type>, quantity_tag<typename Left::unit_type, typename Left::value_type> >{
		typedef 
			//quantity<power_typeof_helper<
			//	atomic::energy, //angstrom_unit, 
			//	static_rational<2> 
			//>::type> 
			typename multiply_typeof_helper<
				Left,
				Right
			>::type
				result_type;
		result_type operator()(Left & left, Right & right) const{
			return left*right;
		}
	};
}}}

namespace boost{namespace units{
namespace ba=boost::accumulators;
template<class Quantity, class Stats = ba::stats<ba::tag::mean, ba::tag::variance> > 
class accumulator_set{
	public:
	ba::accumulator_set<double, Stats> impl_;
	void operator()(Quantity const& q){
		return impl_(q.value());
	}
};
template<class Quantity, class Stats>
Quantity mean(accumulator_set<Quantity, Stats> const& self){
	return Quantity::from_value(ba::mean(self.impl_));
}
template<class Quantity> 
#define RET_TYPE \
	typename \
	power_typeof_helper< \
		Quantity, \
		static_rational<2> \
	>::type
RET_TYPE variance(
	accumulator_set<Quantity> const& self
){
	typedef RET_TYPE ret_type; //necesary to avoid quantity.hpp:828: error: invalid use of incomplete type in g++ 4.4.0
	return ret_type::from_value(ba::variance(self.impl_));
}
#undef RET_TYPE

}}
#endif

#ifdef _BOOST_UNITS_ACCUMULATORS_TEST
#include "boost/units/systems/atomic.hpp"
int main(){
	using namespace boost::units;
	accumulator_set< quantity<atomic::energy> /*, stats<tag::mean >*/ > acc;

	// push in some data ...
	acc(1.2*atomic::hartree);
	acc(2.3*atomic::hartree);
	//acc(3.4*atomic::hartree);
	//acc(4.5*atomic::hartree);

	// Display the results ...
	std::cout << "Mean:   " << boost::units::mean(acc) << std::endl;
	std::cout << "Var:   " << sqrt(boost::units::variance(acc)) << std::endl;
	//std::cout << "Moment: " << moment<2>(acc) << std::endl;
}
#endif
// astyle --brackets=attach --indent=tab --indent-col1-comments --delete-empty-lines --add-brackets --keep-one-line-statements --convert-tabs --align-pointer=type --unpad-paren  qbox.hpp                 
// Editor modelines  -  http://www.wireshark.org/tools/modelines.html                                                                                                                                      
// Local variables:                                                                                                                                                                                        
// c-basic-offset: 4                                                                                                                                                                                       
// tab-width: 4                                                                                                                                                                                            
// indent-tabs-mode: t                                                                                                                                                                                     
// truncate-lines: 1                                                                                                                                                                                       
// End:                                                                                                                                                                                                    
/* vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent: */

