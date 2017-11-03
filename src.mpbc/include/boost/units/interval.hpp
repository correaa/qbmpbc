#ifdef compile_instructions
ln -sf $0 .$0.cpp && c++ -I$HOME/prj -D_TEST_BOOST_UNITS_INTERVAL_HPP .$0.cpp -o ./$0.x && ./$0.x $@
rm -f .$0.x .$0.cpp; exit
#endif
#ifndef UNITS_INTERVAL_HPP_
#define UNITS_INTERVAL_HPP_
#include<boost/units/unit.hpp>
#include<boost/units/quantity.hpp>
#include<boost/units/limits.hpp>
#include <boost/numeric/interval.hpp>
//! The checking policy for boost::interval, if it uses boost::unit
//  model code taken from boost/numeric/interval/checking.hpp
namespace boost {namespace numeric {namespace interval_lib {

template<class Unit>
struct checking_base<boost::units::quantity<Unit> >{
	private:
	typedef boost::units::quantity<Unit> T; //to minimize code change
	public:
	static T pos_inf(){
		assert(std::numeric_limits<T>::has_infinity);
		return std::numeric_limits<T>::infinity();
	}
	static T neg_inf(){
		assert(std::numeric_limits<T>::has_infinity);
		return -std::numeric_limits<T>::infinity();
	}
	static T nan(){
		assert(std::numeric_limits<T>::has_quiet_NaN);
		return std::numeric_limits<T>::quiet_NaN();
	}
	static bool is_nan(const T& x){
		return std::numeric_limits<T>::has_quiet_NaN && (x != x);
	}
	static T empty_lower(){
		return (
			std::numeric_limits<T>::has_quiet_NaN ?
			std::numeric_limits<T>::quiet_NaN() : T::from_value(1) // <-- only change, was static_cast<T>(1)
		);
	}
	static T empty_upper(){
		return (
			std::numeric_limits<T>::has_quiet_NaN ?
			std::numeric_limits<T>::quiet_NaN() : static_cast<T>(0) // this magically still works
		);
	}
	static bool is_empty(const T& l, const T& u){
		return !(l <= u); // safety for partial orders
	}
};
}

//< overloads boost/numeric/interval/utility.hpp:80 for division by integer
template<class Unit>
boost::units::quantity<Unit> median(const interval<boost::units::quantity<Unit> >& iv){
	return (iv.lower()+iv.upper())/2.;
}

using namespace boost::units;
// declaration to allow 
#define DECLARE_OPERATOR( SymboL, NamE) \
template<typename T, class Unit> \
interval<typename NamE##_typeof_helper<T, quantity<Unit> >::type> \
operator SymboL(interval<T> const& iv, quantity<Unit> const& q){ \
	return interval<typename multiply_typeof_helper<T, quantity<Unit> >::type>( \
		iv.lower() SymboL q, iv.upper() SymboL q \
	); \
} \
template<typename T, class Dim,class System> \
interval<typename NamE##_typeof_helper<T, unit<Dim, System> >::type> \
operator SymboL(interval<T> const& iv, unit<Dim, System> const& q){ \
	return interval<typename NamE##_typeof_helper<T, unit<Dim, System> >::type>( \
		iv.lower() SymboL q, iv.upper() SymboL q \
	); \
}
DECLARE_OPERATOR(* , multiply);
DECLARE_OPERATOR(/ , divide);
#undef DECLARE_OPERATOR

}}
#endif
#ifdef _TEST_BOOST_UNITS_INTERVAL_HPP
#include<iostream>
#include<boost/numeric/interval/io.hpp>
#include<boost/units/systems/si.hpp>
#include "boost/numeric/make_interval.hpp"
int main(){
	using namespace boost::units;
	std::clog << boost::numeric::interval<double>(0.,1.)*2. << std::endl;
	std::clog << boost::numeric::interval<double>(0.,1.)*si::meter << std::endl;
	std::clog << (boost::numeric::make_interval(0.,1.)/si::meter).lower() << std::endl;
	return 0;
}
#endif

