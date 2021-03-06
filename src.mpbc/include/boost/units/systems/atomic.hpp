#ifdef compile_instructions
ln -sf $0 .$0.cpp && c++ `#-Wfatal-errors` `#-std=c++0x` .$0.cpp -o .$0.x -D_ATOMIC_HPP_TEST -I../../.. && ./.$0.x $@; exit;
#endif
#ifndef BOOST_UNITS_SYSTEMS_ATOMIC
#define BOOST_UNITS_SYSTEMS_ATOMIC
//file <unspecified>/units/systems/atomic.hpp

//do not forget to include, so this header is useful
//#include<boost/units/quantity.hpp>??
//#include"boost/units/systems/atomic.hpp" <-this unit system
//#include"boost/units/systems/atomic/si_conversion.hpp" <- to make conversions to SI system

#include <boost/units/base_unit.hpp>

//equivalent to includes in units/systems/si/base.hpp
#include <boost/units/static_constant.hpp>
#include <boost/units/unit.hpp>
#include <boost/units/make_system.hpp>
 
//#include <boost/units/io.hpp> //do not remove: '...has initializer but incomplete type...' error
 
//equivalent to includes in //equivalent to ./base_units/si/*.hpp
#include <boost/units/config.hpp>
#include <boost/units/base_unit.hpp>
#include <boost/units/scaled_base_unit.hpp>
#include <boost/units/physical_dimensions/mass.hpp>
#include <boost/units/physical_dimensions/length.hpp>
#include <boost/units/physical_dimensions/time.hpp>
#include <boost/units/physical_dimensions/velocity.hpp>
#include <boost/units/physical_dimensions/current.hpp>
#include <boost/units/physical_dimensions/electric_charge.hpp>
#include <boost/units/physical_dimensions/angular_momentum.hpp>
#include <boost/units/physical_dimensions/force.hpp>
#include <boost/units/physical_dimensions/area.hpp>
#include <boost/units/physical_dimensions/volume.hpp>
#include <boost/units/physical_dimensions/pressure.hpp>
#include <boost/units/physical_dimensions/wavenumber.hpp>
#include <boost/units/physical_dimensions/frequency.hpp>
#include <boost/units/physical_dimensions/energy.hpp>
#include <boost/units/physical_dimensions/action.hpp>
#include <boost/units/physical_dimensions/temperature.hpp>
#include <boost/units/physical_dimensions/heat_capacity.hpp>
#include <boost/units/derived_dimension.hpp>
//begin equivalent to ./base_units/si/*.hpp
//for example boost/units/base_units/atomic/electron_mass.hpp
namespace boost{
namespace units{
/// Atomic system of units, based in the definitions in http://en.wikipedia.org/wiki/Atomic_units#Fundamental_units
namespace atomic{   //define namespace au=boost::units::atomic if necessary 
	using std::string;
//mass
/*1*/	
	struct electron_mass_base_unit : base_unit<electron_mass_base_unit, mass_dimension, 1>{
		static string name()   {return ("electron mass");}
		static string symbol() {return ("\\mathit{m}_e");}
	};
//charge
/*2*/ 	
	struct elementary_charge_base_unit : base_unit<
		elementary_charge_base_unit, 
		electric_charge_dimension,//derived_dimension<time_base_dimension,1, current_base_dimension,1>::type,
		2
	>{
 		static string name()   {return ("elementary charge"/*wiki friendly*/);} // "electron charge"/*wolfram friendly*/
 		static string symbol() {return ("\\mathit{e}");}
 	};
	//action (not angular momentum)
/*3*/ 	struct reduced_planck_constant_base_unit : base_unit<reduced_planck_constant_base_unit, 
		action_dimension, //angular_momentum_dimension, <- in boost.units action!=angular momentum by plane-angle unit
		//derived_dimension<time_base_dimension, -1, length_dimension, 2, mass_dimension, 1>::type,
		3
	>{
 		static string name(){return ("reduced Planck constant");}
 		static string symbol() {return ("\\mathit{\\hbar}");} // â
 	};
	//electric constant
	struct coulomb_force_constant_base_unit : base_unit<
		coulomb_force_constant_base_unit,
		derived_dimension<
			length_base_dimension,   3,
			mass_base_dimension,     1,
			time_base_dimension,    -4,
			current_base_dimension, -2
		>::type,
		4
	>{
 		static string name()   {return ("(Coulomb constant value)"/*wolfram friendly*/);} //("coulomb force constant"/*wiki friendly*/);}
 		static string symbol() {
			return 
				//"1/(4 đÂ đÂâ)"// wiki friendly unicode // not friendly with acrobat annotations
				"1/(4 \\pi \\varepsilon_0)" //wiki friendly latex
				// "k_1" /*jackson friendly*/ 
				// "đÂ"/*wolfram friendly*/
		;
		}
	};
	//other options for base units
 	//struct hartree_energy_base_unit : base_unit<hartree_energy_base_unit, energy_dimension, 4>{
 	//	static string name()   {return ("hartree_energy");}
 	//	static string symbol() {return ("E_h");}
 	//};
	//equivalent to ./base_units/si/meter.hpp>
	//struct bohr_radius_base_unit : base_unit<bohr_radius_base_unit, length_dimension, 5>{
	//	static string name()	  {return ("bohr_radius");}
 	//	static string symbol() {return ("a_0");}
 	//};
	struct 
		boltzman_constant_base_unit : base_unit<
		boltzman_constant_base_unit, 
			derived_dimension<
				length_base_dimension, 2,
				mass_base_dimension, 1,
				time_base_dimension, -2, // == ,energy_dimension, 1,
				temperature_base_dimension, -1
			>::type,
			5
		>{
		static string name()  {return ("boltzman constant");}
		static string symbol(){return ("\\mathit{k}_B");}
	};    
}}}//nss
//end equivalent to ./base_units/si/*.hpp

//begin boost/units/system/atomic/base.hpp (similar to ./systems/si/base.hpp)
namespace boost{
namespace units{
namespace atomic{
 	typedef make_system<
 		electron_mass_base_unit,
 		elementary_charge_base_unit,
 		reduced_planck_constant_base_unit,
		coulomb_force_constant_base_unit,
 		//hartree_energy_base_unit,
 		//bohr_radius_base_unit,
	    boltzman_constant_base_unit
 	>::type system;
}}}
//end boost/units/system/atomic/base.hpp (similar to ./systems/si/base.hpp)

namespace boost{
namespace units{
template<class System>
struct dimensionless{
	typedef unit<dimensionless_type, System> unit_type; 
	//static unit<dimensionless_type, System> unit_;
};
}
}


//begin equivalent to boost/units/systems/si/*.hpp
//for example boost/units/system/atomic/mass|length|electric_charge.hpp
namespace boost{
namespace units{
namespace atomic{
 	typedef unit<dimensionless_type, system> dimensionless;
/*1*/ 	
	typedef unit<mass_dimension,             system> mass;
/*2*/
	typedef unit<electric_charge_dimension,  system> electric_charge;
/*3*/
	typedef unit<action_dimension, system>           action;
/*4*/ // http://www.bipm.org/en/si/si_brochure/chapter4/4-1.html doesn't like this as a base unit
	typedef 
		unit<derived_dimension<
			length_base_dimension,  3,
			mass_base_dimension,    1,
			time_base_dimension,   -4,
			current_base_dimension,-2
		>::type, system> 
		electric_constant;
/*5*/
	typedef unit<heat_capacity_dimension, system> heat_capacity;
	typedef heat_capacity entropy;
	//boltzman_constant;
	//typedef boltzman_constant heat_capacity;
/*1*/ 	static const mass              electron_rest_mass,      electron_mass, m_e;
/*2*/ 	static const electric_charge   elementary_charge,       eminus,        e; 
/*3*/ 	static const action            reduced_planck_constant, hbar;
/*4*/	static const electric_constant coulomb_force_constant,  k_1 /*name taken from Jackson p. 779*/, k /*Wolfram Alpha name is kappa*/;
/*5*/	static const heat_capacity     boltzman_constant,       k_B,           kB;
}}}
//end equivalent to units/systems/si/*.hpp

namespace boost{
namespace units{
namespace atomic{
	typedef unit<energy_dimension,      system> energy;
	typedef unit<length_dimension,      system> length;
	typedef unit<force_dimension,       system> force;
	typedef unit<area_dimension,        system> area;
	typedef unit<volume_dimension,      system> volume;
	typedef unit<pressure_dimension,    system> pressure;
	typedef unit<time_dimension,        system> time;
 	typedef unit<velocity_dimension,    system> velocity;
	typedef unit<temperature_dimension, system> temperature;
	typedef unit<wavenumber_dimension,  system> wavenumber; typedef wavenumber inverse_length; 
	typedef unit<frequency_dimension,   system> frequency;
	typedef unit<temperature_dimension, system> temperature;
	typedef	
		unit<
			derived_dimension<
				length_base_dimension, -3
			>::type, system
		> 
		number_density
	;
	typedef multiply_typeof_helper<electric_charge, number_density>::type electric_charge_density;
	typedef	
		power_typeof_helper<number_density, static_rational<1,2> >::type
		wavefunction_amplitude
	;
	static const time     time_unit;
	static const velocity velocity_unit;
 	static const energy   energy_unit,  hartree_energy, hartree, E_h; 
	static const length   length_unit,  bohr_radius,    bohr,    a_0, B;
}}}

#ifndef BOOST_LATEX_HPP
#include <boost/units/io.hpp>
namespace boost{
namespace units{
//template<class Dimension>//, class System>
//std::string symbol_latex(const typename reduce_unit<Dimension>::type& u){return "{"+symbol_string(u)+"}";}
//}}
namespace detail {
struct format_raw_latex_impl {
    template<class Units>
    void append_units_to(std::string& str) {
        detail::symbol_string_impl<Units::size::value>::template apply<Units>::value(str);
    }
    template<class Scale>
    void append_scale_to(std::string& str) {
        detail::scale_symbol_string_impl<Scale::size::value>::template apply<Scale>::value(str);
    }
    template<class Unit>
    std::string operator()(const Unit& unit) {
        return(to_string_impl(unit, *this));
    }
    template<class Unit>
    bool is_default_string(const std::string&, const Unit&) {
        return(true);
    }
};
struct format_latex_impl : format_raw_latex_impl {
    template<class Unit>
    std::string operator()(const Unit& unit) {
        return(symbol_latex(unit));
    }
    template<class Unit>
    bool is_default_string(const std::string& str, const Unit& unit) {
        return(str == to_string_impl(unit, format_raw_symbol_impl()));
    }
};
// These two overloads of symbol_string and name_string will
// will pick up homogeneous_systems.  They simply call the
// appropriate function with a heterogeneous_system.
template<class Dimension,class System, class SubFormatter>
inline std::string
to_latex_impl(const unit<Dimension,System>&, SubFormatter f){
    return f(typename reduce_unit<unit<Dimension, System> >::type());
}
}
template<class Dimension,class System>
inline std::string symbol_latex(const unit<Dimension, System>&){
    return detail::to_latex_impl(unit<Dimension,System>(), detail::format_latex_impl());
}
}}

#endif

//begin boost/system/atomic/io.hpp
#include <boost/units/io.hpp>
#include <boost/units/reduce_unit.hpp>
// nice names for commonly used units
namespace boost {
namespace units {
namespace atomic{
	using std::string;
	//it is useless (confusing?) to add strings for base_units, do it only for derived units
	inline string name_string  (const reduce_unit<time>::type&                   ){return "atomic unit of time";       } //gives error with boost 1.37
	inline string symbol_string(const reduce_unit<time>::type&                   ){return "â/E_h";                     }
	inline string name_string  (const reduce_unit<energy>::type&                 ){return "hartree";                   }
	inline string symbol_string(const reduce_unit<energy>::type&                 ){return "E_h" /*"E_\\mathrm{h}"*/;   }
	inline string name_string  (const reduce_unit<length>::type&                 ){return "Bohr radius";               } 
	inline string symbol_string(const reduce_unit<length>::type&                 ){return "a_0"; /*"aâ";*/             }
	inline string name_string  (const reduce_unit<volume>::type&                 ){return "atomic unit of volume";     } 
	inline string symbol_string(const reduce_unit<volume>::type&                 ){return symbol_string(reduce_unit<length>::type()) + "^3"  /*"aâ^3"*/;                 }
	inline string name_string  (const reduce_unit<number_density>::type&         ){return "atomic unit of number density";}
	inline string symbol_string(const reduce_unit<number_density>::type&         ){return symbol_string(reduce_unit<length>::type()) + "^-3" /*"aâ^{-3}"*/;                }
	inline string name_string  (const reduce_unit<velocity>::type&               ){return "atomic unit of velocity";                                     }
	inline string symbol_string(const reduce_unit<velocity>::type&               ){return symbol_string(reduce_unit<length>::type()) + " E_h/â" /*"aâ E_h/â"*/;         }
	inline string name_string  (const reduce_unit<temperature>::type&            ){return "atomic unit of energy per boltzman constant";} //"atomic unit of temperature";}
	inline string symbol_string(const reduce_unit<temperature>::type&            ){return "\\mathit{E}_h/\\mathit{k}_B";                }
	inline string symbol_string(const reduce_unit<wavenumber>::type&             ){return symbol_string(reduce_unit<length>::type()) + "^-1";}
	inline string symbol_string(const reduce_unit<frequency>::type&              ){return "E_h/â";}
	inline string symbol_string(const reduce_unit<force>::type&                  ){return "E_h/" + symbol_string(reduce_unit<length>::type()) /*"E_h/aâ"*/;}
	inline string name_string  (const reduce_unit<wavefunction_amplitude>::type& ){return "atomic quantum-mechanical wave function unit";}
	inline string symbol_string(const reduce_unit<wavefunction_amplitude>::type& ){return symbol_string(reduce_unit<length>::type()) + "^-3/2";}
	inline string name_string  (const reduce_unit<electric_charge_density>::type&){return "atomic unit of electric charge density";}
	inline string symbol_string(const reduce_unit<electric_charge_density>::type&){return "e/" + symbol_string(reduce_unit<length>::type()) + "^3";}
}}
}
//end boost/system/atomic/io.hpp
//begin boost/system/si/iox.hpp
#ifndef SI_IOX
#define SI_IOX
#include<boost/units/systems/si/resistivity.hpp>
namespace boost{
namespace units{
namespace si{
	using std::string;
	//it is useless (confusing?) to add strings for base_units, do it only for derived units
	inline string name_string  (const reduce_unit<resistivity>::type&       ){return "si unit of resistivity";}
	inline string symbol_string(const reduce_unit<resistivity>::type&       ){return "âŚm"; /*was "đş m";*/                    }
}}}
#endif

namespace boost{namespace units{
// Non-atomic units, commonly used along with atomic units but not part of the system
namespace nonatomic{ 
	struct rydberg_base_unit : base_unit<rydberg_base_unit, energy_dimension, 76345>{
		static const std::string name() {return "rydberg";}
		static const std::string symbol(){return "Ry";}
	};
	typedef rydberg_base_unit::unit_type rydberg_unit;
	static const rydberg_unit rydberg, Ry;
	//inline std::string symbol_string(const reduce_unit<rydberg_unit>::type&){return "Ry";}
}
}}
BOOST_UNITS_DEFINE_CONVERSION_FACTOR(
	boost::units::nonatomic::rydberg_base_unit,
	boost::units::atomic::energy, //hartree
	double, 0.5 //exact
);
BOOST_UNITS_DEFAULT_CONVERSION(
	boost::units::nonatomic::rydberg_base_unit, 
	boost::units::atomic::energy  //hartree
);

//begin of boost/units/systems/nonsi.hpp
#include <boost/units/base_units/si/meter.hpp>
#include"boost/units/systems/si/io.hpp"
//#include<boost/units/base_units/metric/angstrom.hpp> //for "A" output
namespace boost{namespace units{ //taken from http://groups.google.com/group/boostusers/msg/7ecfddc90900717c?
/// Non-SI units but accepted for use with the SI http://www.bipm.org/en/si/si_brochure/chapter4/4-1.html
namespace nonsi{
	struct bar_base_unit : base_unit<bar_base_unit, pressure_dimension, 94320>{
		static const std::string name()   { return "bar";}
		static const std::string symbol() { return "bar";}
	};
	typedef bar_base_unit::unit_type bar_unit;
	static const bar_unit bar; 

	struct electron_volt_base_unit : base_unit<electron_volt_base_unit, energy_dimension, 43289> {
		static const std::string name()   { return "electronvolt";}
		static const std::string symbol() { return "eV"; }
	};
	typedef electron_volt_base_unit::unit_type electron_volt_unit;
	typedef electron_volt_base_unit::unit_type electronvolt_unit;
	static const electron_volt_unit electron_volt, electronvolt, eV;
	inline std::string symbol_string(const reduce_unit<electron_volt_unit>::type&){return "eV";}

	struct angstrom_base_unit : base_unit<angstrom_base_unit, length_dimension, 33485> {
		static const std::string name()   {return "angstrom";}
		static const std::string symbol() {return "âŤ";}
	};
	//typedef scaled_base_unit<boost::units::si::meter_base_unit, scale<10, static_rational<-10> > > angstrom_base_unit;
	typedef angstrom_base_unit::unit_type angstrom_unit;
	static const angstrom_unit angstrom, A;

	struct dalton_base_unit : base_unit<dalton_base_unit, mass_dimension, 429458> {
		static const std::string name()   { return "dalton";}
		static const std::string symbol() { return "Da"; }
	};
	typedef dalton_base_unit::unit_type dalton_unit;
	static const dalton_unit dalton, daltons, unified_atomic_mass_unit, Da, u; //notation from http://www.bipm.org/en/si/si_brochure/chapter4/table7.html 
	inline std::string symbol_string(const reduce_unit<dalton_unit>::type&){return "Da";}
}
namespace si{
	typedef scaled_base_unit<boost::units::si::second_base_unit, scale<10, static_rational<-15> > > femtosecond_base_unit;
	typedef femtosecond_base_unit::unit_type femtosecond_unit;
	static const femtosecond_unit femtosecond, femtoseconds, fs;
	typedef scaled_base_unit<boost::units::si::second_base_unit, scale<10, static_rational<-12> > > picosecond_base_unit;
	typedef picosecond_base_unit::unit_type picosecond_unit;
	static const picosecond_unit picosecond, picoseconds, ps;

}
}}
BOOST_UNITS_DEFINE_CONVERSION_FACTOR(
	boost::units::nonsi::bar_base_unit,
	boost::units::si::pressure, 
	double, 100000. 
);
BOOST_UNITS_DEFAULT_CONVERSION(
	boost::units::nonsi::bar_base_unit,
	boost::units::si::pressure 
);
BOOST_UNITS_DEFINE_CONVERSION_FACTOR(
	boost::units::nonsi::electron_volt_base_unit,
	boost::units::atomic::energy, 
	double, 0.03674932 //should be defined from available codata, for example electron_volt.value()/atomic::hartree or something
);
BOOST_UNITS_DEFAULT_CONVERSION(
	boost::units::nonsi::electron_volt_base_unit, 
	boost::units::atomic::energy
);
BOOST_UNITS_DEFINE_CONVERSION_FACTOR( //since the conversion is defined for si::length, to convert from angstrom to atomic::borh #include "atomic/si_conversion.hpp"
	boost::units::nonsi::angstrom_base_unit,
	boost::units::si::length, 
	double, 1.e-10 //should be defined from available codata, for example electron_volt.value()/atomic::hartree or something
);
BOOST_UNITS_DEFINE_CONVERSION_FACTOR(
	boost::units::nonsi::dalton_base_unit,
	boost::units::si::mass, 
	double, 1.6605388628e-27 //should be defined from available codata, for example electron_volt.value()/atomic::hartree or something
);

//end of boost/units/systems/nonsi.hpp

//begin of boost/units/power.hpp
namespace boost{
namespace units{
template<class BaseUnit, class Rational>
struct power{
  typedef boost::units::unit<
    typename boost::units::static_power<typename BaseUnit::dimension_type, Rational >::type, 
    typename boost::units::make_system<BaseUnit>::type
  > unit;
};
template<class BaseUnit, int N>
struct power_int{
	typedef typename power<BaseUnit, typename static_rational<N,1>::type >::unit unit;
};
}}
//end of boost/units/power.hpp

//begin of boost/units/generic_programming.hpp
#ifndef BOOST_UNITS_GENERIC_PROGRAMMING_HPP
#define BOOST_UNITS_GENERIC_PROGRAMMING_HPP
namespace boost{namespace units{
//value(q) is more generic than q.value() since it works for not unit types
template<class Unit, class T> T const& value(quantity<Unit, T> const& q){return q.value();}
template<            class T> T const& value(               T  const& v){return v        ;}
}}
#endif
//end of boost/units/generic_programming.hpp

//end of boost/units/static_power.hpp

//begin of boost/units/scale_type.hpp
//#include <boost/units/scaled_base_unit.hpp>
//#include <boost/units/scale.hpp>
#include<boost/units/make_scaled_unit.hpp>
namespace boost{
namespace units{
namespace si{
template<class Unit> struct femto_scaled{typedef typename make_scaled_unit<Unit, scale<10, static_rational<-15> > >::type unit;};
template<class Unit> struct pico_scaled {typedef typename make_scaled_unit<Unit, scale<10, static_rational<-12> > >::type unit;};
template<class Unit> struct nano_scaled {typedef typename make_scaled_unit<Unit, scale<10, static_rational< -9> > >::type unit;};
template<class Unit> struct mega_scaled {typedef typename make_scaled_unit<Unit, scale<10, static_rational< +6> > >::type unit;};
template<class Unit> struct giga_scaled {typedef typename make_scaled_unit<Unit, scale<10, static_rational< +9> > >::type unit;};
}
}}

namespace boost{
namespace units{
	template<class U1, class U2>
	struct divides{
		typedef typename divide_typeof_helper<U1, U2>::type type;
	}; 
}
}


namespace boost{
namespace units{
namespace auto_conversion_operators{
	template<class LeftDimension, class LeftSystem, class RightUnit>
	#define RET_TYPE \
	quantity< \
		typename multiply_typeof_helper< \
			boost::units::unit<LeftDimension , LeftSystem>, \
			boost::units::unit<typename RightUnit::dimension_type, LeftSystem> \
		>::type \
	>
	RET_TYPE operator*(
		unit<LeftDimension, LeftSystem> const& u1,
		quantity<RightUnit> const& q2
	){
		return RET_TYPE(boost::units::operator*(u1,q2));
	}
	#undef RET_TYPE

	template<class LeftUnit, class RightDimension, class RightSystem>
	#define RET_TYPE \
	quantity< \
		typename multiply_typeof_helper< \
			boost::units::unit<typename LeftUnit::dimension_type, RightSystem>, \
			boost::units::unit<RightDimension, RightSystem> \
		>::type \
	>
	RET_TYPE operator*(
		quantity<LeftUnit> const& q1,
		unit<RightDimension, RightSystem> const& u2
	){
		return RET_TYPE(boost::units::operator*(q1, u2));
	}
	#undef RET_TYPE

	template<class LeftUnit, class RightDimension, class RightSystem>
	#define RET_TYPE \
	quantity< \
		typename divide_typeof_helper< \
			boost::units::unit<typename LeftUnit::dimension_type, RightSystem>, \
			boost::units::unit<RightDimension, RightSystem> \
		>::type \
	>
	RET_TYPE operator/(
		quantity<LeftUnit> const& q1,
		unit<RightDimension, RightSystem> const& u2
	){
		return RET_TYPE(boost::units::operator/(q1, u2));
	}
	#undef RET_TYPE

// only useful for * (multiply), and / (divide)
#define DECLARE_AUTO_CONVERSION_OPERATOR(preferredSide, operatorSymbol, operatorName) \
template<class LeftUnit, class RightUnit> \
quantity< \
	typename operatorName##_typeof_helper< \
		boost::units::unit<typename LeftUnit::dimension_type , typename preferredSide::system_type>, \
		boost::units::unit<typename RightUnit::dimension_type, typename preferredSide::system_type>  \
	>::type \
> operator operatorSymbol ( \
	quantity<LeftUnit>  const& t1, \
	quantity<RightUnit> const& t2  \
){ \
	return \
		quantity< \
			typename operatorName##_typeof_helper< \
				boost::units::unit<typename LeftUnit::dimension_type , typename preferredSide::system_type>, \
				boost::units::unit<typename RightUnit::dimension_type, typename preferredSide::system_type> \
			>::type \
		>( \
			boost::units::operator operatorSymbol (t1,t2) \
		) \
	; \
}
namespace to_left{
	DECLARE_AUTO_CONVERSION_OPERATOR(LeftUnit, *, multiply)
	DECLARE_AUTO_CONVERSION_OPERATOR(LeftUnit, /, divide)
	template<class Dimension, class LeftSystem, class RightSystem>
	quantity< unit<Dimension, LeftSystem> > operator+(
		quantity< unit<Dimension, LeftSystem> > const& qleft, 
		quantity< unit<Dimension, RightSystem> > const& qright
	){
		return boost::units::operator + (qleft, quantity< unit<Dimension, LeftSystem> >(qright));
	}
	template<class Dimension, class LeftSystem, class RightSystem>
	quantity< unit<Dimension, LeftSystem> > operator - (
		quantity< unit<Dimension, LeftSystem> > const& qleft, 
		quantity< unit<Dimension, RightSystem> > const& qright
	){
		return boost::units::operator - (qleft, quantity< unit<Dimension, LeftSystem> >(qright));
	}
}
namespace to_right{
	DECLARE_AUTO_CONVERSION_OPERATOR(RightUnit, *, multiply)
	DECLARE_AUTO_CONVERSION_OPERATOR(RightUnit, /, divide)
	DECLARE_AUTO_CONVERSION_OPERATOR(RightUnit, +, add)
	template<class Dimension, class LeftSystem, class RightSystem>
	quantity< unit<Dimension, RightSystem> > operator+(
		quantity< unit<Dimension, LeftSystem> > const& qleft, 
		quantity< unit<Dimension, RightSystem> > const& qright
	){
		return boost::units::operator + (quantity< unit<Dimension, RightSystem> >(qleft), qright);
	}
	template<class Dimension, class LeftSystem, class RightSystem>
	quantity< unit<Dimension, RightSystem> > operator-(
		quantity< unit<Dimension, LeftSystem> > const& qleft, 
		quantity< unit<Dimension, RightSystem> > const& qright
	){
		return boost::units::operator - (quantity< unit<Dimension, RightSystem> >(qleft), qright);
	}
}
#undef DECLARE_AUTO_CONVERSION_MULTIPLY
}
}}



#ifdef _ATOMIC_HPP_TEST 

#include"boost/units/systems/si/io.hpp"
#include<boost/units/systems/si/length.hpp>
#include<boost/units/base_units/metric/angstrom.hpp>
namespace boost { namespace units
{
	//typedef metric::angstrom_base_unit::unit_type angstrom_unit;
	//static const angstrom_unit angstrom; //BOOST_UNITS_STATIC_CONSTANT(angstrom, angstrom_unit); 
}}

#include<iostream>
//#include<boost/units/quantity.hpp>
//#include"boost/units/systems/atomic.hpp"
#include<boost/units/systems/si/prefixes.hpp>
#include<boost/units/systems/si/pressure.hpp>
#include<boost/units/systems/si/time.hpp>
#include"boost/units/systems/si/io.hpp"
#include"boost/units/systems/atomic/si_conversion.hpp"
#include<boost/units/systems/si/codata/physico-chemical_constants.hpp> //for m_u (dalton unit, http://en.wikipedia.org/wiki/Atomic_mass_unit) and k_B (boltzman constant)

#include<boost/units/make_scaled_unit.hpp>
#include<boost/units/cmath.hpp>
#include <boost/units/scaled_base_unit.hpp>
#include <boost/units/static_rational.hpp>
#include <boost/units/scale.hpp>
#include <boost/units/units_fwd.hpp>
#include <boost/units/base_units/si/meter.hpp>

//#include<boost/units/base_units/metric/angstrom.hpp>
using std::cout;
using std::clog;
using std::cerr;
using std::endl;


namespace boost{
namespace units{
namespace atomic{
namespace si_conversion{

template<class Dimension, class Y>
quantity<unit<Dimension, atomic::system>, Y> operator+(
	quantity<unit<Dimension, atomic::system>, Y> const& e1, 
	quantity<unit<Dimension, si::system    >, Y> const& e2
){
	return e1+quantity<unit<Dimension, atomic::system>, Y>(e2); // or e1 + (quantity<atomic::energy, Y>)e2; ???
}
template<class Dimension, class Y>
quantity<unit<Dimension, atomic::system>, Y> operator+(
	quantity<unit<Dimension, si::system    >, Y> const& e1,
	quantity<unit<Dimension, atomic::system>, Y> const& e2
){ 
	return quantity<unit<Dimension, atomic::system>, Y>(e1)+e2;
}
template<class Dimension, class Y>
quantity<unit<Dimension, atomic::system>, Y> operator-(
	quantity<unit<Dimension, atomic::system>, Y> const& e1, 
	quantity<unit<Dimension, si::system    >, Y> const& e2
){
	return e1-quantity<unit<Dimension, atomic::system>, Y>(e2);
}
template<class Dimension, class Y>
quantity<unit<Dimension, atomic::system>, Y> operator-(
	quantity<unit<Dimension, si::system>, Y> const& e1,
	quantity<unit<Dimension, atomic::system>, Y> const& e2
){ 
	return quantity<unit<Dimension, atomic::system>, Y>(e1)-e2;
}
}}}}
#undef DECLARE_AUTO_CONVERSION_MULTIPLY

#include<boost/mpl/divides.hpp>
#include<boost/mpl/times.hpp>

#include<boost/numeric/interval/io.hpp>
#include "boost/units/interval.hpp"
#include "boost/lambda/lambda.hpp"
#include "boost/units/lambda.hpp"
#include<boost/spirit/home/phoenix.hpp>
#include<boost/spirit/home/phoenix/core/argument.hpp>
#include<boost/units/systems/cgs.hpp>
#include "boost/units/systems/si/heat_capacity.hpp"
#include <boost/units/systems/cgs.hpp>

int main(){
	using namespace boost::units;

	{
		//quantity<divide_typeof_helper<si::mass, si::length>::type> ld;
		quantity<divides<si::mass, si::length>::type> ld;
		cout << "ld = "<< ld <<endl;
	}
	
	quantity<si::energy> een1(1.e-17*si::joules);
	quantity<atomic::energy> een2(2.*atomic::hartree);
	
	
	cout<< (boost::lambda::_1+boost::lambda::_1)(een1) <<endl;
	cout<< (boost::phoenix::arg_names::arg1 + boost::phoenix::arg_names::arg1)(een1) <<endl;

	boost::numeric::interval<quantity<atomic::energy> > ivl(2.*atomic::hartree, 3.*atomic::hartree);
	clog<<ivl<<endl;

	using namespace boost::units::atomic::si_conversion;	
	quantity<atomic::energy> een3(een1+een2);
	
	cout<<een3<<endl;
	//std::cout<<(31.*atomic::length_unit)<<std::endl;
	quantity<atomic::length> l(1.*atomic::length_unit);
	std::clog<<l<<std::endl;
	quantity<atomic::energy> en(1.*atomic::energy_unit);
	clog<<atomic::hbar/en<<endl;

	quantity<atomic::time> tau(1.*atomic::hbar/atomic::hartree); //tau = 1.*atomic::hbar/atomic::hartree;
	clog<<tau<<endl;
	std::cout<<((atomic::hbar/atomic::E_h)/atomic::time_unit)<<std::endl;

	quantity<atomic::time> dt = 2.*atomic::time_unit;
	quantity<atomic::length> dl = 2.*atomic::a_0;
	std::cout<<dt<<", "<<dl<<", "<<dl/dt<<std::endl;
	quantity<atomic::velocity> v(dl/dt);
	std::cout<< v <<std::endl;
	quantity<atomic::mass> m_u(si::constants::codata::m_u.value());
	std::cout<<"m_u in au"<<m_u/atomic::m_e<<std::endl;
	//std::cout<<"kinetic "<< v*v*m_u / si::constants::codata::k_B.value() <<std::endl;
	
	quantity<si::mass> m_e_si(1.*atomic::m_e); std::clog<<m_e_si<<std::endl;
	quantity<si::temperature> ttsi(4.*si::kelvin); 	std::cout<<"ttsi "<<ttsi<<std::endl;
	quantity<boltzman_constant_si_unit> dd(1.*atomic::k_B);
	std::clog<<"k_B = "<<3.*atomic::k_B<< " "<<dd<<std::endl;
	quantity<atomic::electric_charge> ee = 1.*atomic::eminus; std::cout<<"ee = "<<ee<<std::endl;
	quantity<si::electric_charge> eesi(ee); cout<< eesi << endl;

	quantity<atomic::energy> eh = 1.*atomic::E_h; cout<<eh<<endl;
	quantity<si::energy> ehsi(eh); cout<< ehsi <<endl;

	quantity<atomic::temperature> at(1.*si::kelvin); cout<<at<<endl;
	quantity<atomic::temperature> tt(v*v*m_u / si::constants::codata::k_B.value()); cout<<"tt "<<tt<<endl;
	quantity<si::temperature> ttSI(tt);cout<<"ttSI "<<ttSI<<endl;

	quantity<si::temperature> ret(eh / si::constants::codata::k_B.value()); cout<<"ret "<<ret<<endl;

	//	std::cout<<"dd "<<twm_u<<std::endl;
	quantity<si::action> hbarsi(1.*atomic::hbar); cout<<"hbarsi "<<hbarsi<<endl;

	quantity<electric_constant_si_unit> one_over_shit(1.*atomic::k_1); cout<<one_over_shit<<endl;

	quantity<si::pressure> pr(1.*atomic::E_h/atomic::a_0/atomic::a_0/atomic::a_0); cout<<pr<<endl;

	quantity<atomic::length> ldd=1.*atomic::a_0;
	clog<<ldd*ldd*ldd<<endl;
	//quantity<si::energy> ehhh(1.*atomic::E_h); clog<<"ehhh = "<<ehhh<<endl;
	clog<< quantity<si::temperature>(2.*0.138259*1.*atomic::E_h/105./3. / si::constants::codata::k_B.value()) <<endl;

	typedef boost::units::make_scaled_unit<boost::units::si::pressure, boost::units::scale<10, boost::units::static_rational<  9> > >::type gpa_unit;

	typedef scaled_base_unit<boost::units::si::meter_base_unit, scale<10, static_rational<-10> > > angstrom_base_unit;
	typedef metric::angstrom_base_unit::unit_type angstrom_unit;
	static const angstrom_unit angstrom;

	static const gpa_unit gigapascal, GPa; //maybe define conversion factor to pascal or atomic::pressure
	typedef boost::units::quantity<gpa_unit> pressure;
	quantity<gpa_unit> lang(1.*si::giga*si::pascal);
	clog<<lang<<endl;

	quantity<atomic::length> angs(1.*angstrom);
	quantity<angstrom_unit> bohhr(1.*atomic::bohr);
	quantity<atomic::length> lag(10.*atomic::bohr);
	clog<< (double)(quantity<angstrom_unit>(1.*atomic::bohr)/angstrom) <<endl;

	//	quantity<nonsi::electron_volt_unit> ene = 1.*nonsi::eV;
	//cout<<(double)quantity<atomic::dimensionless>(ene/atomic::hartree)<<endl;
	cout<<simplify_typename(4.*si::nano*si::meter)<<std::endl;
	//cout<<ene<<std::endl;
	
//	cout<<simplify_typename(4.*si::nano*atomic::eV)<<std::endl;	
//	quantity<atomic::electron_volt_unit> ene2(1.*si::mega*atomic::eV);
//	cout<< 1.*si::mega*atomic::eV  <<std::endl;
	
//	quantity<atomic::energy> one_ev(1.*atomic::eV);
	quantity<atomic::energy> one_MeV(1.*si::mega*nonsi::eV);

	//cout<< _one_ev<<std::endl;
	quantity<si::femto_scaled<si::time>::unit> one_fs(1.*si::femto*si::second);
	quantity<si::femto_scaled<si::time>::unit> one_s(1.*si::second);

	quantity<si::nano_scaled<atomic::energy>::unit> one_nanohartree(1.*si::nano*atomic::hartree);
	quantity<si::nano_scaled<atomic::energy>::unit> one_hartree(1.*atomic::hartree);
	
	quantity<si::mega_scaled<nonsi::electron_volt_unit>::unit> one_Mev(1.*si::mega*nonsi::eV);
	
	
//	quantity<si::mega_scaled<atomic::electron_volt_unit>::unit> one_ev(1.*atomic::eV); //does not compile

	cout<< one_fs<<" "<<one_nanohartree<<" "<< one_Mev <<std::endl;
	cout<< 3.*angstrom<<endl;
	
	quantity<atomic::volume> vv=123.*pow<3>(atomic::a_0);	
	cout<< vv<<endl;
	//using boost::mpl::times; using boost::mpl::divides;
	//cout<< quantity<times<angstrom_unit, times<angstrom_unit, angstrom_unit>::type>::type>(vv)<<endl;
	//auto vvv = pow<3>(3.*angstrom);
	
//	cout<< pow<3>(3.*angstrom) <<" "<<(typeof(pow<3>(1.*angstrom)))(vvv)<<endl;
	cout<< angstrom_unit()<<endl;
	cout<<simplify_typename(pow<3>(3.*angstrom))<<endl;
	quantity<angstrom_unit> L(1.*angstrom);
	typedef	unit<volume_dimension, make_system<angstrom_base_unit>::type> a3_u;

	quantity<a3_u
	> V = pow<3>(L);
	clog<<quantity<si::volume>(V)<<endl;
	
	clog<<simplify_typename(2.*angstrom)<<endl;
	//power_typeof_helper<angstrom_base_unit, static_rational<3,1> >::type ggg ;//=pow<3>(3.*angstrom);
	clog<< fabs(-1.*angstrom) <<endl;
	
	typedef unit<volume_dimension, make_system<angstrom_base_unit>::type> ang3_unit;
	
//	quantity<ang3_unit> vV = pow<3>(3.*angstrom);
/*	quantity<power_typeof_helper<
		quantity<angstrom_unit>::unit_type, //angstrom_unit, 
		static_rational<3> 
	>::type> */
	power_typeof_helper<
		quantity<angstrom_unit>, 
		static_rational<3>
	>::type
	vV =  pow<3>(3.*angstrom);
	clog<<vV<<endl;
	
	quantity<power_typeof_helper<
		atomic::energy, //angstrom_unit, 
		static_rational<2> 
	>::type> eE =  pow<2>(3.*atomic::hartree);
	clog<<eE<<endl;
	typedef typeof(eE) ddd;
	clog<<eE.value()<<"."<<ddd::unit_type()<<std::endl;
    
   	quantity<atomic::energy> a686(1.*nonatomic::rydberg);
   	quantity<nonatomic::rydberg_unit> aRy(a686);
   	std::cout<< a686 << nonatomic::rydberg_unit() << std::endl;
	cout<< quantity<atomic::dimensionless>((1.*nonsi::electron_volt)/(atomic::hartree)) <<endl;
   
	std::cout<< boost::units::value(3.)<<std::endl;
	
	std::cout << pow<2>(atomic::e)/pow<2>(atomic::bohr)*atomic::k_1 <<std::endl;
	{
//		quantity<si::energy> esi(1.e-17*si::joules);
		quantity<atomic::time> tat(2.*atomic::time_unit);
		quantity<si::time> tsi(tat);
		quantity<nonsi::angstrom_unit> l(3.*nonsi::angstrom);
		cout << tsi << " "<< tat << " "<<l<<endl;
	}
	{
		cout << "ddd "<<10.811*nonsi::daltons << endl;
		quantity<si::mass> mb(10.811*nonsi::Da);
		quantity<si::mass_density> mbcgs(mb/(si::meter*si::meter*si::meter));
		cout << mbcgs << endl;
	}
//	{
//		clog << name<boost::units::energy_dimension>() << endl;
//	}
	{
		quantity<divide_typeof_helper<si::mass, si::length>::type> ld;
		clog << "ld = "<< ld <<endl;
	}
	{
		clog << "test heat capacity, atomic/si" << endl;
		quantity<atomic::heat_capacity> cv = 3.*atomic::k_B;
		//quantity<si::heat_capacity> cvsi; 
		quantity<boltzman_constant_si_unit> cvsi(cv);
		quantity<si::energy> en(cv*3.*si::kelvin);
		clog 
			<< "cv = "<<cv <<endl
			<< "cvsi = "<<cvsi << endl  
			<< "en = "<<en <<endl
			//<< cv*4.*si::kelvin << endl
			//<< atomic::energy(cv*4.*si::kelvin) << endl
		;
		clog << std::flush;
		clog << quantity<si::temperature>(1.*atomic::temperature()) << endl;
		clog << quantity<atomic::energy>(1.*nonsi::electron_volt) <<endl;
		clog << quantity<atomic::temperature>(1.*si::kelvin) << endl;
		clog << quantity<atomic::heat_capacity>((1.*nonsi::electron_volt)/quantity<atomic::temperature>(1.*si::kelvin)) << endl;
		quantity<si::temperature> t = 200.*si::kelvin;
		clog << quantity<si::energy>(t*cv) << endl;
	}
	{
		clog << "... " << quantity<atomic::length>(5.*nonsi::angstrom) << endl;
		clog << "... " << quantity<atomic::force>(10.*atomic::hartree/atomic::bohr) << endl;
	}
	{
		clog << (1.*atomic::k_B)*(1000.*si::kelvin) << endl;
		clog << quantity<si::heat_capacity>(1.*atomic::k_B) << endl;
		{
			//using namespace auto_right_conversion_operators;
			using namespace auto_conversion_operators;
			clog << (1.*atomic::k_B)*(1000.*si::kelvin)*si::dimensionless() << endl;
		}
		{
			using namespace auto_conversion_operators::to_right;
			cout << 1.*nonsi::angstrom + 1.*atomic::bohr << endl; 
			//clog << auto_conversion_operators::to_left::operator+(1.*atomic::bohr, 1.*nonsi::angstrom) << std::endl;
		}
	}
	clog << "end test" << endl;
	double atob(quantity<atomic::length>(1.*nonsi::angstrom)/atomic::bohr);
	clog << "atob " << atob << std::endl;

	cout << 3.*si::meter/si::femtosecond << endl;
	cout << 3.*si::meter/(si::femto*si::second) << endl;

	/*
	typedef typeof(si::femto*si::second)    femtosecond;
	cout << 3.*atomic::bohr/femtosecond() << endl; 
	cout << "bohr/fs" << 2.*atomic::bohr/(si::femto_scaled<si::time>::unit()) << endl;
	cout << 3.*cgs::second/cgs::centimeter << endl;
	cout << 3.*femtosecond_base_unit() << endl;
	*/
	quantity<si::temperature> c;
	cout << "******************" << c << std::endl;
	quantity<atomic::temperature> ca ; 
	cout << "******************" << ca << std::endl;
	quantity<atomic::temperature> ca2(c);
	cout << "******************" << ca2 << std::endl;
	return 0;
}
#endif
#endif
// ~/usr/bin-hera/astyle --brackets=attach --indent=tab --indent-col1-comments --pad-oper --delete-empty-lines --add-brackets --keep-one-line-statements --convert-tabs --align-pointer=type --unpad-paren  qbox.hpp
// Editor modelines  -  http://www.wireshark.org/tools/modelines.html
// Local variables:
// c-basic-offset: 4
// tab-width: 4
// indent-tabs-mode: t
// truncate-lines: 1
// End:
/* vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent: */

