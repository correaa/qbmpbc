#ifdef COMPILE_INSTRUCTIONS
ln -sf $0 .$0.cpp && c++ -Wfatal-errors `#-std=c++0x` .$0.cpp -o .$0.x -D_BOOST_UNITS_SYSTEMS_ATOMIC_SI_CONVERSION_TEST -I../../../.. && ./.$0.x $@; exit;
#endif
#ifndef BOOST_UNITS_SYSTEMS_ATOMIC_SI_CONVERSION
#define BOOST_UNITS_SYSTEMS_ATOMIC_SI_CONVERSION
//#include"boost/units/systems/atomic.hpp"
#include<boost/units/systems/si/action.hpp>
#include<boost/units/systems/si/codata/atomic-nuclear_constants.hpp>
#include<boost/units/systems/si/codata/electromagnetic_constants.hpp> // for codata::e
#include<boost/units/systems/si/codata/physico-chemical_constants.hpp> //for m_u (dalton unit, http://en.wikipedia.org/wiki/Atomic_mass_unit) and k_B (boltzman constant)
#include<boost/units/physical_dimensions/heat_capacity.hpp>
#include "../../systems/atomic.hpp"

//conversion factors using library codata
/* 1 - mass */
BOOST_UNITS_DEFINE_CONVERSION_FACTOR(
	boost::units::atomic::electron_mass_base_unit,
	boost::units::si::mass,
	double, boost::units::si::constants::codata::m_e.value().value() //9.10938215e-31
);
BOOST_UNITS_DEFAULT_CONVERSION(
	boost::units::atomic::electron_mass_base_unit,
	boost::units::si::mass
);

/* 2 - electric charge */
typedef boost::units::unit<
		boost::units::electric_charge_dimension, 
		boost::units::si::system
	> electrinc_charge_si_unit;
BOOST_UNITS_DEFINE_CONVERSION_FACTOR(
	boost::units::atomic::elementary_charge_base_unit,
	electrinc_charge_si_unit,
	double, boost::units::si::constants::codata::e.value().value() //1.60217653e-19 /*http://en.wikipedia.org/wiki/Atomic_units#Fundamental_atomic_units*/
);
BOOST_UNITS_DEFAULT_CONVERSION(
	boost::units::atomic::elementary_charge_base_unit,
	electrinc_charge_si_unit
);

/* 3 - action */
typedef boost::units::si::action si_action;
BOOST_UNITS_DEFINE_CONVERSION_FACTOR(
	boost::units::atomic::reduced_planck_constant_base_unit,
	si_action,
	double, 1.05457168e-34 /*http://en.wikipedia.org/wiki/Atomic_units#Fundamental_atomic_units*/
);
BOOST_UNITS_DEFAULT_CONVERSION(
	boost::units::atomic::reduced_planck_constant_base_unit,
        si_action
);

/* 4 - electric constant */
typedef boost::units::unit<
		boost::units::derived_dimension<
			boost::units::length_base_dimension,   3,
			boost::units::mass_base_dimension,     1,
			boost::units::time_base_dimension,    -4,
			boost::units::current_base_dimension, -2
		>::type,
		boost::units::si::system
	> electric_constant_si_unit;
BOOST_UNITS_DEFINE_CONVERSION_FACTOR(
	boost::units::atomic::coulomb_force_constant_base_unit,
	electric_constant_si_unit,
	double, 8.9875517873681e9 /*http://en.wikipedia.org/wiki/Atomic_units#Derived_atomic_units*/
);
BOOST_UNITS_DEFAULT_CONVERSION(
	boost::units::atomic::coulomb_force_constant_base_unit,
	electric_constant_si_unit
);
/* 5 - boltzman_constant */
typedef 
	boost::units::unit<
		boost::units::heat_capacity_dimension, 
		boost::units::si::system
	> 
	boltzman_constant_si_unit
;
BOOST_UNITS_DEFINE_CONVERSION_FACTOR(
	boost::units::atomic::boltzman_constant_base_unit,
	boltzman_constant_si_unit,
	double, 1.3806504e-23 /*http://en.wikipedia.org/wiki/Atomic_units#Derived_atomic_units*/
);
BOOST_UNITS_DEFAULT_CONVERSION(
	boost::units::atomic::boltzman_constant_base_unit,
	boltzman_constant_si_unit
);
/*
BOOST_UNITS_DEFINE_CONVERSION_FACTOR(
	boost::units::si::temperature,
	boost::units::atomic::temperature,
	double, 22.
);
BOOST_UNITS_DEFAULT_CONVERSION(
	boost::units::si::temperature,
	boost::units::atomic::temperature
);*/

/*
BOOST_UNITS_DEFINE_CONVERSION_FACTOR(
	boost::units::atomic::energy,
	boost::units::si::energy,
	double, si::constants::codata::E_h.value().value() //4.35974394e-18
);
BOOST_UNITS_DEFAULT_CONVERSION(
	boost::units::atomic::energy,
	boost::units::si::energy
);*/


/* commented because hartree is not a base unit
BOOST_UNITS_DEFINE_CONVERSION_FACTOR(
//	atomic::hartree_energy_base_unit,
	atomic::energy_unit,
	si::energy,
	double, si::constants::codata::E_h.value().value() //4.35974394e-18
);
BOOST_UNITS_DEFAULT_CONVERSION(
//	atomic::hartree_energy_base_unit,
	atomic::energy_unit,
	si::energy
);
*/
/*
BOOST_UNITS_DEFINE_CONVERSION_FACTOR(
	atomic::boltzman_constant_base_unit,
	(unit<
		derived_dimension<
			length_base_dimension, 2,
			mass_base_dimension, 1,
			time_base_dimension, -2, // == ,energy_dimension, 1,
			temperature_base_dimension, -1
		>::type,
		si::system
	>),
	double, si::constants::codata::k_B.value().value() //1.3806504e-23
);
BOOST_UNITS_DEFAULT_CONVERSION(
	atomic::boltzman_constant_base_unit,
	(unit<
		derived_dimension<
			length_base_dimension, 2,
			mass_base_dimension, 1,
			time_base_dimension, -2, // == ,energy_dimension, 1,
			temperature_base_dimension, -1
		>::type,
		si::system
	>)
);
*/

//#endif
#endif
#ifdef _BOOST_UNITS_SYSTEMS_ATOMIC_SI_CONVERSION_TEST
#include<iostream>
using namespace std;
using namespace boost::units;
int main(){
	cout << "hola" <<endl;
	cout << quantity<boltzman_constant_si_unit>() << " " << boltzman_constant_si_unit() <<endl;
	cout << quantity<si::energy>() << " "<< 1.*boltzman_constant_si_unit()*si::kelvin << endl;
	cout << quantity<atomic::heat_capacity>(1.*boltzman_constant_si_unit()) << endl;
//	cout << quantity<si::temperature>(atomic::temperature()) << endl;
	return 0;
}
#endif

