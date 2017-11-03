#ifdef compile_instructions
ln -sf $0 .$0.cpp && c++ -I$HOME/prj .$0.cpp `#-Wall` `#-Wfatal-errors` -L$HOME/lib `pkg-config --libs gsl` -D_TEST_BOOST_UNITS_PHOENIX_HPP -o ./.$0.x && ./.$0.x $@; exit
#endif
#ifndef _BOOST_UNITS_PHOENIX_HPP
#define _BOOST_UNITS_PHOENIX_HPP
#include<boost/spirit/home/phoenix.hpp>
#include<boost/units/unit.hpp>

// phoenix2 units bridge, result_of_plus/minus don't seem to be necesary, todo: *=, /=, an overload operator*(unit,actor<T>) to resolve ambiguity
//namespace boost{ 
//namespace phoenix{
//}}

namespace boost{
namespace phoenix{
using namespace boost::units;
/*
template<class X, class YUnit, class TY>
struct result_of_multiplies<X&, quantity<YUnit, TY>& >{
	typedef typename multiply_typeof_helper<X, quantity<YUnit, TY> >::type type;
};

template<class XUnit, class YUnit, class TX, class TY>
struct result_of_multiplies<quantity<XUnit, TX>&, quantity<YUnit,TY>& >{
	typedef typename multiply_typeof_helper<quantity<XUnit, TX>, quantity<YUnit, TY> >::type type;
};
*/
#define RESULT_OF_QUANTITY_GEN( PhoenixnamE, UnitsnamE, RefQ1, RefQ2 ) \
	template<class XUnit, class YUnit, class TX, class TY> \
	struct result_of_##PhoenixnamE <quantity<XUnit,TX> RefQ1, quantity<YUnit,TY> RefQ2>{ \
		typedef typename UnitsnamE##_typeof_helper<quantity<XUnit,TX>&,quantity<YUnit,TY>& >::type type; \
	};
	RESULT_OF_QUANTITY_GEN(multiplies, multiply, &,  )
	RESULT_OF_QUANTITY_GEN(multiplies, multiply, &, &)
	RESULT_OF_QUANTITY_GEN(multiplies, multiply,  ,  )
	RESULT_OF_QUANTITY_GEN(multiplies, multiply,  , &)
	RESULT_OF_QUANTITY_GEN(divides   , divide  , &,  )
	RESULT_OF_QUANTITY_GEN(divides   , divide  , &, &)
	RESULT_OF_QUANTITY_GEN(divides   , divide  ,  ,  )
	RESULT_OF_QUANTITY_GEN(divides   , divide  ,  , &)
#undef RESULT_OF_QUANTITY_GEN

/*
#define RESULT_OF_QUANTITY_UNITS_GEN( PhoenixopnamE, UnitsopnamE ) \
template< \
	typename X, \
	class Dim, class System> \
struct result_of_##PhoenixopnamE < \
	X&, \
	boost::units::unit<Dim, System> \
>{ \
	typedef \
		typename boost::units::UnitsopnamE##_typeof_helper< \
			X, \
			boost::units::unit<Dim, System> \
		>::type type; \
};
	RESULT_OF_QUANTITY_UNITS_GEN(multiplies, multiply)
	RESULT_OF_QUANTITY_UNITS_GEN(divides   , divide  )
#undef RESULT_OF_QUANTITY_UNITS_GEN
*/
// this refines the operator* so the phoenix::operator* is used instead of units::operator*
#define REFINE_ARITHMETIC_OPERATORS(PhoenixopnamE, CsymboL) \
template< \
	class PhoenixExpr, \
	class System, \
    class Dim \
> \
actor< \
	composite< \
		PhoenixopnamE##_eval, \
		boost::fusion::vector< \
			PhoenixExpr, \
			value< \
				boost::units::unit<Dim, System> \
			> \
		> \
	> \
> \
operator CsymboL (const actor<PhoenixExpr>& lhs, const unit<Dim,System>& rhs){ \
    return compose<PhoenixopnamE##_eval>(lhs, rhs); \
}
//REFINE_ARITHMETIC_OPERATORS(multiplies, *)
//REFINE_ARITHMETIC_OPERATORS(divides, /)

#undef REFINE_ARITHMETIC_OPERATORS
}}

#include <boost/units/pow.hpp>
namespace boost{
namespace phoenix{
	template<long N>
	struct pow_eval{ 
		template <typename Env, typename Arg_> struct result{
			typedef typename power_typeof_helper<typename boost::mpl::at_c<typename Env::args_type, 0>::type, static_rational<N> >::type type;
			//typedef typename boost::mpl::at_c<typename Env::args_type, 0>::type type;
		};
		template <typename RT, typename Env, typename Arg_> static RT eval(Env const& env, Arg_ const& arg_){
			return pow<N>(arg_.eval(env));
		}
	};
	template<long N, typename Arg_> \
	typename boost::enable_if_c<is_actor<Arg_>::value, \
		actor<typename as_composite<pow_eval<N>, Arg_>::type> \
	>::type \
	pow(Arg_ const& arg_){ \
		return compose<pow_eval<N> >(arg_); \
	}

}}

#endif

#ifdef _TEST_BOOST_UNITS_PHOENIX_HPP
#include<typeinfo>
#include<boost/units/systems/si.hpp>

using namespace boost::phoenix;
using namespace boost::phoenix::arg_names;
using namespace boost::units;

template<class LambdaExp>
void expression_analizer(actor<LambdaExp> l){
	std::clog
		    <<"type of argument expression: "
		    <<typeid(LambdaExp).name()<<std::endl;
	std::clog<<"type of return on double: "
		    <<typeid(
		            typename LambdaExp::template result<
		                    basic_environment<double>
		            >::type
		    ).name()<<std::endl;
	std::clog
		    <<"type of return on quantity<length>: "
		    <<typeid(
		            typename LambdaExp::template result<
		                    basic_environment<
		                            quantity<si::length>
		                    > >::type
		    ).name()<<std::endl;
}
#include<boost/function.hpp>

namespace boost{
namespace phoenix{
template<>
struct result_of_multiplies<quantity<si::length>, quantity<si::length>& >{
	typedef quantity<si::area> type;
};
}}

using namespace std;
int main(){
	cout<<boost::units::detail::demangle(typeid(arg1*arg1).name())<<endl;
	expression_analizer( (arg1+arg1)+arg1
		// /(arg1 - 3.*si::hertz)
	);
	//using boost::phoenix::operator*;
	quantity<si::area> x = pow<2>(5.*si::meter);
	cout << x <<endl;
	using boost::phoenix::pow;
	boost::function<quantity<si::area>(quantity<si::length>)> f = pow<2>(arg1);
	cout << f(2.*si::meter) << endl;
	quantity<si::area> ar = 3.*pow<2>(si::meter);
	boost::function<quantity<si::area>
		(quantity<si::length>)
		//(quantity<si::length> const&) // <-- never do this, it must be (quantity<si::length>)
	> f2 = 
		//arg1*arg1;
		//arg1*si::meter; //todo
		arg1*(1.*si::meter) + ar;
	cout << f2(2.*si::meter) << endl;
	return 0;
}
#endif

