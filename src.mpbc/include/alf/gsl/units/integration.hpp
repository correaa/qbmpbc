#if 0
ln -sf $0 .$0.cpp && c++ -std=c++0x -I$HOME/prj .$0.cpp -Wall `#-Wfatal-errors` -L$HOME/lib `pkg-config --libs gsl` -D_TEST_INTEGRATION_BOOST_UNITS_HPP -o ./.$0.x && ./.$0.x $@; exit
#endif
#ifndef INTEGRATION_BOOST_UNITS_HPP
#define INTEGRATION_BOOST_UNITS_HPP

#include "../../gsl/integration.hpp"
#include "../../gsl/units/adimensional.hpp"
#include<boost/units/quantity.hpp>
#include "boost/units/interval.hpp"
#include "boost/units/phoenix.hpp"
#include<boost/type_traits/remove_reference.hpp>
#include<iostream> //debug
#include<boost/units/detail/utility.hpp> //debug
#include<boost/units/systems/si.hpp>//debug
#include<boost/units/lambda.hpp>//if lots of error messages about lambda appear is because we need to make lambda play nicely with units!!
#include<boost/lambda/bind.hpp>//this ensures that lambda::bind is used!!
#ifndef BOOST_LAMBDA_DEBUG_HPP
#define BOOST_LAMBDA_DEBUG_HPP
#include<boost/lexical_cast.hpp>
#include<boost/just.hpp>
namespace boost{
namespace lambda{
using namespace boost::tuples;

template<int N>
std::string to_string(const lambda_functor<boost::lambda::placeholder<N> >&){
	return "_"+boost::lexical_cast<std::string>(N);
}

//template<class T>
//std::string to_string(const T& t){
//	return boost::lexical_cast<std::string>(t);
//}

#define BOOST_LAMBDA_TO_STRING_GEN(action_name, symbol_string)\
template<class LExp1, class LExp2> \
std::string to_string(const \
	lambda_functor<lambda_functor_base<arithmetic_action<action_name>, tuple< \
		LExp1, \
		LExp2  \
	> > >& action_name##_expr){ \
	return to_string(get<0>(action_name##_expr.args)) + #symbol_string +to_string(get<1>(action_name##_expr.args)); \
}
BOOST_LAMBDA_TO_STRING_GEN(plus_action, +);
BOOST_LAMBDA_TO_STRING_GEN(divide_action, /);
BOOST_LAMBDA_TO_STRING_GEN(minus_action, -);
BOOST_LAMBDA_TO_STRING_GEN(multiply_action, *);

template<class FunctionType, class LExpr>
std::string to_string(const 
	lambda_functor<lambda_functor_base<action<2, function_action<2> >, 
		tuple<
			FunctionType, //double (* const)(double), 
			LExpr
		>
	> >& bind_expr){
	return std::string("bind(")+boost::units::detail::demangle(typeid(FunctionType).name()) +", "+to_string(get<1>(bind_expr.args)) + ")";
}
}}
#endif

namespace gsl{
namespace integration{

using std::clog; using std::endl;
using namespace boost::units;
using namespace boost::tuples;
using boost::remove_reference;

/// QAG integration for functions with units

//struct qag_{
template<class UnitIntegrand, class UnitIntegrandDomain>
#define RET_TYPE typename multiply_typeof_helper<quantity<UnitIntegrand>, quantity<UnitIntegrandDomain> >::type //decltype(quantity<UnitIntegrand>()*quantity<UnitIntegrandDomain>()) //needs c++0x
RET_TYPE //operator()
qag (
	boost::function<quantity<UnitIntegrand>(quantity<UnitIntegrandDomain> /*const&*/)> const& f,
	boost::numeric::interval<quantity<UnitIntegrandDomain> > const& iv,
	  epsilon<RET_TYPE> const& eps //epsilon<RET_TYPE>()
	= epsilon<RET_TYPE>(
			std::numeric_limits<RET_TYPE>::epsilon()*100., 
			std::numeric_limits<double>::epsilon()*100.
	)
) //const
{
	boost::units::adimensionalize_function<UnitIntegrandDomain, UnitIntegrand> ad(f);
	double ret = qag(ad, boost::numeric::interval<double>(lower(iv).value(), upper(iv).value()), gsl::epsilon(eps.absolute().value(), eps.relative() ));
	typedef RET_TYPE ret_type;
	return ret_type::from_value(ret);
}
#undef RET_TYPE
using namespace boost::phoenix;

template<class Expression, class UnitIntegrandDomain>
#define RET_TYPE typename multiply_typeof_helper<typename boost::remove_reference<typename Expression::template result<basic_environment<quantity<UnitIntegrandDomain> > >::type>::type, quantity<UnitIntegrandDomain> >::type
RET_TYPE qag(
	boost::phoenix::actor<Expression> const& f_expr,
	boost::numeric::interval<quantity<UnitIntegrandDomain> > const& iv
){
	boost::phoenix::function<adimensionalized_<UnitIntegrandDomain> > const adimensionalized = adimensionalized_<UnitIntegrandDomain>();
	boost::function<double(double)> fad = adimensionalized(boost::phoenix::arg_names::arg1, boost::phoenix::lambda[
		f_expr //e.g boost::phoenix::arg_names::arg1
	]);
	typedef RET_TYPE ret_type;
	double ret = gsl::integration::qag(
		fad,
		boost::numeric::interval<double>( lower(iv).value() , upper(iv).value() )
	);
	return ret_type::from_value(ret) ;
}
#undef RET_TYPE
// QAWC integration for functions with units
template<class LambdaExp, class UnitIntegrandDomain>
#define RET_TYPE typename remove_reference<typename LambdaExp::template sig<tuple<quantity<UnitIntegrandDomain> > >::type>::type
RET_TYPE //lambda version doesn't work with bind in the expression template
qawc(
	lambda_functor<lambda_functor_base<
		arithmetic_action<divide_action>, 
			tuple<
				lambda_functor<//
					LambdaExp
				>,//
				lambda_functor<lambda_functor_base<
						arithmetic_action<minus_action>, 
						tuple<
							lambda_functor<
								placeholder<1>
							>, 
							quantity<UnitIntegrandDomain> const
						> 
				> >
			> 
	> > const& 
	f_expr, 
	boost::numeric::interval<quantity<UnitIntegrandDomain> > const&
	iv
){
	typedef RET_TYPE ret_type;
	double const c = get<1>(get<1>(f_expr.args).args).value();
	double ret = 
	gsl::integration::qawc(
		boost::lambda::bind(
			&ret_type::value,
			get<0>(f_expr.args)(boost::lambda::bind(&quantity<UnitIntegrandDomain>::from_value, _1))
		),
		boost::numeric::interval<double>( lower(iv).value() , upper(iv).value() ),
		c
	);
	return ret_type::from_value(ret);
}
#undef RET_TYPE


template<class LambdaExp, class UnitIntegrandDomain>
#define RET_TYPE typename remove_reference<typename LambdaExp::template result<basic_environment<quantity<UnitIntegrandDomain> > >::type>::type
RET_TYPE qawc(
	actor<composite<
		divides_eval, vector<
			LambdaExp, composite<
			minus_eval, vector<
				argument<0>,
				boost::phoenix::value<
					quantity<UnitIntegrandDomain> 
				> 
			> 
			> 
	> > > 
		f_expr,
	boost::numeric::interval<boost::units::quantity<UnitIntegrandDomain> > const& iv
){
	typedef RET_TYPE ret_type;
	double c = at_c<1>(at_c<1>(f_expr)).val.value();	
	boost::phoenix::function<adimensionalized_<UnitIntegrandDomain> > const adimensionalized = adimensionalized_<UnitIntegrandDomain>();
	boost::function<double(double)> fad = adimensionalized(boost::phoenix::arg_names::arg1, boost::phoenix::lambda[
		actor<LambdaExp>(at_c<0>(f_expr)) //e.g boost::phoenix::arg_names::arg1
	]);
	double ret = gsl::integration::qawc(
		fad,
		boost::numeric::interval<double>( lower(iv).value() , upper(iv).value() ),
		c
	);
	return ret_type::from_value(ret);
}
#undef RET_TYPE

}}

#endif
#ifdef _TEST_INTEGRATION_BOOST_UNITS_HPP
#include<iostream>
#include<boost/units/systems/si.hpp>
#include<boost/units/systems/si/io.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/lambda/bind.hpp>
#include "boost/numeric/make_interval.hpp"
using namespace std;
using namespace boost::units;

double vl(quantity<si::frequency> const& d){return d.value();}

using namespace std;
double g(double const& x){
	return x*x;
}
quantity<si::area> G(quantity<si::length> const& l){
	return l*l;
}

int main(){
	quantity<si::frequency> wp = 3.*si::hertz;
	{
		using namespace boost::lambda;
		cout<<"lambda "
			<< gsl::integration::qawc((_1+1.)/(_1 - 1.), gsl::integration::interval(0.,2.)) << endl;
		cout<<gsl::integration::qawc((_1)/(_1 - wp), gsl::integration::interval(1.*si::hertz, 4.*si::hertz)) 
			<<endl;
	}
	{
		using namespace boost::phoenix;
		using namespace boost::phoenix::arg_names;
		cout<<"phoenix "<< gsl::integration::qawc((arg1)/(arg1 - 3.), gsl::integration::interval(2.,4.)) << endl;
		using namespace gsl::integration;
		quantity<si::frequency> wp = 3.*si::hertz;
		//cout<<"ph: " << 
		qawc( ((arg1*arg1/arg1 + arg1)*arg1)/(arg1 - 3.*si::hertz), gsl::integration::interval(2.*si::hertz, 4.*si::hertz))
		//<<endl
		;
	}
	{
		gsl::function f(&g);
		gsl::result rqag = gsl::integration::qag(f, gsl::integration::interval(0.,1.) );
		cout << rqag << endl;
	}
	{
		boost::function<quantity<si::area>(quantity<si::length>)> F(&G);
		//dont: boost::function<quantity<si::area>(quantity<si::length> const&)> F(&G);
		//adimensionalize_function<si::length, si::area > AF(F);
		//cout << AF(3.) <<endl;
		cout << gsl::integration::qag(
			F, 
			boost::numeric::interval<quantity<si::length> >(0.*si::meter, 1.*si::meter), 
			gsl::integration::epsilon<quantity<si::volume> >(1.e-4)
//			gsl::integration::epsilon<quantity<si::volume> >::relative(1.e-4) 
		) <<endl;
	}
	{
		using namespace boost::phoenix;
		using namespace boost::phoenix::arg_names;
		std::cout << gsl::integration::qag(
			arg1/(2.*si::meter),
			boost::numeric::make_interval(0.*si::meter, 1.*si::meter)
		) << std::endl;
	}
	return 0;
}
#endif
/* vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent: */

