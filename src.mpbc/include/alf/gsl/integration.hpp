#if compilation_instructions
ln -sf $0 .$0.cpp && c++ .$0.cpp `#-Wfatal-errors` -I$HOME/prj -L$HOME/lib `pkg-config --libs gsl` -D_TEST_INTEGRATION_HPP -o ./.$0.x && ./.$0.x $@; exit
#endif
#ifndef INTEGRATION_HPP
#define INTEGRATION_HPP

#include <gsl/gsl_integration.h>

#include "../gsl/function.hpp"
#include "../gsl/result.hpp"
#include "../gsl/epsilon.hpp"
#include<boost/numeric/interval.hpp>
#include<boost/lambda/lambda.hpp> //for smart qawc
#include<boost/spirit/home/phoenix.hpp>//for very smart qawc
#include "boost/just.hpp"
#include"boost/phoenix/cmath.hpp"
#include<iostream>
using std::clog; using std::endl;

namespace gsl{
/// Quadrature integration of a function in one dimension http://www.gnu.org/software/gsl/manual/html_node/Numerical-Integration.html
namespace integration{

	template<class T> class epsilon;
	class relative_error;
	typedef boost::numeric::interval<double> interval_type; //or just interval?
	using gsl::result;
	class workspace;
	class qawo_table;

	struct relative_error : boost::just<double>{
		relative_error(double const& d) : boost::just<double>(d){if(d<0) throw std::domain_error("relative error cannot be negative");}
		template<class T> operator epsilon<T>();
	};
	/// Integration rule, one of id
	struct rule{
		enum id{
			gauss15 = GSL_INTEG_GAUSS15, //  (key = 1)
			gauss21 = GSL_INTEG_GAUSS21, //  (key = 2)
			gauss31 = GSL_INTEG_GAUSS31, //  (key = 3)
			gauss41 = GSL_INTEG_GAUSS41, //  (key = 4)
			gauss51 = GSL_INTEG_GAUSS51, //  (key = 5)
			gauss61 = GSL_INTEG_GAUSS61  //  (key = 6)
		};
	};
	result qng (gsl::function integrand, interval_type const& domain, gsl::epsilon const& eps = gsl::epsilon(std::numeric_limits<double>::epsilon()*100., std::numeric_limits<double>::epsilon()*100.));
	/// Simple adaptive integration with settings http://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html
	result qag (gsl::function integrand, interval_type const& domain, gsl::epsilon const& tolerance, size_t evaluation_limit, rule::id const& r, workspace& w);
	//integrates f(x)/(x-c) principal value, c=pole
	result qawc(gsl::function integrand, interval_type const& domain, double pole, gsl::epsilon const& tolerance, size_t evaluation_limit, workspace& w);
	result qawc(gsl::function integrand, interval_type const& domain, double pole, gsl::epsilon const& tolerance = gsl::epsilon(std::numeric_limits<double>::epsilon()*10000., std::numeric_limits<double>::epsilon()*10000.), size_t evaluation_limit = 100);
    //Boost.Phoenix qawc declaration goes here
	// Adaptive integration for Fourier integrands
	//result qawf(gsl::function integrand, double lower_limit, double tolerance, size_t evaluation_limit, workspace& w, workspace& cicle, qawo_table& table);


	template<class T>
	class epsilon : 
		boost::just<T>::type,
		relative_error
	{
		public:
		epsilon(relative_error const& re) : boost::just<T>::type(0), relative_error(re){}
		epsilon(
			typename boost::just<T>::type absolute_ = std::numeric_limits<T>::epsilon(), 
			relative_error relative_ = std::numeric_limits<double>::epsilon()
		) : 
			boost::just<T>::type(absolute_), 
			relative_error(relative_)
		{
			if(absolute_< T(0)) throw std::domain_error("absolute error cannot be negative");
		}
		typename boost::just<T>::type const& absolute() const{return *this;}
		relative_error const& relative() const{return *this;}
		static epsilon<T> relative(double relative_){return epsilon<T>(0, relative_error(relative_));}
		static epsilon<T> absolute(T absolute_){return epsilon<T>(absolute_, 0);}
	};
	template<class T>
	boost::numeric::interval<T> interval(T const& a, T const& b){
		return boost::numeric::interval<T>(a,b);
	}
	/// Non-adaptive Gauss-Kronrod integration http://www.gnu.org/software/gsl/manual/html_node/QNG-non_002dadaptive-Gauss_002dKronrod-integration.html
	result qng(
		gsl::function f, 
		interval_type const& iv, 
		gsl::epsilon const& eps /*= gsl::epsilon(std::numeric_limits<double>::epsilon()*100., std::numeric_limits<double>::epsilon()*100.)*/
	){
		double ret;
		double abserr;
		size_t neval;
		gsl_integration_qng(&f, lower(iv), upper(iv), eps.absolute(), eps.relative(), &ret, &abserr, &neval);
		clog<<neval<<endl;
		return gsl::result(ret, abserr);
	}
	class workspace{
		gsl_integration_workspace* pimpl_;
		size_t size_;
		public:
		workspace(size_t n) : pimpl_(gsl_integration_workspace_alloc(n)), size_(n){}
		~workspace(){gsl_integration_workspace_free(pimpl_);}
		operator gsl_integration_workspace*(){return pimpl_;}
		size_t const& size() const{return size_;}
	};
	struct qawo_function{
		enum id{
			cosine = GSL_INTEG_COSINE,
			sine = GSL_INTEG_SINE
		};
	};
	class qawo_table{
		gsl_integration_qawo_table* pimpl_;
		public:
		qawo_table(double w, double L, qawo_function::id const& sineorcosine, size_t n) : 
			pimpl_(gsl_integration_qawo_table_alloc(w, L, (gsl_integration_qawo_enum)sineorcosine, n)){}
		void set(double omega, double length, qawo_function::id const& sineorcosine){assert(
			gsl_integration_qawo_table_set(pimpl_, omega, length, (gsl_integration_qawo_enum)sineorcosine)
		);}
		void set_length(double length){assert(
			gsl_integration_qawo_table_set_length(pimpl_, length)
		);}
		operator gsl_integration_qawo_table*(){return pimpl_;}
		~qawo_table(){gsl_integration_qawo_table_free(pimpl_);}
	};
	result qawo(
		gsl::function f, 
		double a, //b = a+L (given by workspace
		gsl::epsilon const& eps, 
		const size_t limit, \
		workspace& w, 
		qawo_table& wf
	){
		double ret;
		double abserr=-1;
		gsl_integration_qawo(&f, a, eps.absolute(), eps.relative(), limit, w, wf, &ret, &abserr);
		return result(ret, abserr);
	}
	result qawo(
		gsl::function f, 
		interval_type const& ab, 
		qawo_function::id const& sineorcosine,
		double omega,
		gsl::epsilon const& eps = gsl::epsilon(
			std::numeric_limits<double>::epsilon()*340.65735903, /*50 log(2)*/
			std::numeric_limits<double>::epsilon()*340.65735903
		), 
		const size_t limit = 100, 
		const size_t levels = 10 //number of levels of coefficients that are computed
	){
		qawo_table wf(omega, width(ab), sineorcosine, levels);
		workspace w(limit);
		return qawo(f, ab.lower(), eps, limit, w, wf);
	}
	result qag(
		gsl::function f, 
		interval_type const& iv,
		gsl::epsilon const& eps, 
		size_t limit,
		rule::id const& r,
		workspace& w
	){
		assert(limit <= w.size() );
		double ret;
		double abserr;
		//manual: http://www.gnu.org/software/gsl/manual/html_node/QAG-adaptive-integration.html
		//	int gsl_integration_qag(const gsl_function * f, double a , double b , double epsabs , double epsrel , size_t limit, int key, gsl_integration_workspace * workspace, double * result, double * abserr)
		gsl_integration_qag(&f                    , lower(iv), upper(iv), eps.absolute(), eps.relative(), limit       , r      , w                                    , &ret           , &abserr        );
		return gsl::result(ret, abserr);
	}
	/// Simple adaptive integration with default settings
	gsl::result qag(
		gsl::function f, 
		interval_type const& iv,
		gsl::epsilon const& eps = gsl::epsilon(
			std::numeric_limits<double>::epsilon()*1000., 
			std::numeric_limits<double>::epsilon()*1000.
		),
		size_t limit = 100, //infinite? no, it occupies memory
		rule::id const& r = rule::gauss51
	){
		workspace w(limit);
		return qag(f, iv, eps, limit, r, w);
	}
	gsl::result qawc(
		gsl::function f, 
		interval_type const& iv,
		double c, //integrates f(x)/(x-c) principal value
		gsl::epsilon const& eps, 
		size_t limit, 
		workspace& w){
		double ret;
		double abserr;
		// manual: http://www.gnu.org/software/gsl/manual/html_node/QAWC-adaptive-integration-for-Cauchy-principal-values.html
//	int gsl_integration_qawc(gsl_function * f, double a , double b , double c, double epsabs , double epsrel , size_t limit, gsl_integration_workspace * workspace, double * result, double * abserr)
		gsl_integration_qawc(&f              , lower(iv), upper(iv), c       , eps.absolute(), eps.relative(), limit       , w                                    , &ret           , &abserr        );
		return gsl::result(ret, abserr);
	}
	result qawc(
		gsl::function f, 
		interval_type const& iv,
		double c, //integrates f(x)/(x-c) principal value, or f(x)/x if c=0
		gsl::epsilon const& eps /*= gsl::epsilon(
			std::numeric_limits<double>::epsilon()*10000., 
			std::numeric_limits<double>::epsilon()*10000.
		)*/,
		size_t limit /*= 100*/){
		workspace w(limit);
		return qawc(f, iv, c, eps, limit, w);
	}
	//smart qawc support for phoenix
	using boost::phoenix::actor; using boost::phoenix::composite; using boost::phoenix::divides_eval; using boost::phoenix::minus_eval; using boost::phoenix::argument; using boost::phoenix::value;
	//using boost:::fusion::vector;
	//using namespace boost::phoenix;
	using namespace boost::fusion;

	/*
	void caca(double x){}
	template<class PhoenixExpression>
	void qawotest(
		actor<composite<
			boost::phoenix::multiplies_eval, vector<
				PhoenixExpression, //eg value<double>, 
				composite<boost::phoenix::sin_eval, vector<
						composite<boost::phoenix::multiplies_eval, vector<
							value<double>, 
							argument<0>
				> > > > 
			> 
		> > f_expr
	){
		//caca(at_c<0>(at_c<1>(f_expr)));
		clog<< actor<value<double> >(at_c<0>(at_c<0>(at_c<1>(f_expr))))() <<endl;
		return;
	} */

	template<class PhoenixExpression>
	result qawo(
		actor<
			composite<
				boost::phoenix::multiplies_eval, vector<
					PhoenixExpression, //eg value<double>, 
					composite<
						boost::phoenix::sin_eval, 
						vector<
							composite<
								boost::phoenix::multiplies_eval, 
								vector<
									value<double>, 
									argument<0>
								> 
							> 
						> 
					> 
				> 
			> 
		> f_expr, 
		interval_type const& domain
	){
		return qawo(
			actor<PhoenixExpression>(at_c<0>(f_expr)),
			domain,
			qawo_function::sine,
			actor<value<double> >(at_c<0>(at_c<0>(at_c<1>(f_expr))))() //or val
			// expression = blabla * sin( 0.1 * arg1); 0.1 is expression[1][0][0]
		);
	}
	
	template<class LambdaExp>
	gsl::result qawc(
		actor<composite<
			divides_eval, vector< 
				LambdaExp, composite< 
				minus_eval, vector<
					argument<0>, 
					value<double>
					> > > > > f_expr,
		interval_type const& iv,
		gsl::epsilon const& eps = gsl::epsilon(
			std::numeric_limits<double>::epsilon()*34.65735903, /*50 log(2)*/
			std::numeric_limits<double>::epsilon()*34.65735903
		),
		size_t limit = 100
	){
		workspace w(limit);
		return qawc(
			actor<LambdaExp>(at_c<0>(f_expr)),
			iv, 
			actor<value<double> >(at_c<1>(at_c<1>(f_expr)))(), //or .val
			eps, limit, w
		);
	}
	using namespace boost::lambda;
	using namespace boost::tuples;
	template<class LambdaExp>
	result qawc(
		//typename qawc_function<LambdaExp>::type const&  //doesn't work in gcc 4.4
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
								double const // <-- attention adds const to avoid unreadable errors
							> 
					> >
				> 
		> >
			f_expr,
		interval_type const& iv,
		gsl::epsilon const& eps = gsl::epsilon(
			std::numeric_limits<double>::epsilon()*34.65735903, /*50 log(2)*/
			std::numeric_limits<double>::epsilon()*34.65735903
		),
		size_t limit = 100){
		workspace w(limit);
		return qawc(
			get<0>(f_expr.args), 
			iv, 
			get<1>(get<1>(f_expr.args).args), 
			eps, 
			limit, 
			w
		);
	}
}
}

#endif

//test

#ifdef _TEST_INTEGRATION_HPP
#include<boost/numeric/interval/io.hpp>
#include<boost/lambda/lambda.hpp>
#include<iomanip>
#include<boost/units/io.hpp>
#include<boost/spirit/home/phoenix.hpp>
#include"boost/phoenix/cmath.hpp" //boost::phoenix::sin, etc

using std::cout;
using std::endl;
using namespace std;
double g(double const& x){
	return x*x*x*x*exp(x)+sin(x);
}
int main(){
	using namespace boost::lambda;
	gsl::function f(&g);
	cout<< std::numeric_limits<double>::epsilon()<<endl;
	double r = gsl::integration::qng(f, gsl::integration::interval(0.,1.) );
	//gsl::integration::interval_type iv = gsl::integration::qng(f, gsl::integration::interval(0.,1.) );
	gsl::result rt = gsl::integration::qng(f, gsl::integration::interval(0.,1.) );
	clog<< setprecision(18) << r << " "<< /*iv <<*/ " "<< rt << endl;
	gsl::result rqag = gsl::integration::qag(f, gsl::integration::interval(0.,1.) );
	clog << rqag <<endl;
	clog << boost::units::simplify_typename(_1/(_1-20.))<<endl;
	clog << gsl::integration::qawc(constant(1.), gsl::integration::interval(-1.,1.),0.)<<endl;

	clog << gsl::integration::qawc(boost::lambda::constant(1.), gsl::integration::interval(0.,2.),1.)<<endl;
	clog << gsl::integration::qawc(
			boost::phoenix::val(1.)/(boost::phoenix::arg_names::arg1-1.), 
			gsl::integration::interval(0.,2.))<<endl;
	clog << gsl::integration::qawc(boost::lambda::constant(1.)/(_1 - 1.), gsl::integration::interval(0.,100000.)) << endl;

	{
		using namespace gsl::integration;
		clog << qawo(f, interval(0.,5.), qawo_function::cosine, 0.1) << endl;
		double dos = 2.0;
		clog << boost::phoenix::sin(0.1*boost::phoenix::arg_names::arg1)(dos) << endl;
		using namespace boost::phoenix; using namespace arg_names;
		double k=0.1;
		clog << qawo(
			arg1*arg1*sin(k*arg1),
			interval(0.,5.)
		) << endl;
	}
	return 0;
}
#endif
// ~/usr/bin-hera/astyle --brackets=attach --indent=tab --indent-col1-comments --pad-oper --delete-empty-lines --add-brackets --keep-one-line-statements --convert-tabs --align-pointer=type --unpad-paren  qbox.hpp
// Editor modelines  -  http://www.wireshark.org/tools/modelines.html
// Local variables:
// c-basic-offset: 4
// tab-width: 4
// indent-tabs-mode: t
// truncate-lines: 1
// End:
/* vim:set ft=cpp ts=4 sw=4 sts=4 cindent: */

