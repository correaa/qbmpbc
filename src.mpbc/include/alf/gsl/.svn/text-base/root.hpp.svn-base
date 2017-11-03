#ifdef COMPILATION_INSTRUCTIONS
ln -sf $0 .$0.cpp && c++ `#-Wfatal-errors` -Wall .$0.cpp -L$HOME/lib `pkg-config --libs gsl` -D_TEST_GSL_ROOT_HPP -o .$0.x && ./.$0.x $@
rm -f $0.x $0.cpp
exit
#endif
#ifndef GSL_ROOT_HPP
#define GSL_ROOT_HPP
#include<gsl/gsl_roots.h>

#include<boost/numeric/interval.hpp>
#include<iostream> //for debug
#include<boost/operators.hpp> //dereferenceable
#include<boost/utility.hpp> //noncopyable
#include"../gsl/function.hpp"
#include"../gsl/function_fdf.hpp"
#include"../gsl/result.hpp"
#include<boost/iterator/iterator_facade.hpp>


namespace gsl{
/// Numerical root finders in one dimension with no derivatives http://www.gnu.org/software/gsl/manual/html_node/One-dimensional-Root_002dFinding.html
namespace root{
	typedef boost::numeric::interval<double> interval;
	/// Generic bracket solver class, http://www.gnu.org/software/gsl/manual/html_node/Root-Bracketing-Algorithms.html
	class solver : 
		public boost::iterator_facade<solver, result const, boost::forward_traversal_tag, result const /*noref*/ >,
		boost::noncopyable
	{
		function f_; //should own it or not?
		gsl_root_fsolver* pimpl_;
		public:
		typedef const gsl_root_fsolver_type* dynamic_type;
		static dynamic_type const bisection; //< Bisection algorithm
		static dynamic_type const false_position; //< False Position/Regula Falsi position
		static dynamic_type const brent; //< Brent-Dekker method
		//solver(type T=bisection) : fpimpl_(gsl_root_fsolver_alloc(T)){} //todo: allow other types of solvers to be constructed
		solver(gsl::function const /*const&*/ f, interval const /*const&*/ bracket) : 
			f_(f),
			pimpl_(gsl_root_fsolver_alloc(brent)){
			set(f_, bracket);
		}
		void set(gsl::function const& f, boost::numeric::interval<double> const& bracket){
			f_=f;
			assert(not f.empty());
			gsl_root_fsolver_set(pimpl_, const_cast<gsl_function*>((const gsl_function*)&f), lower(bracket), upper(bracket));
			return;
		}
		std::string name() const{return gsl_root_fsolver_name(pimpl_);}
		~solver(){gsl_root_fsolver_free(pimpl_);}
		solver& iterate(){
			gsl_root_fsolver_iterate(pimpl_);
			return *this;
		}
		double x_lower() const{return gsl_root_fsolver_x_lower(pimpl_);}
		double x_upper() const{return gsl_root_fsolver_x_upper(pimpl_);}
		interval x() const{return interval(x_lower(), x_upper());}
		double root() const{return gsl_root_fsolver_root(pimpl_);}
		//gsl::result 
		result dereference() const{ //result operator*() const{
			return result(root(), (x_upper()-x_lower())/2.); 
			//return ret;
			//return root();
		}
		void increment() {this->iterate();}//solver& operator++(){return this->iterate();}
	};
	solver::dynamic_type const solver::bisection      = gsl_root_fsolver_bisection; 
	solver::dynamic_type const solver::false_position = gsl_root_fsolver_falsepos, regula_falsi=gsl_root_fsolver_falsepos; 
	solver::dynamic_type const solver::brent          = gsl_root_fsolver_brent; 
#if 1
	/// Generic root finder with derivatives
	class fdfsolver : 
		public boost::iterator_facade<fdfsolver, result const, boost::forward_traversal_tag, result /*noref*/>,
		boost::noncopyable
	{
		function_fdf fdf_;
		gsl_root_fdfsolver* pimpl_;
		double last_;
		public:
		typedef gsl_root_fdfsolver_type const* dynamic_type;
		static dynamic_type const newton; //< Newton method http://en.wikipedia.org/wiki/Newton%27s_method
		static dynamic_type const secant; //< Secant method http://en.wikipedia.org/wiki/Secant_method
		static dynamic_type const steffenson; static dynamic_type const steffensen; //< Seffensen method http://en.wikipedia.org/wiki/Steffensen%27s_method
		fdfsolver(function_fdf fdf, double initial_guess) : 
			fdf_(fdf), 
			pimpl_(gsl_root_fdfsolver_alloc(newton)){
			set(fdf_, initial_guess);
		}
		void set(function_fdf fdf, double guess){
			last_ = std::numeric_limits<double>::infinity();
			fdf_=fdf;
			gsl_root_fdfsolver_set(pimpl_, (gsl_function_fdf*)&fdf_, guess);
		}
		std::string name() const{return gsl_root_fdfsolver_name(pimpl_);}
		~fdfsolver(){gsl_root_fdfsolver_free(pimpl_);}
		fdfsolver& iterate(){
			last_ = root();
			gsl_root_fdfsolver_iterate(pimpl_);
			return *this;
		}
		double root() const{return gsl_root_fdfsolver_root(pimpl_);}
		protected:
		friend class boost::iterator_core_access;
		/// access current result with operator*(), warning, if called in the first step error is 'inf'
		result dereference() const{ //result operator*() const{
			return result(root(), fabs(root()-last_));
		}
		void increment() {this->iterate();}//solver& operator++(){return this->iterate();}
	};
	// definition of static member fields
	fdfsolver::dynamic_type const fdfsolver::newton     = gsl_root_fdfsolver_newton;
	fdfsolver::dynamic_type const fdfsolver::secant     = gsl_root_fdfsolver_secant;
	fdfsolver::dynamic_type const fdfsolver::steffenson = gsl_root_fdfsolver_steffenson; 
	fdfsolver::dynamic_type const fdfsolver::steffensen = gsl_root_fdfsolver_steffenson; 
#endif
}}
#endif
#ifdef _TEST_GSL_ROOT_HPP
#include<iostream>
#include "error.hpp"
#include<boost/lexical_cast.hpp>
#include<boost/lambda/lambda.hpp>
#include<boost/lambda/bind.hpp>
#include<boost/spirit/home/phoenix.hpp>
#include<boost/lexical_cast.hpp>
using std::cout; using std::endl;
double f(double x){
	return x*x-5.;
}
double g(double x){
	return 3.*x*x-5.;
}

double f2(double x, double a){
	return a*x*x-5.;
}

template<double(* FF)(double)>
std::string c_code(){
//	std::clog<<(void*)FF<<std::endl;
	return "f_"+boost::lexical_cast<std::string>((void*)FF);
}

double timestwo(double x){
	return x*2.;
}

template<>
std::string c_code<&f>(){return "x -> x*x-5";}

void myhandler (const char * reason, const char * file,
                        int line,
                        int gsl_errno){
	throw std::runtime_error(std::string(reason)+" at "+std::string(file)+":"+boost::lexical_cast<std::string>(line)+" error code "+boost::lexical_cast<std::string>(gsl_errno));
}


int main(){
	{
		using namespace boost::lambda;
		try{
			//throw std::runtime_error("aaa");
			gsl::root::solver sr(
				_1*_1-5., 
				gsl::root::interval(5., 10.)
			);
		}catch(...){
			std::clog<<"ja catched"<<std::endl;
		}
		for(
			gsl::root::solver sr(
				_1*_1-5., 
				gsl::root::interval(0., 5.)
			); 
			sr.x_upper() - sr.x_lower() > 0.0001; 
			sr.iterate()
		){
			std::cout<< sr.root()<<" "<<sqrt(5.)<<std::endl;
		}
		boost::function<double(double const&)> f_at_1(bind(&f2, bind(&timestwo, _1), 1.));
		for(gsl::root::solver sr(f_at_1, gsl::root::interval(0.,5.)); sr.x_upper() - sr.x_lower() > 0.0001; sr.iterate()){
			std::cout<< sr.root()<<" "<<sqrt(5.)<<std::endl;
		}
	}
	{
		using namespace boost::phoenix::arg_names;
		gsl::function f(arg1*arg1-5.);
		gsl::root::solver sr(f, gsl::root::interval(0.,5.));
//		++sr; ++sr; ++sr; ++sr; 
//		std::advance(sr, 10);
		gsl::result asr = *sr;
		std::cout << "solver " << *sr << " "<<f(*sr) << endl;

		gsl::function_fdf fdf(arg1*arg1 - 5., 2.*arg1);
		gsl::root::fdfsolver fdfsr(fdf, 5.2);
		++fdfsr; //++fdfsr; ++fdfsr; ++fdfsr;
//		std::advance(fdfsr, 5);
		//cout << *fdfsr << endl;
		//++fdfsr; 		++fdfsr; 		++fdfsr; 		++fdfsr; 		++fdfsr; 
		std::cout << "fdfsolver " << *fdfsr << " "<< fdf(*fdfsr) <<endl;
	}
	return 0;	
}
#endif

