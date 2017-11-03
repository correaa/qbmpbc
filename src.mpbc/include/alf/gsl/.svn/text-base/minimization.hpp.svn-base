#ifdef COMPILATION_INSTRUCTIONS
ln -sf $0 .$0.cpp && c++ -std=c++0x `#-Wfatal-errors` -Wall -Wno-unused-variable .$0.cpp -I$HOME/prj -I$HOME/prj/alf -L$HOME/usr/lib `pkg-config --libs gsl` -lboost_system -lboost_regex -D_TEST_GSL_MIN_HPP -o .$0.x && ./.$0.x $@
rm -f $0.x $0.cpp
exit
#endif
#ifndef GSL_MIN_HPP
#define GSL_MIN_HPP
#include "../gsl/error.hpp"
#include<gsl/gsl_min.h>
#include "boost/just.hpp"
#include<boost/numeric/interval.hpp>
#include"../gsl/function.hpp"
#include<boost/iterator/iterator_facade.hpp>
#include<boost/utility.hpp> //noncopyable
#include<iostream>
using std::clog; using std::endl;
namespace gsl{
namespace minimization{
	typedef const gsl_min_fminimizer_type* type;
	type const golden_section = gsl_min_fminimizer_goldensection;
	typedef boost::numeric::interval<double> interval;
	struct bracket : interval, boost::just<double>{
		bracket(interval const& iv, double guess) : 
			interval(iv), 
			boost::just<double>(guess)
		{}
	};
	class minimizer :
		public boost::iterator_facade<
			minimizer, 
			std::pair<double, double> const, 
			boost::forward_traversal_tag, 
			std::pair<double, double> const /*noref*/ 
		>,
		boost::noncopyable
		{
		gsl_min_fminimizer* const pimpl_;
		function f_;
		public:
		minimizer(function const& f, bracket const& b) : 
			pimpl_(gsl_min_fminimizer_alloc(golden_section)), 
			f_(f){
			//clog << "here" << endl;
			gsl::error::status s = gsl_min_fminimizer_set(pimpl_, &f_, (double const&)b, lower(b), upper(b));
			if(s==gsl::error::failure) std::clog << "warning " << "[" + boost::lexical_cast<std::string>(lower(b)) + "<" + boost::lexical_cast<std::string>(b) + "<" + boost::lexical_cast<std::string>(upper(b)) + "] does not bracket a minimum, minimum find is not guarrantied" << std::endl
				//throw gsl::error::exception(
				//	s, 
				//	"[" + lexical_cast<string>(lower(b)) + "<" + lexical_cast<string>(b) + "<" + lexical_cast<string>(upper(b)) + "] " 
				//	"does not bracket a minimum"
				//)
			;
			//clog << "there" << endl;
		}
		//int gsl_min_fminimizer_set_with_values (gsl_min_fminimizer * s, gsl_function * f, double x_minimum, double f_minimum, double x_lower, double f_lower, double x_upper, double f_upper)
		std::string name() const{return gsl_min_fminimizer_name(pimpl_);}
		gsl::error::status iterate(){
			return gsl::error::code(gsl_min_fminimizer_iterate(pimpl_));
		}
		double x_minimum() const{return gsl_min_fminimizer_x_minimum(pimpl_);}
		double x_lower() const{return gsl_min_fminimizer_x_lower(pimpl_);}
		double x_upper() const{return gsl_min_fminimizer_x_upper(pimpl_);}
		double f_minimum() const{return gsl_min_fminimizer_f_minimum(pimpl_);}
		double f_upper() const{return gsl_min_fminimizer_f_upper(pimpl_);}
		double f_lower() const{return gsl_min_fminimizer_f_lower(pimpl_);}
		bracket x() const{return 
			bracket(
				interval(
					x_lower(),
					x_upper()
				),
				x_minimum()
			);
		}
		std::pair<double, double> dereference() const{
			return std::pair<double, double>(
				x_minimum(), f_minimum()
			);
		}
		void increment() {this->iterate();}
		~minimizer(){gsl_min_fminimizer_free(pimpl_);}
	};
}}
#endif

#ifdef _TEST_GSL_MIN_HPP
#include<boost/assign.hpp>
#include "../gsl/interpolation.hpp"
#include"alf/latex.hpp"
#include<iostream>
using std::cout; using std::endl;

template<class PhoenixFunctor>
struct legendre{
	boost::phoenix::function<PhoenixFunctor> f_;
	//boost::function<double(double)> f_;
	legendre(PhoenixFunctor const& f) : f_(f){}
	template<class Args> struct result{typedef double type;};
	double operator()(double p) const{
		gsl::minimization::bracket bkt(gsl::minimization::interval(-20.,20.),5.);
		using namespace boost::phoenix::arg_names;
		gsl::minimization::minimizer mnr(f_(_1)-p*_1, bkt);
		for(unsigned i=0; i!=100; ++i){
			++mnr;
		}
		return mnr->first;
	}
};

using namespace latex;
ostream lout("minimization.pdf");
int main(){
	gsl::error::unset_handler();
	lout 
		<< *section("Numerical Minimization")
		<< "Suppose we have a function (for example given by a spline interpolation), " << par;
	std::map<double, double> m = boost::assign::map_list_of
		( 1., 8.)
		( 3., 3.)
		( 4., 2.)
		( 6., 8.)
		( 7., 9.)
		(10., 9.)
	;
	gsl::interpolation::spline p(m);
	pgf::plots::axis ax("title = {Spline Interpolation}");
	pgf::plots::coordinates coords("green");
	{
		gsl::minimization::bracket bkt(gsl::minimization::interval(1.,10.),5.5);
		gsl::minimization::minimizer mnr(p, bkt);
		for(unsigned i=0; i!=100; ++i){
			cout 
				<< mnr->first << " " << mnr->second << " "
				<< mnr.x_lower() << " " << mnr.x_upper() << " "
				<< width(mnr.x()) << endl;
			++mnr;
			coords << pgf::plots::pair(mnr->first, mnr->second);
		}
	}
	ax << coords;
	tikz::picture tkzp; tkzp << (ax << p);
	lout << tkzp;
	{
		lout << 
			section("Legendre Transform") <<
			equation("f^\\star(p) = \\inf_x(f(x)- px)");
		std::map<double, double> m = boost::assign::map_list_of
			( 1., 15.)
			( 3., 3.)
			( 4., 2.)
			( 6., 8.)
			( 7., -5.)
			(10., -6.)
		;
		gsl::interpolation::spline p(m);
		{
			pgf::plots::axis ax("title = {Spline Interpolation}");
			tikz::picture tkzp; tkzp << (ax << p);
			lout << tkzp;
		}
		{
			legendre<gsl::interpolation::spline> L(p);
			pgf::plots::axis axL("title = {Legendre Transform}");
			pgf::plots::coordinates coords("no markers");
			for(double pp=0; pp<10; pp+=0.01){
				coords << pgf::plots::pair(pp, L(pp));
			}
			tikz::picture tkzpL; tkzpL << (axL << coords);
			lout << tkzpL;

			legendre<legendre<gsl::interpolation::spline> > LL(L);
			pgf::plots::axis axLL("title = {LLegendre Transform}");
			pgf::plots::coordinates coordsLL("no markers");
			for(double pp=0; pp<10; pp+=0.01){
				coordsLL << pgf::plots::pair(pp, LL(pp));
			}
			tikz::picture tkzpLL; tkzpLL << (axLL << coordsLL);
			lout << tkzpLL;
		}
	}
	return 0;
}
#endif

// astyle --brackets=attach --indent=tab --indent-col1-comments --delete-empty-lines --add-brackets --keep-one-line-statements --convert-tabs --align-pointer=type --unpad-paren fit_linear.hpp
// Editor modelines  -  http://www.wireshark.org/tools/modelines.html
// Local variables:
// c-basic-offset: 4
// tab-width: 4
// indent-tabs-mode: t
// truncate-lines: 1
// End:
/* vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent: */

