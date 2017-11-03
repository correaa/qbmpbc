#ifdef COMPILATION_INSTRUCTIONS
ln -sf $0 .$0.cpp && 
c++ -std=c++0x -Wfatal-errors -Wall -Wextra -Wno-unused-variable -Wno-unused-parameter -Wno-mismatched-tags -Wno-ignored-qualifiers -Wno-sign-compare -Wno-char-subscripts -Wno-unused-function \
.$0.cpp -I$HOME/prj -I$HOME/prj/alf -I$HOME/usr/include -L$HOME/usr/lib `pkg-config --libs gsl` -lboost_system -lboost_regex -laspell -D_TEST_GSL_MULTI_MIN_HPP -o .$0.x && ./.$0.x $@
rm -f $0.x $0.cpp
exit
#endif
#ifndef GSL_MULTI_MIN_HPP
#define GSL_MULTI_MIN_HPP
#include<boost/array.hpp>
#include<boost/function.hpp>
#include<gsl/gsl_multimin.h>
#include<iostream>
#include<boost/iterator/iterator_facade.hpp>
#include<cmath>
using std::clog; using std::endl;
namespace gsl{
namespace multi_min{
	template<size_t NArgs>
	class function : 
		public boost::function<double(boost::array<double, NArgs> const&)>, 
		public gsl_multimin_function{
		public:
		typedef boost::array<double, NArgs> array;
		template<class F>
		function(F function_object) :
			boost::function<double(array const&)>(function_object){
			gsl_multimin_function::f      = &function::free_function_;
			gsl_multimin_function::n      = NArgs; 
			gsl_multimin_function::params = this;
			//clog << "ctor" << endl;
		}
		function(function const& f) : 
			boost::function<double(array const&)>(
				(boost::function<double(array const&)> const&)f
			){
			//clog << "copy" << endl;
			gsl_multimin_function::f      = &function::free_function_;
			gsl_multimin_function::n      = NArgs;
			gsl_multimin_function::params = this;
			assert(not (this->empty()));
		}
		~function(){
			//clog << "destroyed" << endl;
		}
		private:
		static double free_function_(const gsl_vector* x, void* self){
			assert(not ((function*)self)->empty());
			array xarr; std::copy(x->data, x->data+NArgs, xarr.begin());
			return ((function*)self)->operator()(xarr);
		}
	};
	template<size_t NArgs>
	class minimizer : 
		public boost::iterator_facade<
			minimizer<NArgs>, 
			std::pair<boost::array<double, NArgs>, double > const, 
			boost::forward_traversal_tag, 
			std::pair<boost::array<double, NArgs>, double > const /*noref*/ 
		>
	{
		gsl_multimin_fminimizer* const pimpl_;
		function<NArgs> const f_;
		public:
		typedef boost::array<double, NArgs> array;
		minimizer(function<NArgs> const& f, array const& starting_point, array const& step_size) :
			pimpl_(gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, NArgs)),
			f_(f)
		{
			gsl_vector const x_ = gsl_vector_const_view_array(starting_point.data(), NArgs).vector; //check
			gsl_vector const step_size_ = gsl_vector_const_view_array(step_size.data(), NArgs).vector; //check
			gsl_multimin_fminimizer_set(
				pimpl_, 
				&(const_cast<function<NArgs>&>(f_)), 
				&x_, //copied
				&step_size_ //copied
			);
		}
		~minimizer(){
			gsl_multimin_fminimizer_free(pimpl_);
		}

		// thin gsl interface http://www.gnu.org/software/gsl/manual/html_node/Multimin-Iteration.html
		void iterate(){gsl_multimin_fminimizer_iterate(pimpl_);}
		array const& x() const{
			//array const& ret = *(array const*)gsl_multimin_fminimizer_x(pimpl_)->data;
			/*
			std::copy(
				gsl_multimin_fminimizer_x(pimpl_)->data, 
				gsl_multimin_fminimizer_x(pimpl_)->data + NArgs,
				ret.begin()
			);*/
			return *(array const*)gsl_multimin_fminimizer_x(pimpl_)->data;
			//return (array const&)(gsl_multimin_fminimizer_x(pimpl_)->data);
		}
		double minimum() const{return gsl_multimin_fminimizer_minimum(pimpl_);}
		// void restart() ... there is no restart for fminimizer (there is restart only for fdfminimizer
		array const& gradient() const{return *(array const*)(gsl_multimin_fminimizer_x(pimpl_)->data);}
		double characteristic_step_size() /*size()*/ const{return gsl_multimin_fminimizer_size(pimpl_);} // Minimizer specific characteristic size for the minimizer s, name could be confusing characteristic_size()? step_size()? 

		// iterator interface
		protected:
		std::pair<array, double> dereference() const{
				return std::pair<array, double>(x(), minimum());
		}
		//bool equal(minimizer const& other);?
		void increment(){this->iterate();}
    friend class boost::iterator_core_access;
	};
}
}
#endif
#ifdef _TEST_GSL_MULTI_MIN_HPP
#include<iostream>
#include<boost/spirit/home/phoenix.hpp>
using std::cout ; using std::endl;
struct a_class{
	double p0, p1, p2, p3, p4;
	a_class(double p0, double p1, double p2, double p3, double p4) : 
		p0(p0), p1(p1), p2(p2), p3(p3), p4(p4){}
	double operator()(boost::array<double, 2> const& d) const{
		double x = d[0];
		double y = d[1];
		return 
			p2 * (x - p0) * (x - p0) +
			p3 * (y - p1) * (y - p1) + p4
		;
	}
};
struct f_impl{
	template<class Args> struct result{typedef double type;};
	double operator()(double const& x) const{
		return std::pow(x,4)-std::pow(x,2)*3.+x + cos(x*5.)*5.;
	}
};

template<class Functor>
struct convex{
	boost::phoenix::function<Functor> f_;
	convex(Functor const& f) : f_(f){}
	double operator()(double const& x) const{
		using namespace boost::phoenix::arg_names;
		double ret = f_(arg1)(x);
		for(double dx = 0.0001; dx<10.; dx+=1.){ 
			for(double dy = 0.0001; dy<10.; dy+=1.0){ 
				gsl::multi_min::minimizer<2> mnr(
					f_(x - arg1[0]*arg1[0]) 
					+ ((f_(x + arg1[1]*arg1[1]) - f_(x - arg1[0]*arg1[0]))
						/(arg1[1]*arg1[1] + arg1[0]*arg1[0])) * (arg1[0]*arg1[0]) 
					,
					(boost::array<double, 2>){{dx, dy}}, 
					(boost::array<double, 2>){{0.1, 0.1}}
				);
				for(unsigned i = 0; i!=40; ++i){
					++mnr;
				}
				if(mnr->second<ret) ret = mnr->second;
			}
		}
		return ret;
	}
};

double f(boost::array<double, 2> const& d){
		double 
			p0 = 1.0,
			p1 = 2.0, 
			p2 = 10.0, 
			p3 = 20.0, 
			p4 = 30.0
		;
		double x = d[0];
		double y = d[1];
		return 
			p2 * (x - p0) * (x - p0) +
			p3 * (y - p1) * (y - p1) + p4
		;
}

#include "alf/gsl/derivative.hpp"
#include "alf/latex.hpp"
using namespace latex;
latex::ostream lout("multi_min.pdf");
int main(){
	/*
	lout << "Given a function of two variables" << newline;
	a_class a(1.0,2.0, 10.0, 20.0, 30.0);
	pgf::plots::axis ax("title = {$10.  (x - 1.)^2 + 20. (y - 2.)^2 + 30.$}, xlabel = {$x$}, ylabel = {$y$}");
	pgf::plots::coordinates3 coords3("no markers, surf, shader=interp");
	for(double x=4.; x > 0.; x-=0.5){
		for(double y=0.; y < 4.; y+=0.5){
			coords3 << pgfplots::triple(x, y, a((boost::array<double, 2>){{x,y}}));
		}
		coords3 << '\n'; //endr
	}
	ax << coords3;
	gsl::multi_min::minimizer<2> mnr(
		a, 
		(boost::array<double, 2>){{2.0, 3.0}}, 
		(boost::array<double, 2>){{0.01, 0.01}}
	);
	pgf::plots::coordinates3 coords_min;
	unsigned n = 0;
	for(unsigned i = 0; i!=50; ++i){
		mnr.iterate(); ++n;
		cout << mnr->first[0] << ", " << mnr->first[1] <<" -> " << mnr->second << endl;
		coords_min << pgfplots::triple(mnr->first[0], mnr->first[1], mnr->second);
	}
	ax << coords_min;
	tikz::picture tkzp; tkzp << ax;
	lout << tkzp;
	lout << par << "The minimum is found at $x \\to" << mnr->first[0] <<"$, $y \\to " << mnr->first[1] <<"$ with a value of $" << mnr->second<<"$ after "<<(int)n<<" iterations. "; 
	*/

	lout << 
		"Another application of the multivariable minimizer is to convexify a function. "
		"Suppose we have a function. " << 
		equation("f(x) = x^4 - x^2/2 + x/3") << 
		"The following function is the convexified version of it" <<
		equation("F(x) = \\inf_{x_1 < x < x_2}\\left[ f(x_1) + \\frac{f(x_2) - f(x_1)}{x_2 - x_1} (x - x_1)\\right]")
	;
	{
		f_impl f;
		convex<f_impl> f_cvx(f);
		pgf::plots::axis ax("xlabel = {$x$}, ylabel = {$y$}");
		pgf::plots::coordinates coords("no markers");
		pgf::plots::coordinates coords_cvx("no markers");
		for(double x=-3; x<3; x+=0.01){
			coords << pgfplots::pair(x, f(x));
			coords_cvx << pgfplots::pair(x, f_cvx(x));
		}
		ax << coords << coords_cvx;
		tikz::picture tkzp; tkzp << ax;
		lout << tkzp;
		{	
			auto 
				df = gsl::derivative::central(f),
				df_cvx = gsl::derivative::central(f_cvx)
			;
			pgf::plots::axis ax("xlabel = {$x$}, ylabel = {$y$}");
			pgf::plots::coordinates coords("no markers");
			pgf::plots::coordinates coords_cvx("no markers");
			for(double x=-3; x<3; x+=0.01){
				coords << pgfplots::pair(x, df(x));
				coords_cvx << pgfplots::pair(x, df_cvx(x));
			}
			ax << coords << coords_cvx;
			tikz::picture tkzp; tkzp << ax;
			lout << tkzp;
		}
	}
	return 0;
}
#endif

