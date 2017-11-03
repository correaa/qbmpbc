#ifdef COMPILATION_INSTRUCTIONS
	rm -f .$0.cpp.x && ln -sf $0 .$0.cpp && c++ -Wfatal-errors -Wall .$0.cpp -L$HOME/lib `pkg-config --libs gsl` -D_TEST_GSL_FUNCTION_FDF_HPP -o .$0.x && ./.$0.x $1 $2 $3 $4 $5 $6 $7 $8 $9
	rm -f $0.cpp.x $0.cpp
	exit
#endif
#ifndef GSL_FUNCTION_FDF_HPP_
#define GSL_FUNCTION_FDF_HPP_

#include "../gsl/derivative.hpp"
namespace gsl{
class function_fdf : public boost::function<double(double)>, public gsl_function_fdf{
	boost::function<double(double)> df_;
	public: 
	template<class F, class G> 
	function_fdf(F f, G df) : 
		boost::function<double(double)>(f), df_(df){
		gsl_function_fdf::f      = &function_fdf::free_function_;
		gsl_function_fdf::df     = &function_fdf::free_derivative_;
		gsl_function_fdf::fdf    = &function_fdf::free_function_derivative_;
		gsl_function_fdf::params = this;
	}
	function_fdf(function const& f) : 
		boost::function<double(double)>((boost::function<double(double)> const&)f){
		gsl_function_fdf::f      = &function::free_function_;
		gsl_function_fdf::df     = &function_fdf::free_derivative_;
		gsl_function_fdf::fdf    = &function_fdf::free_function_derivative_;
		gsl_function_fdf::params     = this;
		assert(not (this->empty()));
	}
	double cancelation_error(double const& x) const{ return gsl::derivative::central(*this, x) - df_(x);}
	double df(double const& x) const{return df_(x);}
	double f(double const& x) const{return boost::function<double(double)>::operator()(x);}
	static double free_function_(double x, void* self /*void *params*/){
		assert(not ((gsl::function_fdf*)self)->empty());
		return ((function_fdf*)self)->operator()(x);
	}
	static double free_derivative_(double x, void* self /*void *params*/){
		assert(not ((gsl::function_fdf*)self)->empty());
		return ((function_fdf*)self)->df(x);
	}
	static void free_function_derivative_(double x, void* self /**params*/, double *y, double *dy){
		//todo: be smarter here and evaluate both at the same function (if available).
		*y = free_function_(x, self);
		*dy = free_derivative_(x, self);
	}
};
}
#endif

#ifdef _TEST_GSL_FUNCTION_FDF_HPP
#include<iostream> 
using std::cout; using std::endl;
double freef(double d){return d*d-5.+4.*d*d*d;}
double freedf(double d){return 2.*d+12.*d*d;}

int main(){
	gsl::function_fdf fdf(freef, freedf);
	cout << fdf(2.2) <<endl;
	cout << fdf.df(2.2) <<endl;
	cout << fdf.cancelation_error(10.)/fdf.df(10.) <<endl;
	return 0;
}
#endif
