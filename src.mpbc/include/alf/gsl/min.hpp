#ifdef COMPILATION_INSTRUCTIONS
ln -sf $0 .$0.cpp && c++ `#-Wfatal-errors` -Wall .$0.cpp -I$HOME/prj -I$HOME/prj/alf -L$HOME/lib `pkg-config --libs gsl` -lboost_system -D_TEST_GSL_MIN_HPP -o .$0.x && ./.$0.x $@
rm -f $0.x $0.cpp
exit
#endif
#ifndef GSL_MIN_HPP
#define GSL_MIN_HPP
namespace gsl{
namespace minimization{
	typedef const gsl_min_fminimizer_type* type;
	class minimizer{
	};
}
}
#endif

#ifdef _TEST_GSL_MIN_HPP
#include<boost/assign.hpp>
#include "../gsl/interpolation.hpp"
#include"alf/latex.hpp"
using namespace latex;
ostream lout("minimization.pdf");

int main(){
	lout 
		<< *section("Numerical Minimization")
		<< "Suppose we have a function (for example given by a spline interpolation), " << par;
	std::map<double, double> m = boost::assign::map_list_of
		( 1., 4.)
		( 3., 3.)
		( 4., 2.)
		( 6., 8.)
		( 7., 9.)
		(10., 9.)
	;
	gsl::interpolation::spline p(m);
	pgfplots::axis ax("title = {Spline Interpolation}");
	lout << (ax << p);

	gsl::minimization::minimizer
	return 0;
}
#endif

