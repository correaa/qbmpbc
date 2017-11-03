#ifdef COMPILATION_INSTRUCTIONS
ln -sf $0 $0.cpp && c++ -Wl,-rpath=$HOME/usr/lib$CLUSTER -Wall `#-Wfatal-errors` -L$HOME/usr/lib$CLUSTER -I$HOME/prj -I$HOME/usr/include $0.cpp  `pkg-config --libs gsl` -lboost_regex -lboost_filesystem -lboost_system -laspell -D_TEST_GSL_SF_HPP -o ./$0.cpp.x && ./$0.cpp.x $@
rm -f $0.cpp.x $0.cpp
exit
#endif
#ifndef GSL_SF_HPP
#define GSL_SF_HPP

#include<iostream> //for debug
//begin sf/result.hpp
#include "boost/just.hpp"
namespace gsl{
namespace sf{
	struct result : boost::just<double>{
		double error_;
		double const& value() const{return *this;}
		double const& error() const{return error_;}
		result(double value, double error = 0.) : boost::just<double>(value), error_(error){}
	};
}}

//begin fermi_dirac.hpp
#include<gsl/gsl_sf_fermi_dirac.h>
#include<boost/units/static_rational.hpp>
namespace gsl{
namespace sf{
	using boost::units::static_rational;
	inline double fermi_dirac_m1    (double const& x)              {return gsl_sf_fermi_dirac_m1(x);}
	inline double fermi_dirac_0     (double const& x)              {return gsl_sf_fermi_dirac_0 (x);}
	inline double fermi_dirac_1     (double const& x)              {return gsl_sf_fermi_dirac_1 (x);}
	inline double fermi_dirac_2     (double const& x)              {return gsl_sf_fermi_dirac_2 (x);}
	inline double fermi_dirac_int   (int const& j, double const& x){return gsl_sf_fermi_dirac_int(j, x);}
	//inline double fermi_dirac       (int const& j, double const& x){return fermi_dirac_int       (j, x);}
	inline double fermi_dirac_mhalf (double const& x)              {return gsl_sf_fermi_dirac_mhalf(x);}
	inline double fermi_dirac_half  (double const& x)              {return gsl_sf_fermi_dirac_half (x);}
	inline double fermi_dirac_3half(double const& x)               {return gsl_sf_fermi_dirac_3half(x);}
	//so far reproduces gsl
	template<int j> double fermi_dirac_int(double const& x){return fermi_dirac_int(j, x);}
	template<>      double fermi_dirac_int<-1>(double const& x){return fermi_dirac_m1(x);}
	template<>      double fermi_dirac_int< 0>(double const& x){return fermi_dirac_0 (x);}
	template<>      double fermi_dirac_int<+1>(double const& x){return fermi_dirac_1 (x);}
	template<>      double fermi_dirac_int<+2>(double const& x){return fermi_dirac_2 (x);}
	template<class Rational> double fermi_dirac_rational                    (double const& x);	
	template<>               double fermi_dirac_rational<static_rational<-1,2>::type>(double const& x){return fermi_dirac_mhalf(x);}
	template<>               double fermi_dirac_rational<static_rational<+1,2>::type>(double const& x){return fermi_dirac_half(x);}
	template<>               double fermi_dirac_rational<static_rational<+3,2>::type>(double const& x){return fermi_dirac_3half(x);}
	template<int i, int j>   double fermi_dirac_rational(double const& x){return fermi_dirac_rational<typename static_rational<i,j>::type>(x);}
	namespace fermi_dirac{
		template<int j> double F(double const& x){return fermi_dirac_int<j>;}
		template<int i, int j> double F(double const& x){return fermi_dirac_rational<typename static_rational<i,j>::type>(x);}
	}
}namespace special_functions = sf;}
//end fermi_dirac.hpp
//begin debye.hpp
#include<gsl/gsl_sf_debye.h>
#include<boost/math/special_functions/fpclassify.hpp>
namespace gsl{
namespace sf{
	inline double debye_1(double const& x){return gsl_sf_debye_1(x);}
	inline double debye_2(double const& x){return gsl_sf_debye_2(x);}
	inline double debye_3(double const& x){return ((boost::math::isinf)(x) and x>0)?0:gsl_sf_debye_3(x);}
	inline double debye_4(double const& x){return gsl_sf_debye_4(x);}
	inline double debye_5(double const& x){return gsl_sf_debye_5(x);}
	inline double debye_6(double const& x){return gsl_sf_debye_6(x);}
	namespace debye{
		template<unsigned j> double D   (double const& x);
		template<>           double D<1>(double const& x){return debye_1(x);}
		template<>           double D<2>(double const& x){return debye_2(x);}
		template<>           double D<3>(double const& x){return debye_3(x);}
		template<>           double D<4>(double const& x){return debye_4(x);}
		template<>           double D<5>(double const& x){return debye_5(x);}
		template<>           double D<6>(double const& x){return debye_6(x);}
	}
}}
#endif

/* warning: there is a erf and erfc function defined in math.h -lm and therefore in cmath std::erf[c] */
// more info: http://www.opengroup.org/onlinepubs/000095399/functions/erf.html
#include<gsl/gsl_sf_erf.h>
namespace gsl{
namespace sf{
	inline double erf_f(double x){return gsl_sf_erf(x);}
	inline result erf_e(double x){
		gsl_sf_result r;
		gsl_sf_erf_e(x, &r);
		return result(r.val, r.err);
	};
	inline result erf(double x){return erf_e(x);}

	inline double erfc_f(double x){return gsl_sf_erfc(x);}
	inline result erfc_e(double x){
		gsl_sf_result r;
		gsl_sf_erfc_e(x, &r);
		return result(r.val, r.err);
	};
	inline result erfc(double x){return erfc_e(x);}
}}

// more info at http://www.gnu.org/software/gsl/manual/html_node/Associated-Legendre-Polynomials-and-Spherical-Harmonics.html
#include<gsl/gsl_sf_legendre.h>
namespace gsl{
namespace sf{
	namespace spherical_harmonic{
		std::complex<double> Y(int const& l, int const& m, double const& theta, double const& phi){
			// calculated for positive m only, then modified for negative cases
			std::complex<double> ret =  
				gsl_sf_legendre_sphPlm(l, abs(m), cos(theta))
				* std::polar(1., m*phi)
			;
			// case m < 0
			if(m<0){
				if((abs(m) % 2) == 1) ret = -ret;
			}
			return ret;
		}
		template<int l, int m>  //simulate static binding with compile time checking
		std::complex<double> Y(double const& theta, double const& phi){
			BOOST_STATIC_ASSERT(l >= 0);
			BOOST_STATIC_ASSERT(l >= m && m >= -l);
			return Y(l, m, theta, phi);
		}
		/* needs boost proto to implement all this syntactic sugar correctly
		struct Y{
			int l_;
			int m_;
			Y(int l, int m) : l_(l), m_(m){}
			std::complex<double> operator()(double theta, double phi){return gsl::sf::spherical_harmonic::Y(l_, m_, theta, phi);}
		};
		*/
	}
}}

#ifdef _TEST_GSL_SF_HPP
#include<iostream>
#include "alf/latex.hpp"
#include "alf/gsl/error.hpp"
latex::ostream lout("sf.pdf");
using namespace std;

int main(){
	using gsl::sf::fermi_dirac::F;
	cout<<"F_{1/2}(0.5) = "<< F<1,2>(0.5) <<std::endl;
	using gsl::sf::debye::D;
	cout<<"D_3(0.5) = "<<D<3>(0.0001)<<std::endl;
	using namespace latex;
	lout << section("Error Functions");
	lout << subsection("Error Function");
	lout << "\\[ \\Erf{x} = \\frac{2}{\\pi} \\int_0^x e^{-t^2} dt \\]";
	{
		pgfplots::axis ax("xlabel = $x$, ylabel = $\\Erf{x}$");
		pgfplots::coordinates coords("no markers");
		for(double x = -10; x<10; x+=0.1){
			coords << std::pair<double, double>(x, gsl::sf::erf(x) );
		}
		pgfplots::axis ax_err("xlabel = $x$, ylabel = {error}");
		{
			pgfplots::coordinates coords("no markers");
			for(double x = -10; x<10; x+=0.01){
				coords << std::make_pair(x, gsl::sf::erf(x).error() );
			}
			ax_err << coords;
		}
		{
			tikzpicture tkp;
			ax << coords; tkp << ax;
			lout << tkp;
		}
		{
			tikzpicture tkp;
			tkp << ax_err;
			lout << tkp;
		}
	}
	lout << subsection("Complementary Error Function");
	lout << "\\[ \\Erfc{x} = 1 - \\Erf{x}  \\]";
	{
		pgfplots::axis ax("xlabel = $x$, ylabel = $\\Erfc{x}$");
		pgfplots::coordinates coords("no markers");
		for(double x = -10; x<10; x+=0.1){
			coords << std::pair<double, double>(x, gsl::sf::erfc(x) );
		}
		pgfplots::axis ax_err("xlabel = $x$, ylabel = {error}");
		{
			pgfplots::coordinates coords("no markers");
			for(double x = -10; x<10; x+=0.01){
				coords << std::make_pair(x, gsl::sf::erfc(x).error() );
			}
			ax_err << coords;
		}
		{
			tikzpicture tkp;
			ax << coords; tkp << ax;
			lout << tkp;
		}
		{
			tikzpicture tkp;
			tkp << ax_err;
			lout << tkp;
		}
	}
	lout << subsection("Spherical Harmonics Function");
	{
		pgfplots::axis ax("xlabel = $\\phi$, ylabel = {$Y_{\\ell=2 m}(\\theta = 0, \\phi)$}");
		pgfplots::coordinates coords("no markers");
		for(double phi = -3.14/2.; phi<3.14/2.; phi+=0.1){
			std::clog << cos(phi) << std::endl;
			coords << std::pair<double, double>(phi, real(gsl::sf::spherical_harmonic::Y<2,0>(0., phi) ));
		}
		tikz::picture tkzp;
		ax << coords;
		lout << (tkzp << ax) << newline;
//		lout << math(std::complex<double>(1.,2.)) << newline;
		lout << math(gsl::sf::spherical_harmonic::Y(4, 1,0.1234, 0.5678)) << newline;
		cout << gsl::sf::spherical_harmonic::Y(4, -3, 0.1234, 0.5678) << endl;
		cout << gsl::sf::spherical_harmonic::Y<4, -3>(0.1234, 0.5678) << endl;
	}
	return 0;
}
#endif
// astyle --brackets=attach --indent=tab --indent-col1-comments --pad-oper --delete-empty-lines --add-brackets --keep-one-line-statements --convert-tabs --align-pointer=type --unpad-paren  qbox.hpp
// Editor modelines  -  http://www.wireshark.org/tools/modelines.html
// Local variables:
// c-basic-offset: 4
// tab-width: 4
// indent-tabs-mode: t
// truncate-lines: 1
// End:
/* vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent: */

