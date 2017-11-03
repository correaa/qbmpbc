#ifdef compile_instructions
ln -sf $0 $0.cpp && c++ $0.cpp -Wall -Wextra `#-Wfatal-errors` -Wno-unused-variable -Wno-unused-parameter -Wno-ignored-qualifiers \
`pkg-config --libs gsl` -I$HOME/usr/include \
-I$HOME/prj \
-I$HOME/soft/qbox-1.52.3/src \
-Wl,-rpath=$HOME/usr/lib -L$HOME/usr/lib  -lboost_regex -lboost_system -lboost_filesystem -lhunspell \
-L$HOME/soft/qbox-1.52.3/src -lqb_serial -lxerces-c -lfftw \
-D_TEST_QBOX_SPECIES_HPP  -o $0.x && ./$0.x $@ 
rm -rf $0.cpp $0.x
exit
#endif
#ifndef QBOX_SPECIES_HPP
#define QBOX_SPECIES_HPP
#include "Species.h"
#include "SpeciesReader.h"
#include<boost/filesystem.hpp>
//#include<boost/mpi.hpp>
namespace qbox{
using namespace boost::filesystem;
class species : virtual public Species{
	/*
	std::vector< gsl::interpolation::spline > v;
	*/
	public:
	species() : Species(null_context_, "noname") {
		/*
		for(unsigned ell=0; ell<=lmax; ++l){
			assert(rps_.size() == vps_spl_[ell].size());
			std::map<double, double> m;
			for(unsigned xidx=0; xidx < rps_.size(); ++xidx){
				m[rps_[xidx]] = vps_spl_[ell][xidx];
			}
			vps.push_back(gsl::interpolation::spline(m));
		}*/
	}
	double vpsr(int l, double r) const{double ret; const_cast<species*>(this)->Species::vpsr(l, r, ret); return ret;}
	double dvpsr(int l, double r) const{double ignore, ret; const_cast<species*>(this)->Species::dvpsr(l, r, ignore, ret); return ret;}
	double vlocg(double q) const{double ret; const_cast<species*>(this)->Species::vlocg(q, ret); return ret;}
	double dvlocg(double q) const{double ignore, ret; const_cast<species*>(this)->Species::dvlocg(q, ignore, ret); return ret;}
	double vnlg(int l, double q) const{double ret; const_cast<species*>(this)->Species::vnlg(l, q, ret); return ret;}
	double dvnlg(int l, double q) const{double ignore, ret; const_cast<species*>(this)->Species::dvnlg(l,q, ignore, ret); return ret;}
	double phi(int l, double r) const{double ret; const_cast<species*>(this)->Species::phi(l, r, ret); return ret;}
	unsigned nlm() const{return const_cast<species*>(this)->Species::nlm();}
	bool non_local() const{return const_cast<species*>(this)->Species::non_local();}
	double eself() const{return const_cast<species*>(this)->Species::eself();}
	double const& rmax() const{return rps_.back();}
	unsigned rps_size() const{return rps_.size();}
	double const& rps_back() const{return rps_.back();}
	double const& rps_front() const{return rps_.front();}
	path uri() const{return Species::uri();}
	static Context null_context_;
	public:
	static species load(path p) try{
		species ret;
		SpeciesReader reader(null_context_);
		reader.readSpecies(ret, p.string());
		ret.initialize(1.5 /*hardcoded in qbox*/); 
		return ret;
	}catch(std::runtime_error& e){
		throw;
	}
	int const& ndft() const{return ndft_;}
};
Context species::null_context_;
}
#endif

#ifndef QBOX_SPECIES_IO_HPP
#define QBOX_SPECIES_IO_HPP
#include "alf/latex.hpp"
#include "boost/string_cast.hpp" //to_string
#include<boost/units/cmath.hpp>

namespace qbox{
	latex::ostream& operator<<(latex::ostream& lout, species const& s){
		using namespace latex;
		lout << par
			<< " number of non-local projectors $\\sum_{\\ell \\neq \\ell_\\text{local}} 2\\ell+1$ [\\verb+nlm+] : " << s.nlm() << newline
			<< " is non-local [\\verb+non_local+] : " << s.non_local() << newline
			<< " size of radial FFT array (power of 2 $>$ \\verb+rdftmin/deltar_+ ) [\\verb+ndft+] : " << s.ndft() << newline
			<< " name [\\verb+name+] : " << s.name() << newline
			<< " Uniform Resource Identifier [\\verb+uri+] : " << s.uri() << newline
			<< " symbol [\\verb+symbol+] : " << s.symbol() << newline
			<< " atomic number [\\verb+atomic_number+] : " << s.atomic_number() << newline
			<< " mass (in u.a.m.u. carbon = 12.0) [\\verb+mass+] : " << s.mass() << newline
			<< " description [\\verb+description+] : " << s.description() << newline
			<< " valence charge ($Z$) [\\verb+zval+] : " << s.zval() << newline
			<< " largest angular momentum ($\\ell \\leq \\ell_\\text{max}$) [\\verb+lmax+] : " << s.lmax() << newline
			<< " angular momentum taken as local (in general $\\ell_\\text{local} = \\ell_\\text{max}$) [\\verb+llocal+] : " << s.llocal() << newline
			<< " number of semi-local quadrature points [\\verb+nquad+] : " << s.nquad() << newline
			<< " end of semi-local quadrature interval [\\verb+rquad+] : " << s.rquad() << newline
			<< " mesh spacing for potentials and wavefunctions (in $a_0$) [\\verb+deltar+] : " << s.deltar() << newline
			<< " cutoff radius of gaussian pseudocharge (fixed externally(?), hard coded(?)) [\\verb+rcps+] : " << s.rcps() << newline
			<< " self Coulomb energy (?) ($\\frac{zval^2}{\\sqrt{2\\pi}*rcps}$) [\\verb+eself+] : " << s.eself() << newline
			<< " mesh size [\\verb+rps_.size()+] : " << s.rps_size() << newline
			<< " maximum r in mesh [\\verb+rps_.back()+] : " << s.rps_back() << newline
			<< " minimum r in mesh [\\verb+rps_.front()+] : " << s.rps_front() << newline 
		;
		lout << noindent;
		{
			typedef latex::pgfplots::axis::units<atomic::length, atomic::energy> axis; axis ax("no markers, title = {Pseudopotential $V_\\ell(r)$}, xlabel = {$r$}, legend pos = south east");
			typedef axis::coordinates_type coordinates; 
			typedef coordinates::pair_type pair;
			for(int ell=0; ell<= s.lmax(); ++ell){
				coordinates coords;
				for(double x=0; x <= std::min(s.rmax(), 10.); x+=0.01){
					coords << pair(x*atomic::bohr, s.vpsr(ell, x)*atomic::hartree);
				}
				ax << std::make_pair("$\\ell = "+ boost::to_string(ell) +"$ " + ((ell==s.llocal())?"(local)":""), coords);
			}
			coordinates coords;
			for(double x=0.2; x <= std::min(s.rmax(), 10.); x+=0.01){
				coords << pair(x*atomic::bohr, -s.zval()/x*atomic::hartree);
			}
			ax << std::make_pair("$ - Z / r$", coords);
			tikz::picture p; p<<ax;
			lout << p;
		}
		{
			typedef latex::pgfplots::axis::units<atomic::length, atomic::force> axis; axis ax("no markers, title = {Pseudopotential derivative $V'_\\ell(r)$}");
			for(int ell=0; ell<= s.lmax(); ++ell){
				typedef axis::coordinates_type coordinates; coordinates coords;
				typedef coordinates::pair_type pair;
				for(double x=0; x <= std::min(s.rmax(), 10.); x+=0.01){
					coords << pair(x*atomic::bohr, s.dvpsr(ell, x)*atomic::hartree/atomic::bohr);
				}
				ax << std::make_pair("$\\ell = "+ boost::to_string(ell) +"$ " + ((ell==s.llocal())?"(local)":""), coords);
			}
			tikz::picture p; p<<ax;
			lout << p;
		}
		lout << newline;
		{
			typedef latex::pgfplots::axis::units<
				atomic::wavenumber, 
				atomic::energy
			> axis; axis ax("no markers, title = {Pseudopotential $V_\\ell(q)$}, legend pos=south east");
			typedef axis::coordinates_type coordinates;
			for(int ell=0; ell<= s.lmax(); ++ell){
				coordinates coords;
				typedef coordinates::pair_type pair;
				for(double x=0; x <= 20.; x+=0.01){
					coords << pair(x/atomic::bohr, s.vnlg(ell, x)*atomic::hartree);
				}
				ax << std::make_pair("$\\ell = "+ boost::to_string(ell) +"$ " + ((ell==s.llocal())?"(local)":""), coords);
			}
			{
				coordinates coords;
				typedef coordinates::pair_type pair;
				for(double q=0; q <= 20.; q+=0.01){
					coords << pair(q/atomic::bohr, s.vlocg(q)*atomic::hartree);
				}
				ax << std::make_pair("$V_\\text{local}(q)$", coords); 
			}
			tikz::picture p; p<<ax;
			lout << p;
		}
		{
			typedef latex::pgfplots::axis::units<
				atomic::wavenumber, 
				multiply_typeof_helper<atomic::energy, atomic::length>::type
			> axis; axis ax("no markers, title = {Pseudopotential $V'_\\ell(q)$}, legend pos=north east");
			typedef axis::coordinates_type coordinates;
			for(int ell=0; ell<= s.lmax(); ++ell){
				coordinates coords;
				typedef coordinates::pair_type pair;
				for(double x=0; x <= 20.; x+=0.01){
					coords << pair(x/atomic::bohr, s.dvnlg(ell, x)*atomic::hartree*atomic::bohr);
				}
				ax << std::make_pair("$\\ell = "+ boost::to_string(ell) +"$ " + ((ell==s.llocal())?"(local)":""), coords);
			}
			{
				coordinates coords;
				typedef coordinates::pair_type pair;
				for(double q=0; q <= 20.; q+=0.01){
					coords << pair(q/atomic::bohr, s.dvlocg(q)*atomic::hartree*atomic::bohr);
				}
				ax << std::make_pair("$V_\\text{local}(q)$", coords); 
			}
			tikz::picture p; p<<ax;
			lout << p;
		}
		lout << newline;
		{
			typedef latex::pgfplots::axis::units<
				atomic::length, 
				atomic::wavefunction_amplitude
				//atomic::number_density
			> axis; axis ax("no markers, title = {Pseudowavefunction $\\phi_\\ell(r)$, $\\int_0^\\infty |\\phi(r)|^2 r^2 dr = 1$}");
			for(int ell=0; ell<= s.lmax(); ++ell){
				std::map<quantity<atomic::length>, quantity<atomic::number_density> > m;
				typedef axis::coordinates_type coordinates; coordinates coords;
				typedef coordinates::pair_type pair;
				for(double x=0; x <= std::min(s.rmax(), 10.); x+=0.01){
					coords << pair(x*atomic::bohr, s.phi(ell, x)/(1.*pow<static_rational<3,2> >(atomic::bohr)));
					m[x*atomic::bohr] = s.phi(ell, x)*s.phi(ell, x)/(1.*pow<static_rational<3,1> >(atomic::bohr));
				}
				ax << std::make_pair("$\\ell = "+ boost::to_string(ell) +"$ " + ((ell==s.llocal())?"(local)":""), coords);
//				gsl::interpolation::units::spline<atomic::length, atomic::number_density> s(m);
//				boost::function<quantity<atomic::number_density>(quantity<atomic::length>)> F(s);
//				boost::numeric::interval<quantity<atomic::length> > ival(0.*atomic::bohr, 20.*atomic::bohr);
//				std::clog << F(1.*atomic::bohr) << std::endl;
//				boost::phoenix::function<gsl::interpolation::units::spline<atomic::length, atomic::number_density> > sph(s);
//				using boost::phoenix::arg_names::arg1;
//				std::clog << "ell = " << ell << ", ||phi|| = " << gsl::integration::qag(sph(arg1)*arg1*arg1, ival) << std::endl;
			}
			tikz::picture p; p<<ax;
			lout << p;
		}
		return lout;
	}
}
#endif

#ifdef _TEST_QBOX_SPECIES_HPP


#include "alf/gsl/units/interpolation.hpp"
#include "alf/gsl/units/integration.hpp"
#include "alf/gsl/error.hpp"
latex::ostream lout("test/carbon_pbe.pdf");
using namespace latex;
int main(int argc, char* argv[]){
	gsl::error::unset_handler();
	lout << section("Carbon PBE");
	qbox::species s(qbox::species::load("test/carbon_pbe.xml"));
	lout << s;
	lout << par
		<< " number of non-local projectors $\\sum_{\\ell \\neq \\ell_\\text{local}} 2\\ell+1$ [\\verb+nlm+] : " << s.nlm() << newline
		<< " is non-local [\\verb+non_local+] : " << s.non_local() << newline
		<< " size of radial FFT array (power of 2 $>$ \\verb+rdftmin/deltar_+ ) [\\verb+ndft+] : " << s.ndft() << newline
		<< " name [\\verb+name+] : " << s.name() << newline
		<< " Uniform Resource Identifier [\\verb+uri+] : " << s.uri() << newline
		<< " symbol [\\verb+symbol+] : " << s.symbol() << newline
		<< " atomic number [\\verb+atomic_number+] : " << s.atomic_number() << newline
		<< " mass (in u.a.m.u. carbon = 12.0) [\\verb+mass+] : " << s.mass() << newline
		<< " description [\\verb+description+] : " << s.description() << newline
		<< " valence charge ($Z$) [\\verb+zval+] : " << s.zval() << newline
		<< " largest angular momentum ($\\ell \\leq \\ell_\\text{max}$) [\\verb+lmax+] : " << s.lmax() << newline
		<< " angular momentum taken as local (in general $\\ell_\\text{local} = \\ell_\\text{max}$) [\\verb+llocal+] : " << s.llocal() << newline
		<< " number of semi-local quadrature points [\\verb+nquad+] : " << s.nquad() << newline
		<< " end of semi-local quadrature interval [\\verb+rquad+] : " << s.rquad() << newline
		<< " mesh spacing for potentials and wavefunctions (in $a_0$) [\\verb+deltar+] : " << s.deltar() << newline
		<< " cutoff radius of gaussian pseudocharge (fixed externally(?), hard coded(?)) [\\verb+rcps+] : " << s.rcps() << newline
		<< " self Coulomb energy (?) ($\\frac{zval^2}{\\sqrt{2\\pi}*rcps}$) [\\verb+eself+] : " << s.eself() << newline
		<< " mesh size [\\verb+rps_.size()+] : " << s.rps_size() << newline
		<< " maximum r in mesh [\\verb+rps_.back()+] : " << s.rps_back() << newline
		<< " minimum r in mesh [\\verb+rps_.front()+] : " << s.rps_front() << newline 
	;
	lout << noindent;
	{
		typedef latex::pgfplots::axis::units<atomic::length, atomic::energy> axis; axis ax("no markers, title = {Pseudopotential $V_\\ell(r)$}, xlabel = {$r$}, legend pos = south east");
		typedef axis::coordinates_type coordinates; 
		typedef coordinates::pair_type pair;
		for(int ell=0; ell<= s.lmax(); ++ell){
			coordinates coords;
			for(double x=0; x <= std::min(s.rmax(), 10.); x+=0.01){
				coords << pair(x*atomic::bohr, s.vpsr(ell, x)*atomic::hartree);
			}
			ax << std::make_pair("$\\ell = "+ boost::to_string(ell) +"$ " + ((ell==s.llocal())?"(local)":""), coords);
		}
		coordinates coords;
		for(double x=0.2; x <= std::min(s.rmax(), 10.); x+=0.01){
			coords << pair(x*atomic::bohr, -s.zval()/x*atomic::hartree);
		}
		ax << std::make_pair("$ - Z / r$", coords);
		tikz::picture p; p<<ax;
		lout << p;
	}
	{
		typedef latex::pgfplots::axis::units<atomic::length, atomic::force> axis; axis ax("no markers, title = {Pseudopotential derivative $V'_\\ell(r)$}");
		for(int ell=0; ell<= s.lmax(); ++ell){
			typedef axis::coordinates_type coordinates; coordinates coords;
			typedef coordinates::pair_type pair;
			for(double x=0; x <= std::min(s.rmax(), 10.); x+=0.01){
				coords << pair(x*atomic::bohr, s.dvpsr(ell, x)*atomic::hartree/atomic::bohr);
			}
			ax << std::make_pair("$\\ell = "+ boost::to_string(ell) +"$ " + ((ell==s.llocal())?"(local)":""), coords);
		}
		tikz::picture p; p<<ax;
		lout << p;
	}
	{
		typedef latex::pgfplots::axis::units<
			atomic::wavenumber, 
			atomic::energy
		> axis; axis ax("no markers, title = {Pseudopotential $V_\\ell(q)$}, legend pos=south east");
		typedef axis::coordinates_type coordinates;
		for(int ell=0; ell<= s.lmax(); ++ell){
			coordinates coords;
			typedef coordinates::pair_type pair;
			for(double x=0; x <= 20.; x+=0.01){
				coords << pair(x/atomic::bohr, s.vnlg(ell, x)*atomic::hartree);
			}
			ax << std::make_pair("$\\ell = "+ boost::to_string(ell) +"$ " + ((ell==s.llocal())?"(local)":""), coords);
		}
		{
			coordinates coords;
			typedef coordinates::pair_type pair;
			for(double q=0; q <= 20.; q+=0.01){
				coords << pair(q/atomic::bohr, s.vlocg(q)*atomic::hartree);
			}
			ax << std::make_pair("$V_\\text{local}(q)$", coords); 
		}
		tikz::picture p; p<<ax;
		lout << p;
	}
	{
		typedef latex::pgfplots::axis::units<
			atomic::wavenumber, 
			multiply_typeof_helper<atomic::energy, atomic::length>::type
		> axis; axis ax("no markers, title = {Pseudopotential $V'_\\ell(q)$}, legend pos=north east");
		typedef axis::coordinates_type coordinates;
		for(int ell=0; ell<= s.lmax(); ++ell){
			coordinates coords;
			typedef coordinates::pair_type pair;
			for(double x=0; x <= 20.; x+=0.01){
				coords << pair(x/atomic::bohr, s.dvnlg(ell, x)*atomic::hartree*atomic::bohr);
			}
			ax << std::make_pair("$\\ell = "+ boost::to_string(ell) +"$ " + ((ell==s.llocal())?"(local)":""), coords);
		}
		{
			coordinates coords;
			typedef coordinates::pair_type pair;
			for(double q=0; q <= 20.; q+=0.01){
				coords << pair(q/atomic::bohr, s.dvlocg(q)*atomic::hartree*atomic::bohr);
			}
			ax << std::make_pair("$V_\\text{local}(q)$", coords); 
		}
		tikz::picture p; p<<ax;
		lout << p;
	}
	{
		typedef latex::pgfplots::axis::units<
			atomic::length, 
			atomic::wavefunction_amplitude
			//atomic::number_density
		> axis; axis ax("no markers, title = {Pseudowavefunction $\\phi_\\ell(r)$, $\\int_0^\\infty |\\phi(r)|^2 r^2 dr = 1$}");
		for(int ell=0; ell<= s.lmax(); ++ell){
			std::map<quantity<atomic::length>, quantity<atomic::number_density> > m;
			typedef axis::coordinates_type coordinates; coordinates coords;
			typedef coordinates::pair_type pair;
			for(double x=0; x <= std::min(s.rmax(), 10.); x+=0.01){
				coords << pair(x*atomic::bohr, s.phi(ell, x)/(1.*pow<static_rational<3,2> >(atomic::bohr)));
				m[x*atomic::bohr] = s.phi(ell, x)*s.phi(ell, x)/(1.*pow<static_rational<3,1> >(atomic::bohr));
			}
			ax << std::make_pair("$\\ell = "+ boost::to_string(ell) +"$ " + ((ell==s.llocal())?"(local)":""), coords);
			gsl::interpolation::units::spline<atomic::length, atomic::number_density> s(m);
			boost::function<quantity<atomic::number_density>(quantity<atomic::length>)> F(s);
			boost::numeric::interval<quantity<atomic::length> > ival(0.*atomic::bohr, 20.*atomic::bohr);
			std::clog << F(1.*atomic::bohr) << std::endl;
			boost::phoenix::function<gsl::interpolation::units::spline<atomic::length, atomic::number_density> > sph(s);
			using boost::phoenix::arg_names::arg1;
			std::clog << "ell = " << ell << ", ||phi|| = " << gsl::integration::qag(sph(arg1)*arg1*arg1, ival) << std::endl;
		}
		tikz::picture p; p<<ax;
		lout << p;
	}
	return 0;
}
#endif
/* vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent: */

