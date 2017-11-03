#ifndef _ONE_PARTICLE_WAVEFUNCTION_
#define _ONE_PARTICLE_WAVEFUNCTION_
#include "field.hpp"
#include "space.hpp"
#include "alf/mpl_utility.hpp"

typedef boost::rational<int> rational;
using boost::mpl::rational_;
typedef std::complex<double> complex;

template<class RationalFlux>
class wf_mpbc{
	wf_mpbc();
	contravariariant kappa_point_;
	grid real_space_grid;
	field<complex> real_space() const;
	multi_array<complex, 3> ma;
};

template<>
wf_mpbc::wf_mpbc(double energy_cut_off){
	
};

template<>
field<complex> wf_mpbc<rational_<0>::type >::real_space() const{
	const_multi_array_ref<complex, 3> wfXYZ(wf_.sd(ispin, ikp)->c().cvalptr(state*wf_.sd(ispin, ikp)->c().mloc()),
					                                    extents
						                                    [wf_.sd(ispin, ikp)->basis().np(0)]
						                                    [wf_.sd(ispin, ikp)->basis().np(1)]
						                                    [wf_.sd(ispin, ikp)->basis().np(2)]
	);
}

#endif
#ifdef _TEST_b

space const real_space;

vector<&real_space> 

#endif

