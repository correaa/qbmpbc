////////////////////////////////////////////////////////////////////////////////
//
// FourierTransform3.h
//
////////////////////////////////////////////////////////////////////////////////

#ifndef FOURIERTRANSFORM3_H
#define FOURIERTRANSFORM3_H

#include <complex>
#include <vector>
#include <boost/array.hpp>
//#include "fftw.h" //include before fftw3.hpp
#include "alf/fftw3.hpp"
#include "Timer.h"

//class Basis;
//class Context;

using std::complex;
using boost::multi_array;

using boost::multi_array_ref;
using boost::const_multi_array_ref;
using boost::extents;

class FourierTransform{
	private:
		//const Context& ctxt_;
		//const Basis& basis_;
		int nprocs_, myproc_;
		boost::array<unsigned, 3u> np_;
	public:
		//		const Context& context(void) const { return ctxt_; }

		static void backward_pbc(
				const_multi_array_ref<complex<double>, 3> const& wfXYZ, 
				multi_array<complex<double>, 3>& wfxyz){
			fftw3::plan(wfXYZ, wfxyz, FFTW_BACKWARD).execute_normalize();
		}
		static void forward_pbc(
				const_multi_array_ref<complex<double>, 3> const& wfxyz, 
				multi_array<complex<double>, 3>& wfXYZ){
			fftw3::plan(wfxyz, wfXYZ, FFTW_BACKWARD).execute_normalize();
		}
		static void z_backward_mpbc(
				const_multi_array_ref<complex<double>, 2> const& wfXZh, 
				multi_array<complex<double>, 2>& wfXzh){
#ifndef _SHORTFFT
			fftw3::plan(-1, wfXZh, wfXzh, FFTW_BACKWARD).execute_normalize();
#else
			multi_array<complex<double>, 2> wfXZh_short(extents[wfXZh.shape()[0]][wfXZh.shape()[1]/2]);
			multi_array<complex<double>, 2> wfXzh_short(extents[wfXzh.shape()[0]][wfXzh.shape()[1]/2]);

			/* copy wfXZh to wfXZh_short */
			for(unsigned i=0; i < wfXZh_short.shape()[0]; i++){
				for(unsigned j=0; j < wfXZh_short.shape()[1]; j++){
					if (j < wfXZh_short.shape()[1]/2) 
						wfXZh_short[i][j] = wfXZh[i][j];
					else
						wfXZh_short[i][j] = wfXZh[i][j+wfXZh_short.shape()[1]];
				}
			}

			/* anti-aliasing in rcpr space */
			for(unsigned i=0; i < wfXZh_short.shape()[0]; i++){
				for(unsigned j=0; j < wfXZh_short.shape()[1]; j++){
					if ( ( j >= floor(wfXZh_short.shape()[1]/4.0) ) && 
							( j <= ceil(wfXZh_short.shape()[1]*3/4.0) ) ) 
						wfXZh_short[i][j] = 0;
				}
			}

			fftw3::plan(-1, wfXZh_short, wfXzh_short, FFTW_BACKWARD).execute_normalize();

			/* copy wfXzh_short to wfXzh */
			for(unsigned i=0; i < wfXzh.shape()[0]; i++){
				for(unsigned j=0; j < wfXzh.shape()[1]; j++){
					if (j < wfXzh_short.shape()[1]/2) 
						wfXzh[i][j] = wfXzh_short[i][j];
					else if (j >= wfXzh.shape()[1]-wfXzh_short.shape()[1]/2)
						wfXzh[i][j] = wfXzh_short[i][j-wfXzh_short.shape()[1]];
					else
						wfXzh[i][j] = 0;
				}
			}
#endif
		}
		static void z_forward_mpbc(
				const_multi_array_ref<complex<double>,2u> const& dwfXzh, 
				multi_array_ref<complex<double>, 2u>& dwfXZh){
#ifndef _SHORTFFT
			fftw3::plan(-1, dwfXzh, dwfXZh, FFTW_FORWARD).execute_normalize();
#else
			multi_array<complex<double>, 2> dwfXzh_short(extents[dwfXzh.shape()[0]][dwfXzh.shape()[1]/2]);
			multi_array<complex<double>, 2> dwfXZh_short(extents[dwfXZh.shape()[0]][dwfXZh.shape()[1]/2]);

			/* copy dwfXzh to dwfXzh_short */
			for(unsigned i=0; i < dwfXzh_short.shape()[0]; i++){
				for(unsigned j=0; j < dwfXzh_short.shape()[1]; j++){
					if (j < dwfXzh_short.shape()[1]/2) 
						dwfXzh_short[i][j] = dwfXzh[i][j];
					else
						dwfXzh_short[i][j] = dwfXzh[i][j+dwfXzh_short.shape()[1]];
				}
			}

			fftw3::plan(-1, dwfXzh_short, dwfXZh_short, FFTW_FORWARD).execute_normalize();

			/* anti-aliasing in rcpr space */
			for(unsigned i=0; i < dwfXZh_short.shape()[0]; i++){
				for(unsigned j=0; j < dwfXZh_short.shape()[1]; j++){
					if ( ( j >= floor(dwfXZh_short.shape()[1]/4.0) ) && ( j <= ceil(dwfXZh_short.shape()[1]*3/4.0) ) ) 
						dwfXZh_short[i][j] = 0;
				}
			}

			/* copy dwfXZh_short to dwfXZh */
			for(unsigned i=0; i < dwfXZh.shape()[0]; i++){
				for(unsigned j=0; j < dwfXZh.shape()[1]; j++){
					if (j < dwfXZh_short.shape()[1]/2) 
						dwfXZh[i][j] = dwfXZh_short[i][j];
					else if (j >= dwfXZh.shape()[1]-dwfXZh_short.shape()[1]/2)
						dwfXZh[i][j] = dwfXZh_short[i][j-dwfXZh_short.shape()[1]];
					else
						dwfXZh[i][j] = 0;
				}
			}
#endif
		}
		static void xy_backward_mpbc(
				const_multi_array_ref<complex<double>, 3> const& wfXYz, 
				multi_array_ref<complex<double>, 3>& wfxyz){
			fftw3::plan(+2, wfXYz, wfxyz, FFTW_BACKWARD).execute_normalize();
		}
		static void xy_forward_mpbc(
				const_multi_array_ref<complex<double>, 3> const& wfxyz, 
				multi_array_ref<complex<double>, 3>& dwfXYz){
			fftw3::plan(+2, wfxyz, dwfXYz, FFTW_FORWARD).execute_normalize();
		}
		int np0() const { return np_[0]; }
		int np1() const { return np_[1]; }
		int np2() const { return np_[2]; }
		unsigned np012() const { return np_[0] * np_[1] * np_[2]; }

		unsigned np012loc() const{return np012();}
		unsigned np012loc(int /*iproc*/) const{return np012();}
		int np2_first(int&) const{assert(0); return np2();}
		int np2_loc(int&) const{assert(0); return np2();}
		//int np012loc(int iproc) const{assert(0); return 0;}


		//void reset_timers(void);
		//Timer tm_f_map, tm_f_fft, tm_f_pack, tm_f_mpi, tm_f_zero, tm_f_unpack,
		//      tm_b_map, tm_b_fft, tm_b_pack, tm_b_mpi, tm_b_zero, tm_b_unpack;

		~FourierTransform(){}
		FourierTransform (
				//const Basis &basis,
				boost::array<unsigned, 3u> const& np) //int np0, int np1, int np2)      
			:   //ctxt_(basis.context()), //basis_(basis),
				np_(np){
					//assert(ctxt_.npcol() == 1);
					//nprocs_ = ctxt_.size();
					//myproc_ = ctxt_.myproc();
				}
};


#endif
// ~/usr/bin-hera/astyle --brackets=attach --indent=tab --indent-col1-comments --pad-oper --delete-empty-lines --add-brackets --keep-one-line-statements --convert-tabs --align-pointer=type --unpad-paren  qbox.hpp
// Editor modelines  -  http://www.wireshark.org/tools/modelines.html
// Local variables:
// c-basic-offset: 4
// tab-width: 4
// indent-tabs-mode: t
// truncate-lines: 1
// End:
/* vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent: */
/* vi cindent commands: == (current line), or gg=G (entire file) */
