////////////////////////////////////////////////////////////////////////////////
//
// Copyright (c) 2008 The Regents of the University of California
//
// This file is part of Qbox
//
// Qbox is distributed under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 2 of
// the License, or (at your option) any later version.
// See the file COPYING in the root directory of this distribution
// or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// SlaterDet.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SlaterDet.C,v 1.53 2008/09/08 15:56:19 fgygi Exp $
#include "SlaterDet.h"
#include "FourierTransform.h"
#include "Context.h"
#include "blas.h" // daxpy
#include "Base64Transcoder.h"
#include "SharedFilePtr.h"
#include "Timer.h"

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <boost/static_assert.hpp>
#include <complex>

#ifdef _USE_HDF5
#include "alf/hdf5.hpp"
#endif

#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>
#include<fstream>

#ifdef _REGULARGDISTRIB
#include "alf/fftw3.hpp"
#endif

using namespace std;

////////////////////////////////////////////////////////////////////////////////
SlaterDet::SlaterDet(const Context& ctxt, D3vector kpoint) : ctxt_(ctxt),
 c_(ctxt)
{
  //cout << ctxt.mype() << ": SlaterDet::SlaterDet: ctxt.mycol="
  //     << ctxt.mycol() << " basis_->context(): "
  //     << basis_->context();
  my_col_ctxt_ = 0;
  for ( int icol = 0; icol < ctxt_.npcol(); icol++ )
  {
    Context* col_ctxt = new Context(ctxt_,ctxt_.nprow(),1,0,icol);
    ctxt_.barrier();
    if ( icol == ctxt_.mycol() )
      my_col_ctxt_ = col_ctxt;
    else
      delete col_ctxt;
  }
  //cout << ctxt_.mype() << ": SlaterDet::SlaterDet: my_col_ctxt: "
  //     << *my_col_ctxt_;
  basis_ = new Basis(*my_col_ctxt_,kpoint);
}

////////////////////////////////////////////////////////////////////////////////
SlaterDet::SlaterDet(const SlaterDet& rhs) : ctxt_(rhs.context()),
	my_col_ctxt_(new Context(*(rhs.my_col_ctxt_))), 
	basis_(new Basis(*(rhs.basis_))),
	c_(rhs.c_)
{}

////////////////////////////////////////////////////////////////////////////////
SlaterDet::~SlaterDet()
{
  delete basis_;
  delete my_col_ctxt_;
  // cout << ctxt_.mype() << ": SlaterDet::dtor: ctxt=" << ctxt_;
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::resize(const UnitCell& cell, const UnitCell& refcell,
  double ecut, unsigned nst)
{
  // Test in next line should be replaced by test on basis min/max indices
  // to signal change in basis vectors
  //if ( basis_->refcell().volume() != 0.0 && !refcell.encloses(cell) )
  //{
    //cout << " SlaterDet::resize: cell=" << cell;
    //cout << " SlaterDet::resize: refcell=" << basis_->refcell();
    //throw SlaterDetException("could not resize: cell not in refcell");
  //}

  try
  {
    // create a temporary copy of the basis and of the coefficient matrix
    Basis btmp(*basis_);
    ComplexMatrix ctmp(c_);

    // perform normal resize operations, possibly resetting contents of c_
#ifdef _NOANTIALIASING
    basis_->resize(cell,refcell,ecut);
#else
    basis_->resize(cell,refcell,4.0*ecut);
#endif
    basis_->resize(cell,refcell,ecut);
    occ_.resize(nst);
    eig_.resize(nst);

    const int mb = basis_->maxlocalsize();
    const unsigned m = ctxt_.nprow() * mb;
    const int nb = nst/ctxt_.npcol() + (nst%ctxt_.npcol() > 0 ? 1 : 0);

    // Determine if plane wave coefficients must be reset after the resize
    // This is needed if the dimensions of the matrix c_ must be changed
    const bool needs_reset =
		m!=c_.m()   || 
		nst!=c_.n() || 
		mb!=c_.mb() ||
		nb!=c_.nb();

    c_.resize(m,nst,mb,nb);

    if ( needs_reset )
      reset();

    // check if data can be copied from temporary copy
    // It is assumed that nst and ecut are not changing at the same time
    // Only the cases where one change at a time occurs is covered

    // consider only cases where the dimensions are finite
    if ( c_.m() > 0 && c_.n() > 0 )
    {

      // first case: only nst has changed
      if ( c_.m() == ctmp.m() && c_.n() != ctmp.n() )
      {
        //cout << "SlaterDet::resize: c_m/n=   "
        //     << c_.m() << "/" << c_.n() << endl;
        //cout << "SlaterDet::resize: ctmp_m/n=" << ctmp.m()
        //     << "/" << ctmp.n() << endl;
        // nst has changed, basis is unchanged
        // copy coefficients up to min(n_old, n_new)
        if ( c_.n() < ctmp.n() )
        {
          c_.getsub(ctmp,ctmp.m(),c_.n(),0,0);
        }
        else
        {
          c_.getsub(ctmp,ctmp.m(),ctmp.n(),0,0);
        }
        gram();
      }
      // second case: basis was resized, nst unchanged
      if ( btmp.ecut() > 0.0 && basis_->ecut() > 0.0 &&
           c_.m() != ctmp.m() && c_.n() == ctmp.n() )
      {
        // transform all states to real space and interpolate
        int np0 = max(basis_->np(0),btmp.np(0));
        int np1 = max(basis_->np(1),btmp.np(1));
        int np2 = max(basis_->np(2),btmp.np(2));
        cout << " SlaterDet::resize: grid: np0/1/2: "
		     << np0 << " " << np1 << " " << np2 << endl;
        // FourierTransform tf1(oldbasis, new basis grid)
        // FourierTransform tf2(newbasis, new basis grid)

#ifndef _REGULARGDISTRIB
        FourierTransform ft1(btmp,np0,np1,np2);
        FourierTransform ft2(*basis_,np0,np1,np2);
        // allocate real-space grid
        valarray<complex<double> > tmpr(ft1.np012loc());
        // transform each state from old basis to grid to new basis
        for ( int n = 0; n < nstloc(); n++ )
        {
          ft1.backward(ctmp.cvalptr(n*ctmp.mloc()),&tmpr[0]);
          ft2.forward(&tmpr[0], c_.valptr(n*c_.mloc()));
        }
#else
        // need to replace with FFTW3
        clog << "SlaterDet::resize() interpolate wf by FFTW3 not implemented yet!" << endl;
        assert(0);
#endif
      }
    }
  }
  catch ( bad_alloc )
  {
    cout << " bad_alloc exception caught in SlaterDet::resize" << endl;
    throw;
  }
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::reset(void)
{
  // initialize coefficients with lowest plane waves
  if ( c_.n() <= basis_->size() )
  {
    // initialize c_
    c_.clear();
    const double s2i = 1.0 / sqrt(2.0);

    // for each n, find the smallest g vector and initialize
    int ismallest = 0;
    // on each process, basis.isort(ismallest) is the index of the smallest
    // local g vector
    for ( unsigned n = 0; n < c_.n(); n++ )
    {
      double value = 1.0;
      if ( basis().real() && n != 0 )
        value = s2i;

      // find process row holding the smallest g vector
#ifndef _MPBC
      double kpg2 = basis_->kpg2_PBC(basis_->isort(ismallest));
#else
      double kpg2 = basis_->kpg2_MPBC(basis_->isort(ismallest));
#endif
      // cout << "smallest vector on proc " << ctxt_.mype()
      //      << " has norm " << g2 << endl;
      int minrow, mincol;
      ctxt_.dmin('c',' ',1,1,&kpg2,1,&minrow,&mincol,1,-1,-1);

      // find column hosting state n
      int pc = c_.pc(n);
      int pr = minrow;
      if ( pr == ctxt_.myrow() )
      {
        int iii = basis_->isort(ismallest);
        ismallest++; // increment on entire process row
        if ( pc == ctxt_.mycol() )
        {
          // cout << " n=" << n << " on process "
          //      << pr << "," << pc
          //      << " vector " << basis_->idx(3*iii) << " "
          //      << basis_->idx(3*iii+1) << " "
          //      << basis_->idx(3*iii+2) << " norm="
          //      << basis_->g2(iii) << " "
          //      << value << endl;
          int jjj = c_.m(n) * c_.nb() + c_.y(n);
          int index = iii+c_.mloc()*jjj;
          c_[index] = complex<double> (value,0.0);
        }
      }
    }
  }
}

void SlaterDet::init(void)
{
  // initialize coefficients with lowest plane waves
  if ( c_.n() <= basis_->size() )
  {
    // initialize c_
    c_.clear();
    const double s2i = 1.0 / sqrt(2.0);

    // for each n, find the smallest g vector and initialize
    int ismallest = 0;
    // on each process, basis.isort(ismallest) is the index of the smallest
    // local g vector
    for ( int n = 0; n < c_.n(); n++ )
    {
      double value = 1.0;
      if ( basis().real() && n != 0 )
        value = s2i;

      // find process row holding the smallest g vector
      double kpg2 = basis_->kpg2(basis_->isort(ismallest));
      // cout << "smallest vector on proc " << ctxt_.mype()
      //      << " has norm " << g2 << endl;
      int minrow, mincol;
      ctxt_.dmin('c',' ',1,1,&kpg2,1,&minrow,&mincol,1,-1,-1);

      // find column hosting state n
      int pc = c_.pc(n);
      int pr = minrow;
      if ( pr == ctxt_.myrow() )
      {
        int iii = basis_->isort(ismallest);
        ismallest++; // increment on entire process row
        if ( pc == ctxt_.mycol() )
        {
          // cout << " n=" << n << " on process "
          //      << pr << "," << pc
          //      << " vector " << basis_->idx(3*iii) << " "
          //      << basis_->idx(3*iii+1) << " "
          //      << basis_->idx(3*iii+2) << " norm="
          //      << basis_->g2(iii) << " "
          //      << value << endl;
          int jjj = c_.m(n) * c_.nb() + c_.y(n);
          int index = iii+c_.mloc()*jjj;
          c_[index] = complex<double> (value,0.0);
        }
      }
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
void SlaterDet::compute_density(FourierTransform& 
#ifndef _CHANGE_FT_CLASS
ft
#endif
, double weight, double* rho) const
{
	
	assert(basis_->cell().volume() > 0.0);
  const double omega_inv = 1.0 / basis_->cell().volume();

#ifndef _CHANGE_FT_CLASS
  //Timer tm_ft, tm_rhosum;
  // compute density of the states residing on my column of ctxt_
  assert(occ_.size() == c_.n());
  vector<complex<double> > tmp(ft.np012loc());
	const int np012loc = ft.np012loc();


  if ( basis_->real() )
  {
#ifdef _NOANTIALIASING
    assert(0); //do not allow real_ basis if _NOANTIALIASING
#endif
    // transform two states at a time
    for ( int n = 0; n < nstloc()-1; n++, n++ )
    {
      // global n index
      const int nn = ctxt_.mycol() * c_.nb() + n;
      const double fac1 = weight * omega_inv * occ_[nn];
      const double fac2 = weight * omega_inv * occ_[nn+1];

      if ( fac1 + fac2 > 0.0 )
      {
        //tm_ft.start();
        ft.backward(c_.cvalptr(n*c_.mloc()),
                    c_.cvalptr((n+1)*c_.mloc()),&tmp[0]);
        //tm_ft.stop();
        const double* psi = (double*) &tmp[0];
        int ii = 0;
        //tm_rhosum.start();
        for ( int i = 0; i < np012loc; i++ )
        {
          const double psi1 = psi[ii];
          const double psi2 = psi[ii+1];
          rho[i] += fac1 * psi1 * psi1 + fac2 * psi2 * psi2;
          ii++; ii++;
        }
        //tm_rhosum.start();
      }
    }
    if ( nstloc() % 2 != 0 )
    {
      const int n = nstloc()-1;
      // global n index
      const int nn = ctxt_.mycol() * c_.nb() + n;
      const double fac1 = weight * omega_inv * occ_[nn];

      if ( fac1 > 0.0 )
      {
        ft.backward(c_.cvalptr(n*c_.mloc()),&tmp[0]);
        const double* psi = (double*) &tmp[0];
        int ii = 0;
        for ( int i = 0; i < np012loc; i++ )
        {
          const double psi1 = psi[ii];
          rho[i] += fac1 * psi1 * psi1;
          ii++; ii++;
        }
      }
    }
  }
  else
  {
#endif //_CHANGE_FT_CLASS

#ifndef _REGULARGDISTRIB
    // original Qbox code
    // only one transform at a time
    for ( int n = 0; n < nstloc(); n++ )
    {
      // global n index
      const int nn = ctxt_.mycol() * c_.nb() + n;
      const double fac = weight * omega_inv * occ_[nn];

      if ( fac > 0.0 )
      {
        ft.backward(c_.cvalptr(n*c_.mloc()),&tmp[0]);
        for ( int i = 0; i < np012loc; i++ )
          rho[i] += fac * norm(tmp[i]);
      }
    }
#else
    // regular distribution of G vectors

#ifndef _NOANTIALIASING
    BOOST_STATIC_ASSERT(0); //antialiasing must be disabled to use regular distribution of G vectors
#endif

    // use FFTW3 transforms, by-passing qbox rods design
    typedef std::complex<double> complex;
    using boost::const_multi_array_ref;
    using boost::multi_array_ref;
		using boost::multi_array;
    using boost::extents;
    using fftw3::plan;

    int np0=basis_->np(0);
    int np1=basis_->np(1);
    int np2=basis_->np(2);
    //int np012 = np0*np1*np2;

    for ( int n = 0; n < nstloc(); n++ )
    {
      // global n index
      const int nn = ctxt_.mycol() * c_.nb() + n;
      const double fac = weight * omega_inv * occ_[nn];

#ifndef _MPBC /* normal PBC */
      vector<complex> tmp(np012);
      const_multi_array_ref<complex, 3> wfXYZ (c_.cvalptr(n*c_.mloc()), extents[np0][np1][np2]);
      multi_array_ref<complex, 3>       wfxyz (&tmp[0],                 extents[np0][np1][np2]);

      plan wf_backward(wfXYZ, wfxyz, FFTW_BACKWARD, FFTW_ESTIMATE);

      //clog << "n = " << n << "fac = " << fac << "occ = " << occ_[nn] << endl;
      if ( fac > 0.0 )
      {
        wf_backward();
        for ( int i = 0; i < np012; i++ )
          rho[i] += fac * norm(tmp[i]);
      }
#else /* MPBC */
      //vector<complex> tmp_i(np012);
      //multi_array_ref<complex, 2>     wfXzh (&tmp_i[0],               extents[np0][np1* np2]);
      //multi_array_ref<complex, 3>     wfXYz (&tmp_i[0],               extents[np0][np1][np2]);

      const_multi_array_ref<complex, 2> wfXZh (c_.cvalptr(n*c_.mloc()), extents[np0][np1* np2]);
      multi_array<complex, 2>           wfXzh (                         extents[np0][np1* np2]);
      multi_array_ref<complex, 3>       wfXYz (wfXzh.origin(),          extents[np0][np1][np2]);
      multi_array<complex, 3>           wfxyz (extents[np0][np1][np2]);

      /* change these functions to functions defined in FourierTransform3.h */

      /* backward 1D FFT along zhat direction: wfXZh -> wfXzh = wfXYz (tmp) */
      /* 2D array, transform last one dimension (-1) */
      //plan wf_z_backward(-1, wfXZh, wfXzh, FFTW_BACKWARD, FFTW_ESTIMATE); 

      /* backward 2D FFT along x and y directions: wfXYz (tmp) -> wfxyz (tmp_r) */
      /* 3D array, transform first two dimensions (+2) */
      //plan wf_xy_backward(+2, wfXYz, wfxyz, FFTW_BACKWARD, FFTW_ESTIMATE); 

      //clog << "n = " << n << "fac = " << fac << "occ = " << occ_[nn] << endl;
      if ( fac > 0.0 )
      {
        //wf_z_backward();
        //wf_xy_backward();
        /* normalization factor to match that of PBC */
        //for ( int i = 0; i < wfxyz.num_elements(); i++ ) wfxyz.origin()[i] /= sqrt(1.0*np1);
        //clog<<"wfxyz[0][0][0] = "<<wfxyz[0][0][0]<<endl;
        //clog<<"wfxyz[0][0][1] = "<<wfxyz[0][0][1]<<endl;
        //clog<<"wfxyz[0][1][0] = "<<wfxyz[0][1][0]<<endl;

        FourierTransform::z_backward_mpbc(wfXZh, wfXzh);
        FourierTransform::xy_backward_mpbc(wfXYz, wfxyz);
        /* normalization factor to match that of PBC */
        for ( unsigned i = 0; i < wfxyz.num_elements(); i++ ) wfxyz.origin()[i] *= sqrt(np1*np2*np0);

        for ( unsigned i = 0; i < wfxyz.num_elements(); i++ ) rho[i] += fac * norm(wfxyz.origin()[i]);
        //clog << "adding density from state" << n << endl;
      }
#endif
    }
#endif

#ifndef _CHANGE_FT_CLASS
  }
#endif //_CHANGE_FT_CLASS

  if(boost::filesystem::exists("print_density"))
  {
     int np0=basis_->np(0);
     int np1=basis_->np(1);
     int np2=basis_->np(2);
     int np012 = np0*np1*np2;
     typedef std::complex<double> complex;
     boost::multi_array<complex, 3> rhoxyz (boost::extents[np0][np1][np2]);
     clog << "np012 = " << np012; // << "np012loc = " << np012loc << endl;
     for ( int i = 0; i < np012; i++ ) *(rhoxyz.origin()+i) = complex(rho[i],0.0);
#ifdef _USE_HDF5
//     Write_WF_to_HDF5(basis_,rhoxyz,"rho_real",0);
#else
//     Write_WF_to_Matlab(rhoxyz,"rho_real",0);
#endif
  }
  // cout << "SlaterDet: compute_density: ft_bwd time: "
  //      << tm_ft.real() << endl;
  // cout << "SlaterDet: compute_density: rhosum time: "
  //      << tm_rhosum.real() << endl;
}

#ifndef _MPBC  /* regular PBC */
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::rs_mul_add(FourierTransform& ft,
  const double* v, SlaterDet& sdp) const
{
  // transform states to real space, multiply states by v[r] in real space
  // transform back to reciprocal space and add to sdp
  // sdp[n] += v * sd[n]

  vector<complex<double> > tmp(ft.np012loc());
  vector<complex<double> > ctmp(2*c_.mloc());

  const int np012loc = ft.np012loc();
  const int mloc = c_.mloc();
  double* p = (double*) &tmp[0];
  double* dcp = (double*) sdp.c().valptr();

  //clog << "rs_mul_add: np012loc = "<<np012loc<<"  c_.mloc = "<<mloc<<endl;
  //clog << "gvector: " << endl;
  //for ( int i = 0; i < np012; i++ ) 
  //  clog << basis_->gx(i) <<","<< basis_->gx(i+np012) <<","<< basis_->gx(i+np012*2) <<endl;


  if ( basis_->real() )
  {
#ifdef _NOANTIALIASING
    assert(0); //do not allow real_ basis if _NOANTIALIASING
#endif
    // transform two states at a time
    for ( int n = 0; n < nstloc()-1; n++, n++ )
    {
#ifndef _CHANGE_FT_CLASS
      /* functional call using old FT class */
      ft.backward(c_.cvalptr(n*mloc), c_.cvalptr((n+1)*mloc),&tmp[0]);
#else
      // not implemented
      assert(0);
#endif

      int ii = 0;
      for ( int i = 0; i < np012loc; i++ )
      {
        const double psi1 = p[ii];
        const double psi2 = p[ii+1];
        const double vii = v[i];
        p[ii]   = vii * psi1;
        p[ii+1] = vii * psi2;
        ii++; ii++;
      }
#ifndef _CHANGE_FT_CLASS
      /* functional call using old FT class */
      ft.forward(&tmp[0], &ctmp[0], &ctmp[mloc]);
#else
      // not implemented
      assert(0);
#endif

      int len = 4 * mloc;
      int inc1 = 1;
      double alpha = 1.0;
      daxpy(&len,&alpha,(double*)&ctmp[0],&inc1,&dcp[2*n*mloc],&inc1);
    }
    if ( nstloc() % 2 != 0 )
    {
      const int n = nstloc()-1;
#ifndef _CHANGE_FT_CLASS
      /* functional call using old FT class */
      ft.backward(c_.cvalptr(n*mloc),&tmp[0]);
#else
      // not implemented
      assert(0);
#endif
      int ii = 0;
      for ( int i = 0; i < np012loc; i++ )
      {
        const double psi1 = p[ii];
        const double vii = v[i];
        p[ii]   = vii * psi1;
        p[ii+1] = 0.0;
        ii++; ii++;
      }
#ifndef _CHANGE_FT_CLASS
      /* functional call using old FT class */
      ft.forward(&tmp[0], &ctmp[0]);
#else
      // not implemented
      assert(0);
#endif
      int len = 2 * mloc;
      int inc1 = 1;
      double alpha = 1.0;
      daxpy(&len,&alpha,(double*)&ctmp[0],&inc1,&dcp[2*n*mloc],&inc1);
    }
  }
  else
  {
#ifndef _REGULARGDISTRIB
    // original Qbox code
    // only one transform at a time
    for ( int n = 0; n < nstloc(); n++ )
    {
      ft.backward(c_.cvalptr(n*mloc),&tmp[0]);

      typedef std::complex<double> complex;
      using boost::multi_array_ref;
      using boost::extents;

      int np0=ft.np0();
      int np1=ft.np1();
      int np2=ft.np2();
      clog << "ft.np0*np1*np2 = "<<np0*np1*np2<<" np012loc = "<<np012loc<<endl;
      // save wf for HDF5 or Matlab
      if(boost::filesystem::exists("print_wf"))
      {
        multi_array_ref<complex, 3> wfxyz (&tmp[0], extents[np0][np1][np2]);
#ifdef _USE_HDF5
        Write_WF_to_HDF5(basis_,wfxyz,"wf_real",n);
#else
//        Write_WF_to_Matlab(wfxyz,"wf_real",n);
#endif
      }

      for ( int i = 0; i < np012loc; i++ )
        tmp[i] *= v[i];
      ft.forward(&tmp[0], &ctmp[0]);
      int len = 2 * mloc;
      int inc1 = 1;
      double alpha = 1.0;
      daxpy(&len,&alpha,(double*)&ctmp[0],&inc1,&dcp[2*n*mloc],&inc1);
    }
#else

    // regular distribution of G vectors

#ifndef _NOANTIALIASING
    BOOST_STATIC_ASSERT(0); //antialiasing must be disabled to use regular distribution of G vectors
#endif

    // use FFTW3 transforms, by-passing qbox rods design
    std::vector<double> enl_en(nstloc(),0);
    std::vector<double> pot_en(nstloc(),0);
    std::vector<double> kin_en(nstloc(),0);
    std::vector<double> tot_en(nstloc(),0);

    typedef std::complex<double> complex;
    using boost::const_multi_array_ref;
    using boost::multi_array_ref;
    using boost::extents;
    using fftw3::plan;

    int np0=basis_->np(0);
    int np1=basis_->np(1);
    int np2=basis_->np(2);
    int np012 = np0*np1*np2;

    double wfxyz_2sum, wfXYZ_2sum;

    //clog << "rs_mul_add: np0 = "<<np0<<" np1 = "<<np1<<" np2 = "<<np2<<endl;

    for ( int n = 0; n < nstloc(); n++ )
    {
      // bind c_ with multi_array_ref
      const_multi_array_ref<complex, 3> wfXYZ (c_.cvalptr(n*mloc),        extents[np0][np1][np2]);
      multi_array_ref<complex, 3>       wfxyz (&tmp[0],                   extents[np0][np1][np2]);
      multi_array_ref<complex, 3>       dwfXYZ(&ctmp[0],                  extents[np0][np1][np2]);
      multi_array_ref<complex, 3>       dcpXYZ((complex *)&dcp[2*n*mloc], extents[np0][np1][np2]);
      const_multi_array_ref<double,  3> vxyz  (v,                         extents[np0][np1][np2]);

      plan wf_backward(wfXYZ, wfxyz, FFTW_BACKWARD, FFTW_ESTIMATE);
      plan wf_forward (wfxyz,dwfXYZ, FFTW_FORWARD,  FFTW_ESTIMATE);

      //execute (wf_backward);
      wf_backward();

      if(boost::filesystem::exists("print_wf"))
      {
#ifdef _USE_HDF5
        Write_WF_to_HDF5(basis_,wfXYZ,"wf_rcpr",n);
        Write_WF_to_HDF5(basis_,wfxyz,"wf_real",n);
#else
        Write_WF_to_Matlab(wfXYZ,"wf_rcpr",n);
        Write_WF_to_Matlab(wfxyz,"wf_real",n);
#endif
      }

      // compute pot_en
      wfxyz_2sum = 0;
      for ( int i = 0; i < np012; i++ ) wfxyz_2sum +=norm(tmp[i]);
      pot_en[n] = 0;
      for ( int i = 0; i < np012; i++ ) pot_en[n]+=norm(tmp[i])*v[i];
      pot_en[n] /= wfxyz_2sum;

      // compute kin_en
      wfXYZ_2sum = 0;
      for ( int i = 0; i < np012; i++ ) wfXYZ_2sum +=norm(c_.cvalptr(n*mloc)[i]);
      kin_en[n] = 0;
      for ( int i = 0; i < np012; i++ ) kin_en[n]+=norm(c_.cvalptr(n*mloc)[i])*basis_->kpg2_PBC(i)*0.5;
      kin_en[n] /= wfXYZ_2sum;

      // compute non-local potential energy
      enl_en[n] = 0;
      for ( int i = 0; i < np012; i++ ) 
      { complex dpsi(dcp[2*n*mloc+2*i],dcp[2*n*mloc+2*i+1]);
        enl_en[n]+=real(conj(dpsi) * c_.cvalptr(n*mloc)[i]);
      }
      enl_en[n] /= wfXYZ_2sum;
      enl_en[n] -= kin_en[n];

      // apply real space potential v[i]
      for ( int i = 0; i < np012; i++ ) tmp[i] *= v[i] / np012;

      //execute (wf_forward);
      wf_forward();

      for ( int i = 0; i < np012; i++ ) 
      {
         dcp[2*n*mloc+2*i  ] += real(ctmp[i]);
         dcp[2*n*mloc+2*i+1] += imag(ctmp[i]);
      }

      // compute total energy
      tot_en[n] = 0;
      for ( int i = 0; i < np012; i++ ) 
      { complex dpsi(dcp[2*n*mloc+2*i],dcp[2*n*mloc+2*i+1]);
        tot_en[n]+=real(conj(dpsi) * c_.cvalptr(n*mloc)[i]);
      }
      tot_en[n] /= wfXYZ_2sum;
    }

    // print results to screen
    clog<<"FFTW3 enl_en: "; 
    for(int n=0; n!=enl_en.size(); n++)
      clog<<setw(12)<<setprecision(5)<<enl_en[n]<<" "; 
    clog<<endl;

    clog<<"FFTW3 pot_en: "; 
    for(int n=0; n!=pot_en.size(); n++)
      clog<<setw(12)<<setprecision(5)<<pot_en[n]<<" "; 
    clog<<endl;

    clog<<"FFTW3 kin_en: "; 
    for(int n=0; n!=kin_en.size(); n++)
      clog<<setw(12)<<setprecision(5)<<kin_en[n]<<" "; 
    clog<<endl;

//    clog<<"FFTW3 tot_en: "; 
//    for(int n=0; n!=pot_en.size(); n++)
//      clog<<setw(12)<<setprecision(5)<<enl_en[n]+pot_en[n]+kin_en[n]<<" "; 
//    clog<<" (from sum)"<<endl;

    clog<<"FFTW3 tot_en: "; 
    for(int n=0; n!=pot_en.size(); n++){
      clog<<setw(12)<<setprecision(5)<<tot_en[n]<<" "; 
			if(tot_en[n]!=tot_en[n]) throw std::runtime_error("catastrophic error, number is nan"); 
		}
    clog<<" (from <psi|H|psi>)"<<endl;

#endif
  }
}
#endif

#ifdef _MPBC /* Magnetic Periodic Boundary Condition (MPBC) */

void SlaterDet::compute_J_mpbc(FourierTransform& ft, int const state_idx) const
{
  using namespace boost;
  using fftw3::plan;
  typedef std::complex<double> complex;

  boost::array<int, 3> np; 
  np[0]=basis_->np(0);
  np[1]=basis_->np(1);
  np[2]=basis_->np(2);

  boost::array<multi_array<complex, 3> ,3> Jreal;
  Jreal[0].resize(extents[np[0]][np[1]][np[2]]);
  Jreal[1].resize(extents[np[0]][np[1]][np[2]]);
  Jreal[2].resize(extents[np[0]][np[1]][np[2]]);

  const_multi_array_ref<complex, 2> wfXZh(c_.cvalptr(state_idx*c_.mloc()), extents[np[0]][np[1]*np[2]]);
  multi_array<complex, 2>           wfXzh(extents[np[0]][np[1]*np[2]]);
  multi_array_ref<complex, 3>       wfXYz(wfXzh.origin(), extents[np[0]][np[1]][np[2]]);
  multi_array_ref<complex, 3>       wfxyz=wfXYz;

  multi_array<complex, 2>           Jz_wfXZh(extents[np[0]][np[1]*np[2]]);
  multi_array<complex, 2>           Jz_wfXzh(extents[np[0]][np[1]*np[2]]);;
  multi_array_ref<complex, 3>       Jz_wfXYz(Jz_wfXzh.origin(), extents[np[0]][np[1]][np[2]]);
  multi_array<complex, 3>           Jz_wfxyz(extents[np[0]][np[1]][np[2]]); //=Jz_wfXYz;

  multi_array<complex, 2>           Jx_wfXzh(extents[np[0]][np[1]*np[2]]);
  multi_array_ref<complex, 3>       Jx_wfXYz(Jx_wfXzh.origin(), extents[np[0]][np[1]][np[2]]);
  multi_array<complex, 3>           Jx_wfxyz(extents[np[0]][np[1]][np[2]]); //=Jx_wfXYz;

  multi_array<complex, 2>           Jy_wfXzh(extents[np[0]][np[1]*np[2]]);
  multi_array_ref<complex, 3>       Jy_wfXYz(Jy_wfXzh.origin(), extents[np[0]][np[1]][np[2]]);
  multi_array<complex, 3>           Jy_wfxyz(extents[np[0]][np[1]][np[2]]); //=Jy_wfXYz;

  const double * const gx    = basis_->gx_ptr(0);
  const double * const zhat  = basis_->zhat_ptr();
  const double * const gzhat = basis_->gzhat_ptr();

	{std::ofstream ofs("compute_J_z.dat");
  for(unsigned i=0; i<Jz_wfXZh.num_elements(); ++i){
		ofs<<gzhat[i]<<std::endl;
    Jz_wfXZh.origin()[i]=gzhat[i]*wfXZh.origin()[i];
	}}
  //plan(-1, wfXZh   , wfXzh   , fftw3::backward).execute_normalize();
  //plan(-1, Jz_wfXZh, Jz_wfXzh, fftw3::backward).execute_normalize();
  ft.z_backward_mpbc(wfXZh,    wfXzh);
  ft.z_backward_mpbc(Jz_wfXZh, Jz_wfXzh);
	{std::ofstream ofs("compute_J_x.dat");
  for(unsigned i=0; i<Jx_wfXzh.num_elements();++i){
		ofs<<gx[i+0*Jx_wfXzh.num_elements()]<<endl;
    Jx_wfXzh.origin()[i]=gx[i+0*Jx_wfXzh.num_elements()]*wfXzh.origin()[i];
	}}
	{std::ofstream ofs("compute_J_y.dat");
  for(unsigned i=0; i<Jy_wfXzh.num_elements();++i){
		ofs<<zhat[i]<<endl;
    Jy_wfXzh.origin()[i]=zhat[i]*wfXzh.origin()[i];
	}}
  //plan(+2, wfXYz   , wfxyz   , fftw3::backward).execute_normalize();
  //plan(+2, Jx_wfXYz, Jx_wfxyz, fftw3::backward).execute_normalize();
  //plan(+2, Jy_wfXYz, Jy_wfxyz, fftw3::backward).execute_normalize();
  //plan(+2, Jz_wfXYz, Jz_wfxyz, fftw3::backward).execute_normalize();
  ft.xy_backward_mpbc(wfXYz, wfxyz);
  ft.xy_backward_mpbc(Jx_wfXYz, Jx_wfxyz);
  ft.xy_backward_mpbc(Jy_wfXYz, Jy_wfxyz);
  ft.xy_backward_mpbc(Jz_wfXYz, Jz_wfxyz);
	
  //compute J
  for(int i=0; i<np[0]*np[1]*np[2]; ++i){
    Jreal[0].origin()[i]=real(conj(wfxyz.origin()[i])*Jx_wfxyz.origin()[i]);
    Jreal[1].origin()[i]=real(conj(wfxyz.origin()[i])*Jy_wfxyz.origin()[i]);
    Jreal[2].origin()[i]=real(conj(wfxyz.origin()[i])*Jz_wfxyz.origin()[i]);
  }


#ifdef _USE_HDF5
//  std::clog<<"printing J to HDF5 file "<<std::endl;
//  Write_WF_to_HDF5(basis_,Jreal[0],"Jx_real",state_idx);
//  Write_WF_to_HDF5(basis_,Jreal[1],"Jy_real",state_idx);
//  Write_WF_to_HDF5(basis_,Jreal[2],"Jz_real",state_idx);
  boost::multi_array<double, 4> Jreal_cont(extents[np[0]][np[1]][np[2]][3]);
  for(int i=0; i!=np[0]; ++i)
	  for(int j=0; j!=np[1]; ++j)
		  for(int k=0; k!=np[2]; ++k)
			  for(int l=0; l!=3; ++l){
				  Jreal_cont[i][j][k][l]=real(Jreal[l][i][j][k]);
			  }
  hdf5::file("Jreal_cont_"+boost::lexical_cast<std::string>(state_idx)+".hdf5")["array"]<<Jreal_cont;
	//	<<hdf5::make_nvp("array", Jreal_cont);
#else
  std::clog<<"printing J to Matlab file "<<std::endl;
  Write_WF_to_Matlab(wfxyz,"wf_real",state_idx);
  Write_WF_to_Matlab(Jreal[0],"Jx_real",state_idx);
  Write_WF_to_Matlab(Jreal[1],"Jy_real",state_idx);
  Write_WF_to_Matlab(Jreal[2],"Jz_real",state_idx);
#endif

  std::clog<<"computing magnetic field of current of state "<<state_idx<<std::endl;
  
  boost::array<multi_array<complex, 3> ,3> Jrcpr;
  Jrcpr[0].resize(extents[np[0]][np[1]][np[2]]);
  Jrcpr[1].resize(extents[np[0]][np[1]][np[2]]);
  Jrcpr[2].resize(extents[np[0]][np[1]][np[2]]);
 
  plan(Jreal[0], Jrcpr[0], fftw3::forward).execute_normalize();
  plan(Jreal[1], Jrcpr[1], fftw3::forward).execute_normalize();
  plan(Jreal[2], Jrcpr[2], fftw3::forward).execute_normalize();
  //ft.forward_pbc(Jreal[0], Jrcpr[0]);
  //ft.forward_pbc(Jreal[1], Jrcpr[1]);
  //ft.forward_pbc(Jreal[2], Jrcpr[2]);
  
  //In[181]:= Cross[{gx, gy, gz}, {Jx, Jy, Jz}]
  //Out[181]= {-gz Jy + gy Jz, gz Jx - gx Jz, -gy Jx + gx Jy}
  
  const_multi_array_ref<double, 1> gy(basis_->gx_ptr(1), extents[np[0]*np[1]*np[2]]);
  const_multi_array_ref<double, 1> gz(basis_->gx_ptr(2), extents[np[0]*np[1]*np[2]]);
  const_multi_array_ref<double, 1> g2(basis_->g2_ptr(),  extents[np[0]*np[1]*np[2]]);
  
  boost::array<multi_array<complex,3>, 3> Brcpr;
  Brcpr[0].resize(extents[np[0]][np[1]][np[2]]);
  Brcpr[1].resize(extents[np[0]][np[1]][np[2]]);
  Brcpr[2].resize(extents[np[0]][np[1]][np[2]]);

  complex I(0.0,1.0);
  for(unsigned i=0; i!=Brcpr[0].num_elements(); ++i)
  {
    if(g2[i]==0) continue;
    Brcpr[0].origin()[i]=I*((gy[i]*Jrcpr[2].origin()[i])-(gz[i]*Jrcpr[1].origin()[i]))/g2[i];
    Brcpr[1].origin()[i]=I*((gz[i]*Jrcpr[0].origin()[i])-(gx[i]*Jrcpr[2].origin()[i]))/g2[i];
    Brcpr[2].origin()[i]=I*((gx[i]*Jrcpr[1].origin()[i])-(gy[i]*Jrcpr[0].origin()[i]))/g2[i]; //there was an error here
  } 

  boost::array<multi_array<complex,3>, 3> Breal;
  Breal[0].resize(extents[np[0]][np[1]][np[2]]);
  Breal[1].resize(extents[np[0]][np[1]][np[2]]);
  Breal[2].resize(extents[np[0]][np[1]][np[2]]);

  plan(Brcpr[0], Breal[0], fftw3::backward).execute_normalize();
  plan(Brcpr[1], Breal[1], fftw3::backward).execute_normalize();
  plan(Brcpr[2], Breal[2], fftw3::backward).execute_normalize();
  //ft.backward_pbc(Brcpr[0], Breal[0]);
  //ft.backward_pbc(Brcpr[1], Breal[1]);
  //ft.backward_pbc(Brcpr[2], Breal[2]);
  

#ifdef _USE_HDF5
  std::clog<<"printing B to HDF5 file "<<std::endl;
//  Write_WF_to_HDF5(basis_,Breal[0],"Bx_real",state_idx);
//  Write_WF_to_HDF5(basis_,Breal[1],"By_real",state_idx);
//  Write_WF_to_HDF5(basis_,Breal[2],"Bz_real",state_idx);
	boost::multi_array<double, 3> Breal_cont(extents[np[1]][np[2]][2]);
	//for(int i=0; i!=np[0]; ++i)
		for(int j=0; j!=np[1]; ++j)
			for(int k=0; k!=np[2]; ++k)
				for(int l=0; l!=2; ++l){
					Breal_cont[j][k][l]=real(Breal[l][j][k][22]);
				}
	hdf5::file("Breal_cont_"+boost::lexical_cast<std::string>(state_idx)+".hdf5")["array"]<<Breal_cont;
	//	<<hdf5::make_nvp("array", Breal_cont);
#else
	std::clog<<"printing B to Matlab file "<<std::endl;
	Write_WF_to_Matlab(Breal[0],"Bx_real",state_idx);
	Write_WF_to_Matlab(Breal[1],"By_real",state_idx);
	Write_WF_to_Matlab(Breal[2],"Bz_real",state_idx);
#endif        
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::anti_aliasing_mpbc(FourierTransform& /*ft*/)
{
  if ( basis_->real() )
    assert(0); // basis has to be complex

  // only single CPU for now
  int np0=basis_->np(0);
  int np1=basis_->np(1);
  int np2=basis_->np(2);
  //int np012 = np0*np1*np2;
  const int mloc = c_.mloc();

  typedef std::complex<double> complex;

  using boost::multi_array_ref;
  using boost::multi_array;
  using boost::extents;

  for ( int n = 0; n < nstloc(); n++ )
  {
    multi_array_ref<complex, 3> wfXYZ(c_.valptr(n*mloc), extents[np0][np1][np2]);
    multi_array_ref<complex, 2> wfXZh(wfXYZ.origin()   , extents[np0][np1* np2]);
    multi_array<complex, 2>     wfXzh (                  extents[np0][np1* np2]);

#if 0 /* do intermediate space anti-aliasing only in rs_mul_add_mpbc */

    ft.z_backward_mpbc(wfXZh, wfXzh);

    /* anti-aliasing */
    for ( int i = 0; i < np0; i++)
      for ( int j = 0; j < np1*np2; j++)
      {
         if ( ( i >= floor(np0/4.0)     ) && ( i <= ceil(np0*3/4.0)     ) ) wfXzh[i][j] = 0;
         if ( ( j >= floor(np1*np2/4.0) ) && ( j <= ceil(np1*np2*3/4.0) ) ) wfXzh[i][j] = 0;
      }
    ft.z_forward_mpbc(wfXzh,wfXZh);

#endif

    /* anti-aliasing */
    for ( int i = 0; i < np0; i++)
      for ( int j = 0; j < np1*np2; j++)
      {
         if ( ( i >= floor(np0/4.0)     ) && ( i <= ceil(np0*3/4.0)     ) ) wfXZh[i][j] = 0;
#ifndef _SHORTFFT
         if ( ( j >= floor(np1*np2/4.0) ) && ( j <= ceil(np1*np2*3/4.0) ) ) wfXZh[i][j] = 0;
#else
         if ( ( j >= floor(np1*np2/8.0) ) && ( j <= ceil(np1*np2*7/8.0) ) ) wfXZh[i][j] = 0;
#endif
      }

    clog<<"calling ft.z_backward_mpbc from anti_aliasing_mpbc"<<endl;
    FourierTransform::z_backward_mpbc(wfXZh, wfXzh);
    FourierTransform::z_forward_mpbc (wfXzh, wfXZh);

  }
}
////////////////////////////////////////////////////////////////////////////////
void SlaterDet::rs_mul_add_mpbc(
	FourierTransform& ft,
	const double* v, SlaterDet& sdp) const
{
  // transform states to real space, multiply states by v[r] in real space
  // transform back to reciprocal space and add to sdp
  // sdp[n] += v * sd[n]

#ifndef _NOANTIALIASING
  BOOST_STATIC_ASSERT(0);
#endif

#ifndef _REGULARGDISTRIB
  BOOST_STATIC_ASSERT(0);
#endif

  if ( basis_->real() )
    assert(0); // basis has to be complex

  // only single CPU for now
  int np0=basis_->np(0);
  int np1=basis_->np(1);
  int np2=basis_->np(2);
  int np012 = np0*np1*np2;
  const int mloc = c_.mloc();
  double* dcp = (double*) sdp.c().valptr();
  const double *const gx  = basis_->gx_ptr(0);
  const double *const kpg2  = basis_->kpg2_MPBC_ptr();
  const double *const zhat2  = basis_->zhat2_ptr();
  const double *const gzhat  = basis_->gzhat_ptr();

  std::vector<double> enl_en(nstloc(),0);
  std::vector<double> pot_en(nstloc(),0);
  std::vector<double> kin_en_xz(nstloc(),0);
  std::vector<double> kin_en_x (nstloc(),0);
  std::vector<double> kin_en_y (nstloc(),0);
  std::vector<double> kin_en_z (nstloc(),0);

  typedef std::complex<double> complex;

  using boost::const_multi_array_ref;
  using boost::multi_array_ref;
	using boost::multi_array;
	
  using boost::extents;
  using fftw3::plan;

  double 
//	wfxyz_2sum, 
	wfXZh_2sum, 
	wfXzh_2sum;

  multi_array<complex, 2>           dwfXZh(                    extents[np0][np1* np2]);
  multi_array<complex, 2> 	    wfXzh (                    extents[np0][np1* np2]);
  multi_array_ref<complex, 3>       wfXYz (wfXzh.origin()    , extents[np0][np1][np2]);
  multi_array<complex, 3>     	    wfxyz (	 	       extents[np0][np1][np2]);
  multi_array<complex, 2>           dwfXzh(                    extents[np0][np1* np2]);
  multi_array_ref<complex, 3>       dwfXYz(dwfXzh.origin()   , extents[np0][np1][np2]);

  for ( int n = 0; n < nstloc(); n++ )
  {
    // bind c_ with multi_array_ref
    const_multi_array_ref<complex, 3> wfXYZ(c_.cvalptr(n*mloc), extents[np0][np1][np2]);
    const_multi_array_ref<complex, 2> wfXZh(wfXYZ.origin()    , extents[np0][np1* np2]);

    wfXZh_2sum = 0;
    for ( unsigned i = 0; i < wfXZh.num_elements(); i++ ) wfXZh_2sum +=norm(wfXZh.origin()[i]);
	
    ft.z_backward_mpbc(wfXZh, wfXzh);

#if 0 /* debug: see if z_backward_mpbc and z_forward_mpbc are reversible */
    multi_array<complex, 2> wfXZh_tmp (extents[np0][np1* np2]);
    ft.z_forward_mpbc(wfXzh, wfXZh_tmp);
    for (int i=0; i<np0; i++)
        for (int j=0; j<np1*np2; j++)
        {
            clog<<"wfXZh["<<i<<"]["<<j<<"]="<<wfXZh[i][j]<<"  "<<wfXZh_tmp[i][j]<<endl;
        }
#endif

    FourierTransform::xy_backward_mpbc(wfXYz, wfxyz);
		
    if(boost::filesystem::exists("print_wf"))
    {
#ifdef _USE_HDF5
//         Write_WF_to_HDF5(basis_,wfXYZ,"wf_rcpr",n);
//         Write_WF_to_HDF5(basis_,wfXYz,"wf_intr",n);
//         Write_WF_to_HDF5(basis_,wfxyz,"wf_real",n);
#else
//         Write_WF_to_Matlab(wfXYZ,"wf_rcpr",n);
//         Write_WF_to_Matlab(wfXYz,"wf_intr",n);
//         Write_WF_to_Matlab(wfxyz,"wf_real",n);
#endif
    }
	
    if(boost::filesystem::exists("print_J")){
			compute_J_mpbc(ft,n);
    }

    // compute pot_en
    {   pot_en[n] = 0;
        for ( unsigned i = 0; i < wfxyz.num_elements(); i++ ) pot_en[n]+=norm(wfxyz.origin()[i])*v[i];
    }

    // compute kin_en_x, kin_en_z
 	
    kin_en_xz[n] = 0; kin_en_x[n] = 0;  kin_en_z[n] = 0;
    for ( int i = 0; i < np012; i++ )
    {
       kin_en_xz[n]+=norm(c_.cvalptr(n*mloc)[i])*kpg2[i]*0.5;
       kin_en_x[n] +=norm(c_.cvalptr(n*mloc)[i])*gx[i]*gx[i]*0.5;
       kin_en_z[n] +=norm(c_.cvalptr(n*mloc)[i])*gzhat[i]*gzhat[i]*0.5;
    }
    kin_en_xz[n] /= wfXZh_2sum;
    kin_en_x [n] /= wfXZh_2sum;
    kin_en_z [n] /= wfXZh_2sum;

    // compute non-local potential energy
    enl_en[n] = 0;
    for ( int i = 0; i < np012; i++ ) 
    { complex dpsi(dcp[2*n*mloc+2*i],dcp[2*n*mloc+2*i+1]);
      enl_en[n]+=real(conj(dpsi) * c_.cvalptr(n*mloc)[i]);
    }
    enl_en[n] /= wfXZh_2sum;
    enl_en[n] -= kin_en_xz[n];

    // compute kin_en_y
    wfXzh_2sum = 0;
    for ( unsigned i = 0; i < wfXzh.num_elements(); i++ ) wfXzh_2sum +=norm(wfXzh.origin()[i]);
    kin_en_y[n] = 0;
    for ( unsigned i = 0; i < wfXzh.num_elements(); i++ ) kin_en_y[n]+=norm(wfXzh.origin()[i])*zhat2[i]* (0.5); // test
    kin_en_y[n] /= wfXzh_2sum;

    for ( int i = 0; i < np012; i++ ) wfxyz.origin()[i] *= v[i];

    /* forward 2D FFT along x and y directions: wfxyz (tmp_r) -> dwfXYz (tmp_i) */
    /* 3D array, transform first two dimensions (+2) */
    FourierTransform::xy_forward_mpbc(wfxyz,dwfXYz);

    /* apply (zhat)^2 potential to wfXzh (tmp, intermediate space) and add to dwfXzh (tmp_i) */
    for ( int i = 0; i < np012; i++ ) dwfXzh.origin()[i] += wfXzh.origin()[i] * zhat2[i] * (0.5); // test

#if 1
    /* anti-aliasing in intermediate space */
    for ( int i = 0; i < np0; i++)
      for ( int j = 0; j < np1*np2; j++)
      {
         if ( ( i >= floor(np0/4.0)     ) && ( i <= ceil(np0*3/4.0)     ) ) dwfXzh[i][j] = 0;
         if ( ( j >= floor(np1*np2/4.0) ) && ( j <= ceil(np1*np2*3/4.0) ) ) dwfXzh[i][j] = 0;
      }
#endif

    ft.z_forward_mpbc(dwfXzh,dwfXZh);

#if 1
    /* anti-aliasing in rcpr space */
    for ( int i = 0; i < np0; i++)
      for ( int j = 0; j < np1*np2; j++)
      {
         if ( ( i >= floor(np0/4.0)     ) && ( i <= ceil(np0*3/4.0)     ) ) dwfXZh[i][j] = 0;
#ifndef _SHORTFFT
         if ( ( j >= floor(np1*np2/4.0) ) && ( j <= ceil(np1*np2*3/4.0) ) ) dwfXZh[i][j] = 0;
#else
         if ( ( j >= floor(np1*np2/8.0) ) && ( j <= ceil(np1*np2*7/8.0) ) ) dwfXZh[i][j] = 0;
#endif
      }
#endif

    for ( unsigned i = 0; i < dwfXZh.num_elements(); i++ ) {
       dcp[2*n*mloc+2*i  ] += real(dwfXZh.origin()[i]);
       dcp[2*n*mloc+2*i+1] += imag(dwfXZh.origin()[i]);
    }
  }

  // print results to screen
  clog<<"FFTW3 enl_en: "; 
  for(unsigned n=0; n!=enl_en.size(); n++)
    clog<<setw(12)<<setprecision(5)<<enl_en[n]<<" "; 
  clog<<endl;

  clog<<"FFTW3 pot_en: "; 
  for(unsigned n=0; n!=pot_en.size(); n++)
    clog<<setw(12)<<setprecision(5)<<pot_en[n]<<" "; 
  clog<<endl;

  clog<<"FFTW3 kin_en_xz: "; 
  for(unsigned n=0; n!=kin_en_xz.size(); n++)
    clog<<setw(12)<<setprecision(5)<<kin_en_xz[n]<<" "; 
  clog<<endl;

  clog<<"FFTW3 kin_en_x: "; 
  for(unsigned n=0; n!=kin_en_y.size(); n++)
    clog<<setw(12)<<setprecision(5)<<kin_en_x[n]<<" "; 
  clog<<endl;

  clog<<"FFTW3 kin_en_y: "; 
  for(unsigned n=0; n!=kin_en_y.size(); n++)
    clog<<setw(12)<<setprecision(5)<<kin_en_y[n]<<" "; 
  clog<<endl;

  clog<<"FFTW3 kin_en_z: "; 
  for(unsigned n=0; n!=kin_en_y.size(); n++)
    clog<<setw(12)<<setprecision(5)<<kin_en_z[n]<<" "; 
  clog<<endl;

  clog<<"FFTW3 tot_en: "; 
  for(unsigned n=0; n!=pot_en.size(); n++){
		double etot = pot_en[n]+kin_en_xz[n]+kin_en_y[n];
    clog<<setw(12)<<setprecision(12)<<etot<<" "; 
		if(etot!=etot) throw std::runtime_error("Catastrophic error, number is nan"); 
	}
  clog<<endl;

}
#endif

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::gram(void)
{
  if ( basis_->real() )
  {
    // k = 0 case
    // create a DoubleMatrix proxy for c_
    DoubleMatrix c_proxy(c_);
    DoubleMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    s.syrk('l','t',2.0,c_proxy,0.0);
    s.syr('l',-1.0,c_proxy,0,'r');
#ifdef CHOLESKY_REMAP
    // create a square context for the Cholesky decomposition
    // int nsq = (int) sqrt((double) ctxt_.size());
    int nsq = CHOLESKY_REMAP;
    Context csq(nsq,nsq);
    DoubleMatrix ssq(csq,c_.n(),c_.n(),c_.nb(),c_.nb());
    ssq.getsub(s,s.m(),s.n(),0,0);
    ssq.potrf('l'); // Cholesky decomposition: S = L * L^T
    s.getsub(ssq,s.m(),s.n(),0,0);
#else
    s.potrf('l'); // Cholesky decomposition: S = L * L^T
#endif
    // solve triangular system X * L^T = C
    c_proxy.trsm('r','l','t','n',1.0,s);
  }
  else
  {
    // k != 0 case
    ComplexMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    s.herk('l','c',1.0,c_,1.0);
    s.potrf('l'); // Cholesky decomposition: S = L * L^H
    // solve triangular system X * L^H = C
    c_.trsm('r','l','c','n',1.0,s);
  }
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::riccati(const SlaterDet& sd)
{
  if ( basis_->real() )
  {
    // k = 0 case
    DoubleMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix r(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    s.identity();
    r.identity();

    DoubleMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix xm(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());

    // DoubleMatrix proxy for c_ and sd.c()
    DoubleMatrix c_proxy(c_);
    DoubleMatrix sdc_proxy(sd.c());

    // Factor -1.0 in next line: -0.5 from definition of s, 2.0 for G and -G
    s.syrk('l','t',-1.0,c_proxy,0.5); // s = 0.5 * ( I - A )
    // symmetric rank-1 update using first row of c_proxy
    s.syr('l',0.5,c_proxy,0,'r');
    // factor -2.0 in next line: G and -G
    r.gemm('t','n',-2.0,sdc_proxy,c_proxy,1.0); // r = ( I - B )
    // rank-1 update using first row of sdc_proxy() and c_proxy
    r.ger(1.0,sdc_proxy,0,c_proxy,0);

    xm = s;
    xm.symmetrize('l');

    s.syrk('l','t',0.5,r,1.0); // s = s + 0.5 * r^T * r
    s.symmetrize('l');

    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;

    while ( iter < maxiter && diff > epsilon )
    {
      // x = s - 0.5 * ( r - xm )^T * ( r - xm )
      // Note: t and r are not symmetric, x, xm, and s are symmetric

      for ( unsigned i = 0; i < t.size(); i++ )
        t[i] = r[i] - xm[i];

      x = s;
      x.syrk('l','t',-0.5,t,1.0);

      // get full matrix x
      x.symmetrize('l');

      for ( unsigned i = 0; i < t.size(); i++ )
        t[i] = x[i] - xm[i];

      diff = t.nrm2();

      xm = x;
      iter++;
    }
    c_proxy.symm('r','l',1.0,x,sdc_proxy,1.0);
  }
  else
  {
    // k != 0 case
    ComplexMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix r(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    s.identity();
    r.identity();

    ComplexMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix xm(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    ComplexMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());

    // s = 0.5 * ( I - A )
    s.herk('l','c',-0.5,c_,0.5);
    // r = ( I - B )
    r.gemm('c','n',-1.0,sd.c(),c_,1.0);

    xm = s;
    xm.symmetrize('l');

    // s = s + 0.5 * r^H * r
    s.herk('l','c',0.5,r,1.0);
    s.symmetrize('l');

    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;

    while ( iter < maxiter && diff > epsilon )
    {
      // x = s - 0.5 * ( r - xm )^H * ( r - xm )
      // Note: t and r are not hermitian, x, xm, and s are hermitian

      for ( unsigned i = 0; i < t.size(); i++ )
        t[i] = r[i] - xm[i];

      x = s;
      x.herk('l','c',-0.5,t,1.0);
      x.symmetrize('l');

      for ( unsigned i = 0; i < t.size(); i++ )
        t[i] = x[i] - xm[i];

      diff = t.nrm2();

      xm = x;
      iter++;
    }
    c_.hemm('r','l',1.0,x,sd.c(),1.0);
  }
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::lowdin(void)
{
  // Higham algorithm for polar decomposition
  if ( basis_->real() )
  {
    ComplexMatrix c_tmp(c_);
    DoubleMatrix c_proxy(c_);
    DoubleMatrix c_tmp_proxy(c_tmp);
    DoubleMatrix l(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix xp(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());

    l.clear();
    l.syrk('l','t',2.0,c_proxy,0.0);
    l.syr('l',-1.0,c_proxy,0,'r');

    //cout << "SlaterDet::lowdin: A=\n" << l << endl;

    // Cholesky decomposition of A=Y^T Y
    l.potrf('l');
    // The lower triangle of l now contains the Cholesky factor of Y^T Y

    //cout << "SlaterDet::lowdin: L=\n" << l << endl;

    // Compute the polar decomposition of R = L^T

    x.transpose(1.0,l,0.0);
    // x now contains R
    //cout << "SlaterDet::lowdin: R=\n" << x << endl;

    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;

    while ( iter < maxiter && diff > epsilon )
    {
      // t = X^T
      t.transpose(1.0,x,0.0);
      t.inverse();

      // t now contains X^-T

      // xp = 0.5 * ( x + x^-T );
      for ( unsigned i = 0; i < x.size(); i++ )
        xp[i] = 0.5 * ( x[i] + t[i] );


      // Next lines: use t as temporary to compute || x - xp ||_F
      for ( unsigned i = 0; i < t.size(); i++ )
        t[i] = x[i] - xp[i];

      diff = t.nrm2();

      //cout << " SlaterDet::lowdin: diff=" << diff << endl;

      x = xp;
      //cout << "SlaterDet::lowdin: X=\n" << x << endl;

      iter++;
    }

    // x now contains the orthogonal polar factor U of the
    // polar decomposition R = UH

    //cout << " SlaterDet::lowdin: orthogonal polar factor=\n" << x << endl;

    // Compute L^-1
    l.trtri('l','n');
    // l now contains L^-1

    // Form the product L^-T U
    t.gemm('t','n',1.0,l,x,0.0);

    // Multiply c by L^-T U
    c_proxy.gemm('n','n',1.0,c_tmp_proxy,t,0.0);

  }
  else
  {
    // complex case: not implemented
    if ( ctxt_.onpe0() )
      cout << " SlaterDet::lowdin: not implemented, reverting to Gram-Schmidt"
           << endl;
    gram();
  }
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::ortho_align(const SlaterDet& sd)
{
  // Orthogonalize *this and align with sd
  // Higham algorithm for polar decomposition
  if ( basis_->real() )
  {
    ComplexMatrix c_tmp(c_);
    DoubleMatrix c_proxy(c_);
    DoubleMatrix sdc_proxy(sd.c());
    DoubleMatrix c_tmp_proxy(c_tmp);
    DoubleMatrix l(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix xp(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());

#if TIMING
    tmap["syrk"].reset();
    tmap["syrk"].start();
#endif
    l.clear();
    l.syrk('l','t',2.0,c_proxy,0.0);
    l.syr('l',-1.0,c_proxy,0,'r');
#if TIMING
    tmap["syrk"].stop();
#endif

    //cout << "SlaterDet::ortho_align: A=\n" << l << endl;

    // Cholesky decomposition of A=Y^T Y
#if TIMING
    tmap["potrf"].reset();
    tmap["potrf"].start();
#endif
    l.potrf('l');
#if TIMING
    tmap["potrf"].stop();
#endif
    // The lower triangle of l now contains the Cholesky factor of Y^T Y

    //cout << "SlaterDet::ortho_align: L=\n" << l << endl;

    // Compute the polar decomposition of L^-1 B
    // where B = C^T sd.C

    // Compute B: store result in x
#if TIMING
    tmap["gemm"].reset();
    tmap["gemm"].start();
#endif
    // factor -2.0 in next line: G and -G
    x.gemm('t','n',2.0,c_proxy,sdc_proxy,0.0);
    // rank-1 update using first row of sdc_proxy() and c_proxy
    x.ger(-1.0,c_proxy,0,sdc_proxy,0);
#if TIMING
    tmap["gemm"].stop();
#endif

    // Form the product L^-1 B, store result in x
    // triangular solve: L X = B
    // trtrs: solve op(*this) * X = Z, output in Z
#if TIMING
    tmap["trtrs"].reset();
    tmap["trtrs"].start();
#endif
    l.trtrs('l','n','n',x);
#if TIMING
    tmap["trtrs"].stop();
#endif
    // x now contains L^-1 B

    //cout << "SlaterDet::ortho_align: L^-1 B=\n" << x << endl;

    // compute the polar decomposition of L^-1 B
    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;

#if TIMING
    tmap["transpose"].reset();
    tmap["inverse"].reset();
#endif
    while ( iter < maxiter && diff > epsilon )
    {
      // t = X^T
#if TIMING
      tmap["transpose"].start();
#endif
      t.transpose(1.0,x,0.0);
#if TIMING
      tmap["transpose"].stop();
      tmap["inverse"].start();
#endif
      t.inverse();
#if TIMING
      tmap["inverse"].stop();
#endif

      // t now contains X^-T

      // xp = 0.5 * ( x + x^-T );
      for ( unsigned i = 0; i < x.size(); i++ )
        xp[i] = 0.5 * ( x[i] + t[i] );


      // Next lines: use t as temporary to compute || x - xp ||_F
      for ( unsigned i = 0; i < t.size(); i++ )
        t[i] = x[i] - xp[i];

#if TIMING
      tmap["nrm2"].start();
#endif
      diff = t.nrm2();
#if TIMING
      tmap["nrm2"].stop();
#endif

      //cout << " SlaterDet::ortho_align: diff=" << diff << endl;

      x = xp;
      //cout << "SlaterDet::ortho_align: X=\n" << x << endl;

      iter++;
    }

    // x now contains the orthogonal polar factor X of the
    // polar decomposition L^-1 B = XH

    //cout << " SlaterDet::ortho_align: orthogonal polar factor=\n"
    //     << x << endl;

    // Form the product L^-T Q
    // Solve trans(L) Z = X
#if TIMING
    tmap["trtrs2"].reset();
    tmap["trtrs2"].start();
#endif
    l.trtrs('l','t','n',x);
#if TIMING
    tmap["trtrs2"].stop();
#endif

    // x now contains L^-T Q

    // Multiply c by L^-T Q
#if TIMING
    tmap["gemm2"].reset();
    tmap["gemm2"].start();
#endif
    c_proxy.gemm('n','n',1.0,c_tmp_proxy,x,0.0);
#if TIMING
    tmap["gemm2"].stop();
#endif

  }
  else
  {
    // complex case: not implemented
    if ( ctxt_.onpe0() )
      cout << " SlaterDet::lowdin: not implemented, reverting to riccati"
           << endl;
    riccati(sd);
  }
#if TIMING
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    ctxt_.dmin(1,1,&tmin,1);
    ctxt_.dmax(1,1,&tmax,1);
    if ( ctxt_.onpe0() )
    {
      cout << "<timing name=\""
           << setw(15) << (*i).first << "\""
           << " min=\"" << setprecision(3) << setw(9) << tmin << "\""
           << " max=\"" << setprecision(3) << setw(9) << tmax << "\"/>"
           << endl;
    }
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::align(const SlaterDet& sd)
{
  // Align *this with sd
  // Higham algorithm for polar decomposition
  if ( basis_->real() )
  {
    ComplexMatrix c_tmp(c_);
    DoubleMatrix c_proxy(c_);
    DoubleMatrix sdc_proxy(sd.c());
    DoubleMatrix c_tmp_proxy(c_tmp);
    DoubleMatrix x(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix xp(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    DoubleMatrix t(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());


    // Compute the polar decomposition of B
    // where B = C^T sd.C

    // Compute B: store result in x
    // factor -2.0 in next line: G and -G
    x.gemm('t','n',2.0,c_proxy,sdc_proxy,0.0);
    // rank-1 update using first row of sdc_proxy() and c_proxy
    x.ger(-1.0,c_proxy,0,sdc_proxy,0);

    // x now contains B

    //cout << "SlaterDet::align: B=\n" << x << endl;

    // Compute the distance | c - sdc | before alignment
    //for ( int i = 0; i < c_proxy.size(); i++ )
    //  c_tmp_proxy[i] = c_proxy[i] - sdc_proxy[i];
    //cout << " SlaterDet::align: distance before: "
    //     << c_tmp_proxy.nrm2() << endl;

    // compute the polar decomposition of B
    double diff = 1.0;
    const double epsilon = 1.e-10;
    const int maxiter = 20;
    int iter = 0;

    while ( iter < maxiter && diff > epsilon )
    {
      // t = X^T
      t.transpose(1.0,x,0.0);
      t.inverse();

      // t now contains X^-T

      // xp = 0.5 * ( x + x^-T );
      for ( unsigned i = 0; i < x.size(); i++ )
        xp[i] = 0.5 * ( x[i] + t[i] );


      // Next lines: use t as temporary to compute || x - xp ||_F
      //for ( int i = 0; i < t.size(); i++ )
      //  t[i] = x[i] - xp[i];

      //diff = t.nrm2();

      //cout << " SlaterDet::align: diff=" << diff << endl;

      x = xp;
      //cout << "SlaterDet::align: X=\n" << x << endl;

      iter++;
    }

    // x now contains the orthogonal polar factor X of the
    // polar decomposition B = XH

    //cout << " SlaterDet::align: orthogonal polar factor=\n" << x << endl;

    // Multiply c by X
    c_tmp_proxy = c_proxy;
    c_proxy.gemm('n','n',1.0,c_tmp_proxy,x,0.0);

    // Compute the distance | c - sdc | after alignment
    //for ( int i = 0; i < c_proxy.size(); i++ )
    //  c_tmp_proxy[i] = c_proxy[i] - sdc_proxy[i];
    //cout << " SlaterDet::align: distance after:  "
    //     << c_tmp_proxy.nrm2() << endl;

  }
  else
  {
    // complex case: not implemented
    if ( ctxt_.onpe0() )
      cout << " SlaterDet::align: not implemented, alignment skipped"
           << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
complex<double> SlaterDet::dot(const SlaterDet& sd) const
{
  // dot product of Slater determinants at the same kpoint: dot = tr (V^T W)
  assert(basis_->kpoint() == sd.kpoint());
  if ( basis_->real() )
  {
    // DoubleMatrix proxy for c_ and sd.c()
    const DoubleMatrix c_proxy(c_);
    const DoubleMatrix sdc_proxy(sd.c());
    // factor 2.0: G and -G
    double d = 2.0 * c_proxy.dot(sdc_proxy);

    // correct double counting of first element
    double sum = 0.0;
    if ( ctxt_.myrow() == 0 )
    {
      // compute the scalar product of the first rows of c_ and sd.c_
      const double *c = c_proxy.cvalptr(0);
      const double *sdc = sdc_proxy.cvalptr(0);
      int len = c_proxy.nloc();
      // stride of scalar product is mloc
      int stride = c_proxy.mloc();
      sum = ddot(&len,c,&stride,sdc,&stride);
    }
    ctxt_.dsum(1,1,&sum,1);
    return d - sum;
  }
  else
  {
    return c_.dot(sd.c());
  }
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::update_occ(unsigned nel, unsigned nspin)
{
  // compute occupation numbers as 0.0, 1.0 or 2.0
  // if nspin = 1: use 0, 1 or 2
  // if nspin = 2: use 0 or 1;
  //assert (nel >= 0);
  assert (occ_.size() == c_.n());
  if ( nspin == 1 )
  {
    assert (nel <= 2*c_.n());
    int ndouble = nel/2;
    for ( int n = 0; n < ndouble; n++ )
      occ_[n] = 2.0;
    for ( unsigned n = ndouble; n < ndouble+nel%2; n++ )
      occ_[n] = 1.0;
    for ( unsigned n = ndouble+nel%2; n < c_.n(); n++ )
      occ_[n] = 0.0;
  }
  else if ( nspin == 2 )
  {
    assert (nel <= c_.n());
    for ( unsigned n = 0; n < nel; n++ )
      occ_[n] = 1.0;
    for ( unsigned n = nel; n < c_.n(); n++ )
      occ_[n] = 0.0;
  }
  else
  {
    // incorrect value of nspin_
    assert(false);
  }
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::total_charge(void) const
{
  // compute total charge from occ_[i]
  double sum = 0.0;
  for ( unsigned n = 0; n < occ_.size(); n++ )
  {
    sum += occ_[n];
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::update_occ(unsigned nspin, double mu, double temp)
{
  // compute occupation numbers using a Fermi distribution f(mu,temp)
  // and the eigenvalues in eig_[i]
  assert(nspin==1 || nspin==2);
  assert (occ_.size() == c_.n());
  assert (eig_.size() == c_.n());
  if ( nspin == 1 )
  {
    for ( unsigned n = 0; n < eig_.size(); n++ )
    {
      occ_[n] = 2.0 * fermi(eig_[n],mu,temp);
    }
  }
  else if ( nspin == 2 )
  {
    for ( unsigned n = 0; n < eig_.size(); n++ )
    {
      occ_[n] = fermi(eig_[n],mu,temp);
    }
  }
  else
  {
    // incorrect value of nspin_
    assert(false);
  }
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::fermi(double e, double mu, double fermitemp)
{
  // e, mu in Hartree, fermitemp in Kelvin

  if ( fermitemp == 0.0 )
  {
    if ( e < mu ) return 1.0;
    else if ( e == mu ) return 0.5;
    else return 0.0;
  }
  const double kb = 3.1667907e-6; // Hartree/Kelvin
  const double kt = kb * fermitemp;
  double arg = ( e - mu ) / kt;

  if ( arg < -30.0 ) return 1.0;
  if ( arg >  30.0 ) return 0.0;

  return 1.0 / ( 1.0 + exp ( arg ) );
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::entropy(int nspin) const
{
  // return dimensionless entropy
  // the contribution to the free energy is - t_kelvin * k_boltz * wf.entropy()

  assert(nspin==1 || nspin==2);
  const double fac = ( nspin > 1 ) ? 1.0 : 2.0;
  double sum = 0.0;
  for ( unsigned n = 0; n < occ_.size(); n++ )
  {
    const double f = occ_[n] / fac;
    if ( f > 0.0  &&  f < 1.0 )
    {
      sum -= fac * ( f * log(f) + (1.0-f) * log(1.0-f) );
    }
  }
  return sum;
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::ortho_error(void) const
{
  // deviation from orthogonality of c_
  double error;
  if ( basis_->real() )
  {
    // k = 0 case
    // declare a proxy DoubleMatrix for c_
    DoubleMatrix c_proxy(c_);

    DoubleMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());

    // real symmetric rank-k update
    // factor 2.0 in next line: G and -G
    s.syrk('l','t',2.0,c_proxy,0.0); // compute real overlap matrix

    // correct for double counting of G=0
    // symmetric rank-1 update using first row of c_proxy
    s.syr('l',-1.0,c_proxy,0,'r');

    DoubleMatrix id(ctxt_,s.m(),s.n(),s.mb(),s.nb());
    id.identity();

    s -= id; // subtract identity matrix from S

    error = s.nrm2();
  }
  else
  {
    // k != 0 case

    ComplexMatrix s(ctxt_,c_.n(),c_.n(),c_.nb(),c_.nb());
    s.herk('l','c',1.0,c_,0.0);

    ComplexMatrix id(ctxt_,s.m(),s.n(),s.mb(),s.nb());
    id.identity();

    s -= id; // subtract identity matrix from S

    error = s.nrm2();
  }
  return error;
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::randomize(double amplitude)
{
  if ( basis_->size() == 0 )
    return;
  // Note: randomization results depend on the process grid size and shape
  srand48(ctxt_.myproc());
  for ( int n = 0; n < c_.nloc(); n++ )
  {
    complex<double>* p = c_.valptr(c_.mloc()*n);
    for ( int i = 0; i < basis_->localsize(); i++ )
    {
      double re = drand48();
      double im = drand48();
      p[i] += amplitude * complex<double>(re,im);
    }
#ifdef _MPBC
#if 1
    // anti-aliasing
    int np0=basis_->np(0);
    int np1=basis_->np(1);
    int np2=basis_->np(2);
    for ( int i = 0; i < np0; i++ )
      for ( int j = 0; j < np1*np2; j++ )
      {
         if ( ( i >= floor(np0/4.0)     ) && ( i <= floor(np0*3/4.0)     ) ) p[i*np1*np2+j] = 0;
#ifndef _SHORTFFT
         if ( ( j >= floor(np1*np2/4.0) ) && ( j <= floor(np1*np2*3/4.0) ) ) p[i*np1*np2+j] = 0;
#else
         if ( ( j >= floor(np1*np2/8.0) ) && ( j <= floor(np1*np2*7/8.0) ) ) p[i*np1*np2+j] = 0;
#endif
      }
#endif
#endif
  }
  cleanup();
  gram();
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::cleanup(void)
{
  // set Im( c(G=0) ) to zero and
  // set the empty rows of the matrix c_ to zero
  // The empty rows are located between i = basis_->localsize() and
  // c_.mloc(). Empty rows are necessary to insure that the
  // local size c_.mloc() is the same on all processes, while the
  // local basis size is not.
  for ( int n = 0; n < c_.nloc(); n++ )
  {
    complex<double>* p = c_.valptr(c_.mloc()*n);
    // reset imaginary part of G=0 component to zero
    if ( ctxt_.myrow() == 0 )
    {
      // index of G=0 element
      int izero;
      if ( basis_->real() )
        izero = 0;
      else
        izero = basis_->rod_size(0)/2;
      //cout << " izero = " << izero << " G = " << basis_->kv(3*izero) << " "
      //     << basis_->kv(3*izero+1) << " " << basis_->kv(3*izero+2) << endl;
      p[izero] = complex<double> ( p[izero].real(), 0.0);
    }
    // reset values of empty rows of c_ to zero
    for ( int i = basis_->localsize(); i < c_.mloc(); i++ )
    {
      p[i] = 0.0;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
SlaterDet& SlaterDet::operator=(SlaterDet& rhs)
{
  if ( this == &rhs ) return *this;
  assert(ctxt_.ictxt() == rhs.context().ictxt());
  c_ = rhs.c_;
  return *this;
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::memsize(void) const
{
  return basis_->memsize() + c_.memsize();
}

////////////////////////////////////////////////////////////////////////////////
double SlaterDet::localmemsize(void) const
{
  return basis_->localmemsize() + c_.localmemsize();
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::print(ostream& os, string encoding, double weight, int ispin,
  int nspin) const
{	
#ifndef _CHANGE_FT_CLASS
  FourierTransform ft(*basis_,basis_->np(0),basis_->np(1),basis_->np(2));
#else
  FourierTransform ft(basis_->np());
#endif
	
  vector<complex<double> > wftmp(ft.np012loc());
  const bool real_basis = basis_->real();
  const int wftmpr_size = real_basis ? ft.np012() : 2*ft.np012();
  const int wftmpr_loc_size = real_basis ? ft.np012loc() : 2*ft.np012loc();
  vector<double> wftmpr(wftmpr_size);
  Base64Transcoder xcdr;

  if ( ctxt_.onpe0() )
  {
    string spin = (ispin > 0) ? "down" : "up";
    os << "<slater_determinant";
    if ( nspin == 2 )
      os << " spin=\"" << spin << "\"";
    os << " kpoint=\"" << basis_->kpoint() << "\"\n"
       << "  weight=\"" << setprecision(12) <<  weight << "\""
       << " size=\"" << nst() << "\">" << endl;

    os << "<density_matrix form=\"diagonal\" size=\"" << nst() << "\">"
       << endl;
    os.setf(ios::fixed,ios::floatfield);
    os.setf(ios::right,ios::adjustfield);
    for ( int i = 0; i < nst(); i++ )
    {
      os << " " << setprecision(8) << occ_[i];
      if ( i%10 == 9 )
        os << endl;
    }
    if ( nst()%10 != 0 )
      os << endl;
    os << "</density_matrix>" << endl;
  }

  for ( int n = 0; n < nst(); n++ )
  {
    // Barrier to limit the number of messages sent to task 0
    // that don't have a receive posted
    ctxt_.barrier();

    // transform data on ctxt_.mycol()
    if ( c_.pc(n) == ctxt_.mycol() )
    {
      //cout << " state " << n << " is stored on column "
      //     << ctxt_.mycol() << " local index: " << c_.y(n) << endl;
      int nloc = c_.y(n); // local index

#ifndef _REGULARGDISTRIB
      ft.backward(c_.cvalptr(c_.mloc()*nloc),&wftmp[0]);
#else
      // use FFTW3 transforms, by-passing qbox rods design
      typedef std::complex<double> complex;
      using boost::const_multi_array_ref;
      using boost::multi_array_ref;
      using boost::extents;
      using fftw3::plan;
      int np0=basis_->np(0); int np1=basis_->np(1); int np2=basis_->np(2); /*int np012 = np0*np1*np2;*/

#ifndef _MPBC 
      // 3D FFT in conventional PBC

      // use FFTW3 transforms, by-passing qbox rods design
      const_multi_array_ref<complex, 3> wfXYZ (c_.cvalptr(c_.mloc()*nloc), extents[np0][np1][np2]);
      multi_array_ref<complex, 3>       wfxyz (&wftmp[0],                 extents[np0][np1][np2]);

      plan wf_backward(wfXYZ, wfxyz, FFTW_BACKWARD, FFTW_ESTIMATE); wf_backward();
#else
      // FT from reciprocal space to intermediate space to real space in MPBC
      multi_array<complex, 2>           wfXzh (                           extents[np0][np1* np2]);
      multi_array_ref<complex, 3>       wfXYz (wfXzh.origin(),            extents[np0][np1][np2]);

      const_multi_array_ref<complex, 3> wfXYZ (c_.cvalptr(c_.mloc()*nloc),extents[np0][np1][np2]);
      const_multi_array_ref<complex, 2> wfXZh (c_.cvalptr(c_.mloc()*nloc),extents[np0][np1* np2]);
      multi_array_ref<complex, 3>       wfxyz (&wftmp[0],                 extents[np0][np1][np2]);

      /* backward 1D FFT along zhat direction: wfXZh -> wfXzh = wfXYz (tmp) */
      /* 2D array, transform last one dimension (-1) */
      //plan wf_z_backward(-1, wfXZh, wfXzh, FFTW_BACKWARD, FFTW_ESTIMATE);  wf_z_backward();
      /* normalization factor (1) */
      //for ( int i = 0; i < np012; i++ ) wfXzh.origin()[i] /= sqrt(np1*np2);

      /* backward 2D FFT along x and y directions: wfXYz (tmp) -> wfxyz (tmp_r) */
      /* 3D array, transform first two dimensions (+2) */
      //plan wf_xy_backward(+2, wfXYz, wfxyz, FFTW_BACKWARD, FFTW_ESTIMATE);  wf_xy_backward();
      /* normalization factor (2) */
      //for ( int i = 0; i < np012; i++ ) wftmp[i] /= sqrt(np0*np1);

      FourierTransform::z_backward_mpbc(wfXZh, wfXzh);
      FourierTransform::xy_backward_mpbc(wfXYZ, wfxyz);

#endif //_MPBC
#endif //_REGULARGDISTRIB


      if ( real_basis )
      {
#ifdef _NOANTIALIASING
        assert(0);
#endif
        double *a = (double*) &wftmp[0];
        for ( unsigned i = 0; i < ft.np012loc(); i++ )
          wftmpr[i] = a[2*i];
      }
      else
      {
#ifndef _REGULARGDISTRIB
        memcpy((void*)&wftmpr[0],(void*)&wftmp[0],
               ft.np012loc()*sizeof(complex<double>));
#else
        memcpy((void*)&wftmpr[0],(void*)&wftmp[0],
               ft.np012loc()*sizeof(complex));
#endif
      }
    }

    // send blocks of wftmpr to pe0
    for ( int i = 0; i < ctxt_.nprow(); i++ )
    {
      bool iamsending = c_.pc(n) == ctxt_.mycol() && i == ctxt_.myrow();

      // send size of wftmpr block
      int size=-1;
      if ( ctxt_.onpe0() )
      {
        if ( iamsending )
        {
          // sending to self, size not needed
        }
        else
          ctxt_.irecv(1,1,&size,1,i,c_.pc(n));
      }
      else
      {
        if ( iamsending )
        {
          size = wftmpr_loc_size;
          ctxt_.isend(1,1,&size,1,0,0);
        }
      }

      // send wftmpr block
      if ( ctxt_.onpe0() )
      {
        if ( iamsending )
        {
          // do nothing, data is already in place
        }
        else
        {
          int istart = ft.np0() * ft.np1() * ft.np2_first(i);
          if ( !real_basis )
            istart *= 2;
          ctxt_.drecv(size,1,&wftmpr[istart],size,i,c_.pc(n));
        }
      }
      else
      {
        if ( iamsending )
        {
          ctxt_.dsend(size,1,&wftmpr[0],size,0,0);
        }
      }
    }

    // process the data
    if ( ctxt_.onpe0() )
    {
      // wftmpr is now complete on task 0
      // wftmpr contains either a real of a complex array

      const string element_type = real_basis ? "double" : "complex";

      if ( encoding == "base64" )
      {
        #if PLT_BIG_ENDIAN
        xcdr.byteswap_double(wftmpr_size,&wftmpr[0]);
        #endif
        int nbytes = wftmpr_size*sizeof(double);
        int outlen = xcdr.nchars(nbytes);
        char* b = new char[outlen];
        assert(b!=0);
        xcdr.encode(nbytes,(byte*) &wftmpr[0],b);
        // Note: optional x0,y0,z0 attributes not used, default is zero
        os << "<grid_function type=\"" << element_type << "\""
           << " nx=\"" << ft.np0()
           << "\" ny=\"" << ft.np1() << "\" nz=\"" << ft.np2() << "\""
           << " encoding=\"base64\">" << endl;
        xcdr.print(outlen,(char*) b, os);
        os << "</grid_function>\n";
        delete [] b;
      }
      else
      {
        // encoding == "text" or unknown encoding
        // Note: optional x0,y0,z0 attributes not used, default is zero
        os << "<grid_function type=\"" << element_type << "\""
           << " nx=\"" << ft.np0()
           << "\" ny=\"" << ft.np1() << "\" nz=\"" << ft.np2() << "\""
           << " encoding=\"text\">" << endl;
        int count = 0;
        for ( int k = 0; k < ft.np2(); k++ )
          for ( int j = 0; j < ft.np1(); j++ )
            for ( int i = 0; i < ft.np0(); i++ )
            {
							#ifndef _CHANGE_FT_CLASS
              int index = ft.index(i,j,k);							
              if ( real_basis )
                os << " " << wftmpr[index];
              else
                os << " " << wftmpr[2*index] << " " << wftmpr[2*index+1];
              if ( count++%4 == 3)
                os << "\n";
							#else
							assert(0); //non-base64 not supported with new FT class 
							#endif
            }
        if ( count%4 != 0 )
          os << "\n";
        os << "</grid_function>\n";
      }
    }
  }
  if ( ctxt_.onpe0() )
    os << "</slater_determinant>" << endl;
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::write(SharedFilePtr& sfp, string encoding, double weight, int ispin,
  int nspin) const
{
#ifndef _CHANGE_FT_CLASS
  FourierTransform ft(*basis_,basis_->np(0),basis_->np(1),basis_->np(2));
#else
  FourierTransform ft(basis_->np());
#endif
  vector<complex<double> > wftmp(ft.np012loc());
  const bool real_basis = basis_->real();
  const int wftmpr_loc_size = real_basis ? ft.np012loc() : 2*ft.np012loc();
  vector<double> wftmpr(wftmpr_loc_size);
  Base64Transcoder xcdr;

  char* wbuf = 0;
  size_t wbufsize = 0;

  // Segment n on process iprow is sent to row (n*nprow+iprow)/(nprow)
  const Context& colctxt = basis_->context();
  const int nprow = ctxt_.nprow();
  vector<int> scounts(nprow), sdispl(nprow), rcounts(nprow), rdispl(nprow);

  string header;
  if ( ctxt_.onpe0() )
  {
    ostringstream ostr_hdr;
    string spin = (ispin > 0) ? "down" : "up";
    ostr_hdr << "<slater_determinant";
    if ( nspin == 2 )
      ostr_hdr << " spin=\"" << spin << "\"";
    ostr_hdr << " kpoint=\"" << basis_->kpoint() << "\"\n"
       << "  weight=\"" << setprecision(12) << weight << "\""
       << " size=\"" << nst() << "\">" << endl;

    ostr_hdr << "<density_matrix form=\"diagonal\" size=\"" << nst() << "\">"
       << endl;
    ostr_hdr.setf(ios::fixed,ios::floatfield);
    ostr_hdr.setf(ios::right,ios::adjustfield);
    for ( int i = 0; i < nst(); i++ )
    {
      ostr_hdr << " " << setprecision(8) << occ_[i];
      if ( i%10 == 9 )
        ostr_hdr << endl;
    }
    if ( nst()%10 != 0 )
      ostr_hdr << endl;
    ostr_hdr << "</density_matrix>" << endl;
    header = ostr_hdr.str();
  }

  // serialize all local columns of c and store in segments seg[n]

  string seg;
  for ( int n = 0; n < nstloc(); n++ )
  {
    seg.clear();
    if ( n == 0 && ctxt_.myrow() == 0 )
      seg = header;

    ostringstream ostr;
    //cout << " state " << n << " is stored on column "
    //     << ctxt_.mycol() << " local index: " << c_.y(n) << endl;

#ifndef _REGULARGDISTRIB
      ft.backward(c_.cvalptr(c_.mloc()*n),&wftmp[0]); /* must use original Fourier transform class */
#else
      // use FFTW3 transforms, by-passing qbox rods design
      typedef std::complex<double> complex;
      using boost::const_multi_array_ref;
      using boost::multi_array_ref;
      using boost::extents;
      using fftw3::plan;
      int np0=basis_->np(0); int np1=basis_->np(1); int np2=basis_->np(2); /*int np012 = np0*np1*np2;*/

#ifndef _MPBC 
      // 3D FFT in conventional PBC

      // use FFTW3 transforms, by-passing qbox rods design
      const_multi_array_ref<complex, 3> wfXYZ (c_.cvalptr(c_.mloc()*n), extents[np0][np1][np2]);
      multi_array_ref<complex, 3>       wfxyz (&wftmp[0],                 extents[np0][np1][np2]);

      plan wf_backward(wfXYZ, wfxyz, FFTW_BACKWARD, FFTW_ESTIMATE);
      wf_backward();
#else
      // FT from reciprocal space to intermediate space to real space in MPBC
      multi_array<complex, 2>           wfXzh (                           extents[np0][np1* np2]);
      multi_array_ref<complex, 3>       wfXYz (wfXzh.origin(),            extents[np0][np1][np2]);

      const_multi_array_ref<complex, 3> wfXYZ (c_.cvalptr(c_.mloc()*n),extents[np0][np1][np2]);
      const_multi_array_ref<complex, 2> wfXZh (c_.cvalptr(c_.mloc()*n),extents[np0][np1* np2]);
      multi_array_ref<complex, 3>       wfxyz (&wftmp[0],                 extents[np0][np1][np2]);

      /* backward 1D FFT along zhat direction: wfXZh -> wfXzh = wfXYz (tmp) */
      /* 2D array, transform last one dimension (-1) */
      //plan wf_z_backward(-1, wfXZh, wfXzh, FFTW_BACKWARD, FFTW_ESTIMATE);  wf_z_backward();
      /* normalization factor (1) */
      //for ( int i = 0; i < np012; i++ ) wfXzh.origin()[i] /= sqrt(np1*np2);

      /* backward 2D FFT along x and y directions: wfXYz (tmp) -> wfxyz (tmp_r) */
      /* 3D array, transform first two dimensions (+2) */
      //plan wf_xy_backward(+2, wfXYz, wfxyz, FFTW_BACKWARD, FFTW_ESTIMATE);  wf_xy_backward();
      /* normalization factor (2) */
      //for ( int i = 0; i < np012; i++ ) wftmp[i] /= sqrt(np0*np1);

      FourierTransform::z_backward_mpbc(wfXZh,wfXzh);
      FourierTransform::xy_backward_mpbc(wfXYz,wfxyz);

#endif //_MPBC
#endif //_REGULARGDISTRIB


    if ( real_basis )
    {
#ifdef _NOANTIALIASING
      assert(0); //do not allow real_ basis if _NOANTIALIASING
#endif
      double *a = (double*) &wftmp[0];
      for ( unsigned i = 0; i < ft.np012loc(); i++ )
        wftmpr[i] = a[2*i];
    }
    else
    {
#ifndef _REGULARGDISTRIB
      memcpy((void*)&wftmpr[0],(void*)&wftmp[0],
             ft.np012loc()*sizeof(complex<double>));
#else
      memcpy((void*)&wftmpr[0],(void*)&wftmp[0],
             ft.np012loc()*sizeof(complex));
#endif
    }

    // find index of last process holding some data
    int lastproc = ctxt_.nprow()-1;
		#ifndef _CHANGE_FT_CLASS
    while ( lastproc >= 0 && ft.np2_loc(lastproc) == 0 ) lastproc--;
		#endif
		
    assert(lastproc>=0);

    // Adjust number of values on each task to have a number of values
    // divisible by three. This is necessary in order to have base64
    // encoding without trailing '=' characters.
    // The last node in the process column may have a number of values
    // not divisible by 3.

    // data now resides in wftmpr, distributed on ctxt_.mycol()
    // All nodes in the process column except the last have the
    // same wftmpr_loc_size
    // Use group-of-three redistribution algorithm to make all sizes
    // multiples of 3. In the group-of-three algorithm, nodes are divided
    // into groups of three nodes. In each group, the left and right members
    // send 1 or 2 values to the center member so that all three members
    // end up with a number of values divisible by three.

    // Determine how many values must be sent to the center-of-three node
    int ndiff;
    const int myrow = ctxt_.myrow();
    if ( myrow == 0 )
    {
      ndiff = wftmpr_loc_size % 3;
      ctxt_.ibcast_send('c',1,1,&ndiff,1);
    }
    else
    {
      ctxt_.ibcast_recv('c',1,1,&ndiff,1,0,ctxt_.mycol());
    }
    // assume that all nodes have at least ndiff values
    if ( myrow <= lastproc ) assert(wftmpr_loc_size >= ndiff);

    // Compute number of values to be sent to neighbors
    int nsend_left=0, nsend_right=0, nrecv_left=0, nrecv_right=0;
    if ( myrow % 3 == 0 )
    {
      // mype is the left member of a group of three
      // send ndiff values to the right if not on the last node
      if ( myrow < lastproc )
        nsend_right = ndiff;
    }
    else if ( myrow % 3 == 1 )
    {
      // mype is the center member of a group of three
      if ( myrow <= lastproc )
        nrecv_left = ndiff;
      if ( myrow <= lastproc-1 )
        nrecv_right = ndiff;
    }
    else if ( myrow % 3 == 2 )
    {
      // mype is the right member of a group of three
      // send ndiff values to the left if not on the first or last node
      if ( myrow <= lastproc && myrow > 0 )
        nsend_left = ndiff;
    }

    double rbuf_left[2], rbuf_right[2], sbuf_left[2], sbuf_right[2];
    int tmpr_size = wftmpr_loc_size;
    if ( nsend_left > 0 )
    {
      for ( int i = 0; i < ndiff; i++ )
        sbuf_left[i] = wftmpr[i];
      ctxt_.dsend(ndiff,1,sbuf_left,ndiff,ctxt_.myrow()-1,ctxt_.mycol());
      tmpr_size -= ndiff;
    }
    if ( nsend_right > 0 )
    {
      for ( int i = 0; i < ndiff; i++)
        sbuf_right[i] = wftmpr[wftmpr_loc_size-ndiff+i];
      ctxt_.dsend(ndiff,1,sbuf_right,ndiff,ctxt_.myrow()+1,ctxt_.mycol());
      tmpr_size -= ndiff;
    }
    if ( nrecv_left > 0 )
    {
      ctxt_.drecv(ndiff,1,rbuf_left,ndiff,ctxt_.myrow()-1,ctxt_.mycol());
      tmpr_size += ndiff;
    }
    if ( nrecv_right > 0 )
    {
      ctxt_.drecv(ndiff,1,rbuf_right,ndiff,ctxt_.myrow()+1,ctxt_.mycol());
      tmpr_size += ndiff;
    }

    // check that size is a multiple of 3 (except on last node)
    // cout << ctxt_.mype() << ": tmpr_size: " << tmpr_size << endl;
    if ( ctxt_.myrow() != lastproc )
      assert(tmpr_size%3 == 0);
    vector<double> tmpr(tmpr_size);

    // Note: all nodes either receive data or send data, not both
    if ( nrecv_left > 0 || nrecv_right > 0 )
    {
      // this node is a receiver
      int index = 0;
      if ( nrecv_left > 0 )
      {
        for ( int i = 0; i < ndiff; i++ )
          tmpr[index++] = rbuf_left[i];
      }
      for ( int i = 0; i < wftmpr_loc_size; i++ )
        tmpr[index++] = wftmpr[i];
      if ( nrecv_right > 0 )
      {
        for ( int i = 0; i < ndiff; i++ )
          tmpr[index++] = rbuf_right[i];
      }
      assert(index==tmpr_size);
    }
    else if ( nsend_left > 0 || nsend_right > 0 )
    {
      // this node is a sender
      int index = 0;
      int istart=0, iend=wftmpr_loc_size;
      if ( nsend_left > 0 )
        istart = ndiff;
      if ( nsend_right > 0 )
        iend = wftmpr_loc_size - ndiff;
      for ( int i = istart; i < iend; i++ )
          tmpr[index++] = wftmpr[i];
      assert(index==tmpr_size);
    }
    else
    {
      // no send and no recv
      for ( int i = 0; i < wftmpr_loc_size; i++ )
        tmpr[i] = wftmpr[i];
      assert(tmpr_size==wftmpr_loc_size);
    }

    // All nodes (except the last) now have a number of values
    // divisible by 3 in tmpr[]

    if ( ctxt_.myrow()!=lastproc ) assert(tmpr_size%3==0);

    // convert local data to base64 and write to outfile

    // tmpr contains either a real or a complex array

    const string element_type = real_basis ? "double" : "complex";

    if ( encoding == "base64" )
    {
      #if PLT_BIG_ENDIAN
      xcdr.byteswap_double(tmpr_size,&tmpr[0]);
      #endif
      int nbytes = tmpr_size*sizeof(double);
      int outlen = xcdr.nchars(nbytes);
      char* b = new char[outlen];
      assert(b!=0);
      xcdr.encode(nbytes,(byte*) &tmpr[0],b);
      // Note: optional x0,y0,z0 attributes not used, default is zero
      if ( ctxt_.myrow() == 0 )
      {
        // if on first row, write grid function header
        ostr << "<grid_function type=\"" << element_type << "\""
           << " nx=\"" << ft.np0()
           << "\" ny=\"" << ft.np1() << "\" nz=\"" << ft.np2() << "\""
           << " encoding=\"base64\">" << endl;
      }
      xcdr.print(outlen,(char*) b, ostr);
      if ( ctxt_.myrow() == lastproc )
        ostr << "</grid_function>\n";
      delete [] b;
    }
    else
    {
      // encoding == "text" or unknown encoding
      // Note: optional x0,y0,z0 attributes not used, default is zero
      if ( ctxt_.myrow() == 0 )
      {
        // if on first row, write grid function header
        ostr << "<grid_function type=\"" << element_type << "\""
             << " nx=\"" << ft.np0()
             << "\" ny=\"" << ft.np1() << "\" nz=\"" << ft.np2() << "\""
             << " encoding=\"text\">" << endl;
      }
      int count = 0;
      for ( int k = 0; k < ft.np2(); k++ )
        for ( int j = 0; j < ft.np1(); j++ )
          for ( int i = 0; i < ft.np0(); i++ )
          {
						#ifndef _CHANGE_FT_CLASS
            int index = ft.index(i,j,k);
            if ( real_basis )
              ostr << " " << tmpr[index];
            else
              ostr << " " << tmpr[2*index] << " " << tmpr[2*index+1];
            if ( count++%4 == 3)
              ostr << "\n";
						#else
						assert(0);
						#endif
          }
      if ( count%4 != 0 )
        ostr << "\n";
      if ( ctxt_.myrow() == lastproc )
        ostr << "</grid_function>\n";
    }
    // copy contents of ostr stringstream to segment
    seg += ostr.str();
    // cout << ctxt_.mype() << ": segment " << n << " size: " << seg.size()
    //      << endl;

    // seg is defined
#if USE_MPI
    // redistribute segments to tasks within each process column

    for ( int i = 0; i < nprow; i++ )
    {
      scounts[i] = 0;
      sdispl[i] = 0;
      rcounts[i] = 0;
      rdispl[i] = 0;
    }

    int idest = (n*nprow+ctxt_.myrow())/nstloc();
    scounts[idest] = seg.size();

    // send sendcounts to all procs
    MPI_Alltoall(&scounts[0],1,MPI_INT,&rcounts[0],1,MPI_INT,colctxt.comm());

    // dimension receive buffer
    int rbufsize = rcounts[0];
    rdispl[0] = 0;
    for ( int i = 1; i < ctxt_.nprow(); i++ )
    {
      rbufsize += rcounts[i];
      rdispl[i] = rdispl[i-1] + rcounts[i-1];
    }
    char* rbuf = new char[rbufsize];

    int err = MPI_Alltoallv((void*)seg.data(),&scounts[0],&sdispl[0],
              MPI_CHAR,rbuf,&rcounts[0],&rdispl[0],MPI_CHAR,colctxt.comm());

    if ( err != 0 )
       cout << ctxt_.mype()
            << " SlaterDet::write: error in MPI_Alltoallv" << endl;

    if ( rbufsize > 0 )
    {
      // append rbuf to wbuf
      char* tmp = 0;
      try
      {
        tmp = new char[wbufsize+rbufsize];
      }
      catch ( bad_alloc )
      {
        cout << ctxt_.mype() << " bad_alloc in wbuf append "
             << " n=" << n
             << " rbufsize=" << rbufsize
             << " wbufsize=" << wbufsize << endl;
      }
      memcpy(tmp,wbuf,wbufsize);
      memcpy(tmp+wbufsize,rbuf,rbufsize);
      delete [] wbuf;
      wbuf = tmp;
      wbufsize += rbufsize;
    }
    delete [] rbuf;
#else
    sfp.file().write(seg.data(),seg.size());
#endif
  }
#if USE_MPI
  // wbuf now contains the data to be written in the correct order

  ctxt_.barrier();

  // compute offsets
  sfp.sync();

  MPI_Offset off;
  long long int local_offset,current_offset;
  current_offset = sfp.offset();

  // compute local offset of next write
  long long int local_size = wbufsize;
  MPI_Scan(&local_size, &local_offset, 1,
           MPI_LONG_LONG, MPI_SUM, ctxt_.comm());
  // add base and correct for inclusive scan by subtracting local_size
  local_offset += current_offset - local_size;
  off = local_offset;

  MPI_Status status;

  // write wbuf from all tasks using computed offset
  int len = wbufsize;
  int err = MPI_File_write_at_all(sfp.file(),off,(void*)wbuf,len,
                               MPI_CHAR,&status);
  if ( err != 0 )
    cout << ctxt_.mype()
         << " error in MPI_File_write_at_all" << endl;
  sfp.set_offset(local_offset+len);

  sfp.sync();

  delete [] wbuf;

  if ( ctxt_.onpe0() )
  {
    string s("</slater_determinant>\n");
    int err = MPI_File_write_at(sfp.file(),sfp.mpi_offset(),(void*) s.data(),
              s.size(),MPI_CHAR,&status);
    if ( err != 0 )
      cout << ctxt_.mype()
           << " error in MPI_File_write, slater_determinant trailer"
           << endl;
    sfp.advance(s.size());
  }
#else
  sfp.file() << "</slater_determinant>\n";
#endif // USE_MPI
}

////////////////////////////////////////////////////////////////////////////////
void SlaterDet::info(ostream& os) const
{
  if ( ctxt_.onpe0() )
  {
    os << "<slater_determinant kpoint=\"" << basis_->kpoint() << "\""
       << " size=\"" << nst() << "\">" << endl;
    os << " sdcontext: " << ctxt_.nprow() << "x" << ctxt_.npcol() << endl;
    //os << " sdcontext: " << ctxt_ << endl;
    os << " basis size: " << basis_->size() << endl;
    os << " c dimensions: "
       << c_.m() << "x" << c_.n()
       << "   (" << c_.mb() << "x" << c_.nb() << " blocks)" << endl;
    os << " <density_matrix form=\"diagonal\" size=\"" << nst() << "\">"
       << endl;
    os << " </density_matrix>" << endl;
    os << "</slater_determinant>" << endl;
  }
}

////////////////////////////////////////////////////////////////////////////////
ostream& operator<<(ostream& os, SlaterDet& sd)
{
  sd.print(os,"text",0.0,0,1);
  return os;
}
// ~/usr/bin-hera/astyle --brackets=attach --indent=tab --indent-col1-comments --pad-oper --delete-empty-lines --add-brackets --keep-one-line-statements --convert-tabs --align-pointer=type --unpad-paren  qbox.hpp
// Editor modelines  -  http://www.wireshark.org/tools/modelines.html
// Local variables:
// c-basic-offset: 4
// tab-width: 4
// indent-tabs-mode: t
// truncate-lines: 1
// End:
/* vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent: */

