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
// ChargeDensity.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ChargeDensity.C,v 1.19 2008/09/08 15:56:18 fgygi Exp $

#include "ChargeDensity.h"
#include "Basis.h"
#include "Wavefunction.h"
#include "FourierTransform.h"
#include "SlaterDet.h"

#include <iomanip>
using namespace std;

////////////////////////////////////////////////////////////////////////////////
ChargeDensity::ChargeDensity(const Wavefunction& wf) : ctxt_(wf.context()),
wf_(wf), vcontext_(wf.sd(0,0)->basis().context())

{
  vbasis_ = new Basis(vcontext_, D3vector(0,0,0));

#ifdef _NOANTIALIASING
  vbasis_->resize(wf.cell(),wf.refcell(),1.0*wf.ecut());
  const Basis& vb = *vbasis_;

  int np0v = vbasis_->np(0);
  int np1v = vbasis_->np(1);
  int np2v = vbasis_->np(2);
#else
  vbasis_->resize(wf.cell(),wf.refcell(),4.0*wf.ecut());
  const Basis& vb = *vbasis_;

  // define vft_, FT on vbasis context for transforming the density

  // add 2 to grid size to avoid aliasing when using non-zero k-points
  // adding 1 would suffice, but add 2 to keep even numbers
  cout << " vbasis_->np012 " << vbasis_->np(0) << "/" << vbasis_->np(1)
       << "/" << vbasis_->np(2) << endl;

  int np0v = vbasis_->np(0)+2;
  int np1v = vbasis_->np(1)+2;
  int np2v = vbasis_->np(2)+2;
  while (!vbasis_->factorizable(np0v)) np0v += 2;
  while (!vbasis_->factorizable(np1v)) np1v += 2;
  while (!vbasis_->factorizable(np2v)) np2v += 2;
#endif

#if 1
  cout << " ChargeDensity: vbasis: " << endl;
  cout << " idxmin: " << vb.idxmin(0) << "/" << vb.idxmin(1)
       << "/" << vb.idxmin(2) << endl;
  cout << " idxmax: " << vb.idxmax(0) << "/" << vb.idxmax(1)
       << "/" << vb.idxmax(2) << endl;
  cout << " vft grid: " << np0v << "/" << np1v << "/" << np2v << endl;
#endif
  vft_ = new FourierTransform(*vbasis_,np0v,np1v,np2v);

  rhor.resize(wf.nspin());
  rhog.resize(wf.nspin());
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    rhor[ispin].resize(vft_->np012loc());
    rhog[ispin].resize(vbasis_->localsize());
  }
  rhotmp.resize(vft_->np012loc());

  // FT for interpolation of wavefunctions on the fine grid
  ft_.resize(wf.nkp());
  for ( int ispin = 0; ispin < wf.nspin(); ispin++ )
  {
    for ( int ikp = 0; ikp < wf.nkp(); ikp++ )
    {
      const Basis& wb = wf.sd(ispin,ikp)->basis();
#if DEBUG
      cout << " ChargeDensity: ikp=" << ikp << " wbasis: " << endl;
      cout << " idxmin: " << wb.idxmin(0) << "/" << wb.idxmin(1)
           << "/" << wb.idxmin(2) << endl;
      cout << " idxmax: " << wb.idxmax(0) << "/" << wb.idxmax(1)
           << "/" << wb.idxmax(2) << endl;
#endif
      // check that no aliasing error can occur

      assert(2*np0v > 2*wb.idxmax(0)+vb.idxmax(0));
      assert(2*np1v > 2*wb.idxmax(1)+vb.idxmax(1));
      assert(2*np2v > 2*wb.idxmax(2)+vb.idxmax(2));
      assert(2*np0v > -2*wb.idxmin(0)-vb.idxmin(0));
      assert(2*np1v > -2*wb.idxmin(1)-vb.idxmin(1));
      assert(2*np2v > -2*wb.idxmin(2)-vb.idxmin(2));
      ft_[ikp] = new FourierTransform(wb,np0v,np1v,np2v);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
ChargeDensity::~ChargeDensity(void)
{
  delete vbasis_;
  delete vft_;
  for ( int ikp = 0; ikp < ft_.size(); ikp++ )
    delete ft_[ikp];
  for ( TimerMap::iterator i = tmap.begin(); i != tmap.end(); i++ )
  {
    double time = (*i).second.real();
    double tmin = time;
    double tmax = time;
    ctxt_.dmin(1,1,&tmin,1);
    ctxt_.dmax(1,1,&tmax,1);
    if ( ctxt_.myproc()==0 )
    {
      cout << "<timing name=\""
           << setw(15) << (*i).first << "\""
           << " min=\"" << setprecision(3) << setw(9) << tmin << "\""
           << " max=\"" << setprecision(3) << setw(9) << tmax << "\"/>"
           << endl;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
void ChargeDensity::update_density(void)
{
  assert(rhor.size() == wf_.nspin());
  const double omega = vbasis_->cell().volume();

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    assert(rhor[ispin].size() == vft_->np012loc() );
    assert(rhotmp.size() == vft_->np012loc() );
    const int rhor_size = rhor[ispin].size();
    tmap["charge_compute"].start();
    double *p = &rhor[ispin][0];
    for ( int i = 0; i < rhor_size; i++ )
      p[i] = 0.0;
    for ( int ikp = 0; ikp < wf_.nkp(); ikp++ )
    {
      assert(rhor[ispin].size()==ft_[ikp]->np012loc());
      wf_.sd(ispin,ikp)->compute_density(*ft_[ikp],
          wf_.weight(ikp), &rhor[ispin][0]);
    }
    tmap["charge_compute"].stop();

    // sum over kpoints: sum along rows of kpcontext
    tmap["charge_rowsum"].start();
    wf_.kpcontext()->dsum('r',vft_->np012loc(),1,
      &rhor[ispin][0],vft_->np012loc());
    tmap["charge_rowsum"].stop();

    // check integral of charge density
    // compute Fourier coefficients of the charge density
    double sum = 0.0;
    const double *const prhor = &rhor[ispin][0];
    tmap["charge_integral"].start();
    for ( int i = 0; i < rhor_size; i++ )
    {
      const double prh = prhor[i];
      sum += prh;
      rhotmp[i] = complex<double>(omega * prh, 0.0);
    }
    sum *= omega / vft_->np012();

    // sum on all indices except spin: sum along columns of spincontext
    wf_.spincontext()->dsum('c',1,1,&sum,1);
    tmap["charge_integral"].stop();
    if ( ctxt_.onpe0() )
    {
      cout.setf(ios::fixed,ios::floatfield);
      cout.setf(ios::right,ios::adjustfield);
      cout << "  total_electronic_charge: " << setprecision(8) << sum
           << endl;
    }

    tmap["charge_vft"].start();
    vft_->forward(&rhotmp[0],&rhog[ispin][0]);
    tmap["charge_vft"].stop();
  }
}
////////////////////////////////////////////////////////////////////////////////
void ChargeDensity::update_rhor(void)
{
  // recalculate rhor from rhog
  assert(rhor.size() == wf_.nspin());
  const double omega = vbasis_->cell().volume();
  assert(omega!=0.0);
  const double omega_inv = 1.0 / omega;

  for ( int ispin = 0; ispin < wf_.nspin(); ispin++ )
  {
    assert(rhor[ispin].size() == vft_->np012loc() );
    assert(rhotmp.size() == vft_->np012loc() );

    vft_->backward(&rhog[ispin][0],&rhotmp[0]);

    const int rhor_size = rhor[ispin].size();
    double *const prhor = &rhor[ispin][0];
    for ( int i = 0; i < rhor_size; i++ )
    {
      prhor[i] = rhotmp[i].real() * omega_inv;
    }
  }
}
