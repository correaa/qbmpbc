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
// AndersonMixer.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: AndersonMixer.C,v 1.7 2008/09/08 15:56:17 fgygi Exp $

#include "AndersonMixer.h"
#include "blas.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
void AndersonMixer::restart(void)
{
  extrapolate_ = false;
}

////////////////////////////////////////////////////////////////////////////////
void AndersonMixer::update(const double* f, double* theta, double* fbar)
{
    *theta = 0.0;
    if ( extrapolate_ )
    {
      int ione=1;
      valarray<double> tmp0(n_);

      // compute theta = - a / b
      // tmp0 = delta_F = f - flast
      for ( int i = 0; i < n_; i++ )
        tmp0[i] = f[i] - flast_[i];

      // a = delta_F * F
      double a = ddot(&n_, &tmp0[0], &ione, f, &ione);

      // b = delta_F * delta_F
      double b = ddot(&n_, &tmp0[0], &ione, &tmp0[0], &ione);

      if ( pctxt_ != 0 )
      {
        double work[2] = { a, b };
        pctxt_->dsum(2,1,work,2);
        a = work[0];
        b = work[1];
      }

      if ( b != 0.0 )
        *theta = - a / b;
      else
        *theta = 0.0;

      // test if residual has increased
      if ( *theta <= -1.0 )
      {
        *theta = theta_nc_;
      }

      *theta = min(theta_max_,*theta);
    }

    // fbar = f + theta * ( f - flast )
    // flast_ = f_
    for ( int i = 0; i < n_; i++ )
    {
      const double ftmp = f[i];
      fbar[i] = ftmp + *theta * ( ftmp - flast_[i] );
      flast_[i] = ftmp;
    }
    extrapolate_ = true;
}
