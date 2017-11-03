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
// sinft.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: sinft.C,v 1.4 2008/09/08 15:56:20 fgygi Exp $

#include "sinft.h"
#include <math.h>
#include <assert.h>
#include <vector>
#include <complex>
#include<iostream>
#ifdef USE_FFTW
//#include "fftw.h"
//namespace fftw3{
	#include<fftw3.h>
//}
#endif

using namespace std;
//alf: I had to update this library to fftw3 because fftw2 was giving segmentation fault.
void sinft(int n, double *y)
{
  clog<<"enter sinf"<<endl;
#ifdef USE_FFTW
  vector<complex<double> > zin(2*n), zout(2*n);
  zin[0] = 0.0;
  for ( int i = 1; i < n; i++ )
  {
    const double t = y[i];
    zin[i] = t;
    zin[2*n-i] = -t;
  }
  fftw_plan fwplan=fftw_plan_dft_1d(2*n, (fftw_complex*)&zin[0], (fftw_complex*)&zout[0], /* FFTW_FORWARD*/ -1, FFTW_ESTIMATE);
	fftw_execute(fwplan);
	//fftw_plan fwplan;
	//fwplan = fftw_create_plan(2*n,FFTW_FORWARD,FFTW_ESTIMATE);
  //fftw_one(fwplan,(fftw_complex*)&zin[0],(fftw_complex*)&zout[0]);
  for ( int i = 0; i < n; i++ )
  {
    y[i] = -0.5 * imag(zout[i]);
  }
  fftw_destroy_plan(fwplan);
  //fftw_destroy_plan(fwplan);
#else
  fprintf(stderr,"sinft() not implemented yet!\n");
#endif
  clog<<"exit sinf"<<endl;
}

void cosft1(int n, double *y)
{
#ifdef USE_FFTW
  /* Note: the array y contains n+1 elements */
  vector<complex<double> > zin(2*n), zout(2*n);

  zin[0] = y[0];
  for ( int i = 1; i < n+1; i++ )
  {
    const double t = y[i];
    zin[i] = t;
    zin[2*n-i] = t;
  }
  fftw_plan fwplan=fftw_plan_dft_1d(2*n, (fftw_complex*)&zin[0], (fftw_complex*)&zout[0], /* FFTW_FORWARD*/ -1, FFTW_ESTIMATE);
	fftw_execute(fwplan);
//  fftw_plan fwplan;
//  fwplan = fftw_create_plan(2*n,FFTW_FORWARD,FFTW_ESTIMATE);
 // fftw_one(fwplan,(fftw_complex*)&zin[0],(fftw_complex*)&zout[0]);
  y[0] = 0.5 * real(zout[0]);
  for ( int i = 1; i < n; i++ )
  {
    y[i] = 0.5 * real(zout[i]);
  }
  fftw_destroy_plan(fwplan);
  //fftw_destroy_plan(fwplan);
#else
  fprintf(stderr,"cosft1() not implemented yet!\n");
#endif
}
