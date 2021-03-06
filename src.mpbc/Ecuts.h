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
// Ecuts.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Ecuts.h,v 1.4 2008/09/08 15:56:18 fgygi Exp $

#ifndef ECUTS_H
#define ECUTS_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include "Sample.h"

class Ecuts : public Var
{
  Sample *s;

  public:

  char const*name ( void ) const { return "ecuts"; };

  int set ( int argc, char **argv )
  {
    if ( argc != 2 )
    {
      if ( ui->onpe0() )
      cout << " ecuts takes only one value" << endl;
      return 1;
    }

    double v = atof(argv[1]);
    if ( v < 0.0 )
    {
      if ( ui->onpe0() )
        cout << " ecuts must be non-negative" << endl;
      return 1;
    }

    s->ctrl.ecuts = 0.5 * v;

    return 0;
  }

  string print (void) const
  {
     ostringstream st;
     st.setf(ios::left,ios::adjustfield);
     st << setw(10) << name() << " = ";
     st.setf(ios::right,ios::adjustfield);
     st << setw(10) << 2 * s->ctrl.ecuts;
     return st.str();
  }

  Ecuts(Sample *sample) : s(sample) { s->ctrl.ecuts = 0.0; }
};
#endif
