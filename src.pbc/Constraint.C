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
//  Constraint.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Constraint.C,v 1.4 2008/09/08 15:56:18 fgygi Exp $

#include "Constraint.h"
#include <iostream>
using namespace std;

ostream& operator << ( ostream &os, Constraint &c )
{
  return c.print(os);
}

