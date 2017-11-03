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
// ComputeMLWFCmd.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ComputeMLWFCmd.h,v 1.4 2008/09/08 15:56:18 fgygi Exp $

#ifndef COMPUTEMLWFCMD_H
#define COMPUTEMLWFCMD_H

#include <iostream>
#include <stdlib.h>
#include <string>

#include "UserInterface.h"
#include "Sample.h"
#include "MLWFTransform.h"

class ComputeMLWFCmd : public Cmd
{
  private:

  MLWFTransform* mlwft;

  public:
  Sample *s;

  ComputeMLWFCmd(Sample *sample) : s(sample) {};

  char const*name(void) const { return "compute_mlwf"; }

  char const*help_msg(void) const
  {
    return
    "\n compute_mlwf\n\n"
    " syntax: compute_mlwf\n\n"
    "   The compute_mlwf command computes maximally localized \n"
    " Wannier functions.\n\n";
  }

  int action(int argc, char **argv);

  ComputeMLWFCmd();
  ~ComputeMLWFCmd(){}
};
#endif
