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
// Control.h:
//
////////////////////////////////////////////////////////////////////////////////
// $Id: Control.h,v 1.14 2008/09/08 15:56:18 fgygi Exp $

#ifndef CONTROL_H
#define CONTROL_H

#include <string>
#include <vector>

struct Control
{
  // control variables
  std::string debug, timing;
  std::string wf_dyn, atoms_dyn; // dynamics string flags
  int nite;
  double emass;       // electron mass

  double fermi_temp;  // temperature of Fermi distribution
  double ecutprec;

  std::string wf_diag;

  std::string tcp;
  double tcp_rcut;
  double tcp_sigma;

  double gms_mix; // mixing factor for generalized minimum spread functions

  std::string thermostat;
  double th_temp,th_time, th_width; // thermostat control

  std::string stress;
  std::string cell_dyn;
  std::string cell_lock;
  double cell_mass;
  double ecuts;         // confinement potential energy cutoff
  double ext_stress[6]; // external stress tensor: xx,yy,zz,xy,yz,xz

  std::string xc;
  std::string spin;
  int delta_spin;

  double dt;
  int iprint;
  int timeout;

  double charge_mix_coeff;
  double charge_mix_rcut;
};
#endif