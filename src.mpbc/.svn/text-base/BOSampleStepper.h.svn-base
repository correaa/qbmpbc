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
// BOSampleStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: BOSampleStepper.h,v 1.8 2008/09/08 15:56:18 fgygi Exp $

#ifndef BOSAMPLESTEPPER_H
#define BOSAMPLESTEPPER_H

#include "SampleStepper.h"
#include "EnergyFunctional.h"
#include "Sample.h"
#include "Wavefunction.h"

class WavefunctionStepper;
class IonicStepper;

class BOSampleStepper : public SampleStepper
{
  private:

  Wavefunction dwf;
  Wavefunction* wfv;
  int nitscf_;
  int nite_;
  ChargeDensity cd_;			//alf: Steppers own the charge density!
  EnergyFunctional ef_;
  std::ostream& cout; //alf: sorry about this name, it is so, so I don't have to change all the cout<<" <eign " xml-output

  WavefunctionStepper* wf_stepper;
 public:
  IonicStepper* ionic_stepper;
 private:
  // Do not allow construction of BOSampleStepper unrelated to a Sample
  BOSampleStepper(void);

  public:

  mutable TimerMap tmap;

  void step(int niter);

  BOSampleStepper(Sample& s, int nitscf, int nite, std::ostream& os=std::cout);
  ~BOSampleStepper();
  EnergyFunctional const& ef() const{return ef_;}
};
#endif
