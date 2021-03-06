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
// SpeciesCmd.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SpeciesCmd.C,v 1.12 2008/09/08 15:56:19 fgygi Exp $

#include "SpeciesCmd.h"
#include "SpeciesReader.h"
#include "Species.h"
using namespace std;

class Species;

////////////////////////////////////////////////////////////////////////////////
int SpeciesCmd::action(int argc, char **argv)
{
  if ( argc != 3 )
  {
    if ( ui->onpe0() )
      cout << "  Use: species name uri" << endl;
    return 1;
  }

  if ( ui->onpe0() )
    cout << "  SpeciesCmd: defining species " << argv[1]
         << " as " << argv[2] << endl;

  SpeciesReader sp_reader(s->ctxt_);

  Species* sp = new Species(s->ctxt_,argv[1]);

  try
  {
    sp_reader.readSpecies(*sp,argv[2]);
    sp_reader.bcastSpecies(*sp);
    s->atoms.addSpecies(sp,argv[1]);
  }
  catch ( const SpeciesReaderException& e )
  {
    cout << " SpeciesReaderException caught in SpeciesCmd" << endl;
    cout << " SpeciesReaderException: cannot define Species" << endl;
  }
  catch (...)
  {
    cout << " SpeciesCmd: cannot define Species" << endl;
  }

  return 0;
}
