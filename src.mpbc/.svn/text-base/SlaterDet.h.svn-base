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
// SlaterDet.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: SlaterDet.h,v 1.27 2008/09/08 15:56:19 fgygi Exp $

#ifndef SLATERDET_H
#define SLATERDET_H

class FourierTransform;
#include "Context.h"
#include "Basis.h"
#include "Matrix.h"

#include "D3vector.h"
#include <vector>
#include <iostream>
#include "Timer.h"
#include <string>
#include <map>

class SharedFilePtr;

typedef std::map<std::string,Timer> TimerMap;

class SlaterDet
{
  private:

  const Context& ctxt_;
  Context* my_col_ctxt_;
  Basis* basis_;
  ComplexMatrix c_;
#ifdef QBOX_HPP_
 public:
#endif
  std::vector<double> occ_;
  std::vector<double> eig_;

  void byteswap_double(size_t n, double* x);
  double fermi(double e, double mu, double fermitemp);

  public:

  mutable TimerMap tmap;

  SlaterDet(const Context& ctxt, D3vector kpoint);
  SlaterDet(const SlaterDet& rhs);
  ~SlaterDet();
  const Context& context(void) const { return ctxt_; }
  const Basis& basis(void) const { return *basis_; }
  const D3vector kpoint(void) const { return basis_->kpoint(); }
  const ComplexMatrix& c(void) const { return c_; }
  ComplexMatrix& c(void) { return c_; }
  const std::vector<double>& occ(void) const { return occ_; }
  const std::vector<double>& eig(void) const { return eig_; }
  int nst(void) const { return c_.n(); }
  int nstloc(void) const { return c_.nloc(); }
  void resize(const UnitCell& cell, const UnitCell& refcell,
              double ecut, unsigned nst);
  void compute_density(FourierTransform& ft, double weight, double* rho) const;

#ifndef _MPBC
  void rs_mul_add(FourierTransform& ft, const double* v, SlaterDet& sdp) const;
#else
  void rs_mul_add_mpbc(FourierTransform& ft, const double* v, SlaterDet& sdp) const;
  void anti_aliasing_mpbc(FourierTransform& ft);
  void compute_J_mpbc(FourierTransform& ft, int state_idx) const;
#endif

  void randomize(double amplitude);
  void cleanup(void);
  void reset(void);
  void init(void);
  void gram(void);
  void riccati(const SlaterDet& sd);
  void lowdin(void);
  void align(const SlaterDet& sd);
  void ortho_align(const SlaterDet& sd);
  std::complex<double> dot(const SlaterDet& sd) const;
  double total_charge(void) const;
  void update_occ(unsigned nel, unsigned nspin);
  void update_occ(unsigned nspin, double mu, double temp);
  double eig(int i) const { return eig_[i]; };
  const double* eig_ptr(void) const { return &eig_[0]; }
  const double* eig_ptr(int i) const { return &eig_[i]; }
  double occ(int i) const { return occ_[i]; };
  const double* occ_ptr(void) const { return &occ_[0]; }
  const double* occ_ptr(int i) const { return &occ_[i]; }
  void set_occ(std::vector<double>& occ)
    { assert(occ_.size()==occ.size()); occ_ = occ; }
  void set_eig(std::vector<double>& eig)
    { assert(eig_.size()==eig.size()); eig_ = eig; }
  void set_eig(std::valarray<double>& eig)
    { assert(eig_.size()==eig.size());
      for ( unsigned i = 0; i < eig.size(); i++ )
        eig_[i] = eig[i];
    }
  double entropy(int nspin) const;
  double ortho_error(void) const;
  double memsize(void) const;
  double localmemsize(void) const;
  SlaterDet& operator=(SlaterDet& rhs);
  void print(std::ostream& os, std::string encoding, double weight, int ispin,
    int nspin) const;
  void write(SharedFilePtr& fh, std::string encoding, double weight, int ispin,
    int nspin) const;
  void info(std::ostream& os) const;
};
std::ostream& operator << ( std::ostream& os, SlaterDet& sd );

class SlaterDetException
{
  public:
  std::string msg;
  SlaterDetException(std::string s) : msg(s) {}
};

#if 0
//#ifdef _REGULARGDISTRIB
#include<boost/multi_array.hpp>
#include"alf/hdf5.hpp"
/* write wave function matrix into matlab file, wf(Zind,Yind,Xind) */
// alf: duplication avoided by use of const_multi_array_ref const&
void Write_WF_to_Matlab(boost::const_multi_array_ref<std::complex<double>, 3> const& wf, char const* name, int n)
{
        FILE *fp; char filename[100]; int np0, np1, np2;
        assert( wf.num_dimensions() == 3 );
        np0 = wf.shape()[0];
        np1 = wf.shape()[1];
        np2 = wf.shape()[2];
        sprintf(filename,"load%s_%d.m",name,n);
        fp = fopen(filename,"w");
        fprintf(fp,"%s_%d_data = load('%s_%d_data.out');\n",name,n,name,n);
        fprintf(fp,"I = sqrt(-1);\n"); 
        fprintf(fp,"%s_%d = permute(reshape(%s_%d_data(:,1)+%s_%d_data(:,2)*I, %d, %d, %d),[3,2,1]);\n",
                   name,n,name,n,name,n,np2,np1,np0);
        fprintf(fp,"clear %s_%d_data\n",name,n);
        fclose(fp);

        sprintf(filename,"%s_%d_data.out",name,n);
        fp = fopen(filename,"w");
        fprintf(fp,"%% real and imag part of (%d x %d x %d) complex array \n",np0,np1,np2);
        fprintf(fp,"%% I = sqrt(-1);\n"); 
        fprintf(fp,"%% %s_%d = permute(reshape(%s_%d_data(:,1)+%s_%d_data(:,2)*I, %d, %d, %d),[3,2,1]);\n",
                   name,n,name,n,name,n,np0,np1,np2);
        for(unsigned i=0; i<wf.shape()[0]; i++)
          for(unsigned j=0; j<wf.shape()[1]; j++)
            for(unsigned k=0; k<wf.shape()[2]; k++)
            {
              fprintf(fp,"%20.12e     %20.12e\n",real(wf[i][j][k]),imag(wf[i][j][k]));
            }
        fclose(fp);
}

#ifdef _USE_HDF5
void Write_WF_to_HDF5(Basis const * basis, boost::const_multi_array_ref<std::complex<double>, 3> const& wf, char const* name, int n)
{
        boost::multi_array<double, 2> cell(boost::extents[3][3]);
        cell[0]=boost::const_multi_array_ref<double, 1>((double*)&basis->cell().a(0), boost::extents[3]);
        cell[1]=boost::const_multi_array_ref<double, 1>((double*)&basis->cell().a(1), boost::extents[3]);
        cell[2]=boost::const_multi_array_ref<double, 1>((double*)&basis->cell().a(2), boost::extents[3]);
        hdf5::file f(std::string(name)+"_"+boost::lexical_cast<std::string>(n)+".hdf5");
        f["cell"]<<cell;
				f["array"]<<wf; 
//        hdf5::dataset(f, "array", hdf5::type::native< double >::id(), 
//          copy_append(boost::shape_hsize_t(wf),2) 
//        ).write( (double const*)wf.origin() );
}
#endif
//#endif _REGULARGDISTRIB
#endif

#endif


