#include<boost/mpi.hpp>
#include "../../src.wei/driver/qbox.hpp"

using qbox::vector;
using qbox::cell;


int main(int argc, char* argv[]) {
  boost::mpi::environment env(argc, argv);
  qbox::session qb;
  
  qb.set_nrowmax  (1u);

  qb.set_cell(cell( vector(12,  0,  0), vector( 0, 12,  0), vector( 0,  0, 12) ));

  qb.kpoint_delete(vector(0     , 0 , 0));
  qb.kpoint_add   (vector(0.0001, 0., 0), 1.);
  
  qb.species("carbon"  , "carbon_pbe.xml"  );

  qb.atom("C" , "carbon"  , vector( 0.00000000,  0.00000000,  0.00000000));
  //qb.atom("C" , "carbon"  , vector( 6.00000000,  0.00000000,  0.00000000));

  qb.set_ecut(35);

  qb.set_xc("PBE");

  qb.set_wf_dyn("PSD");

  qb.set_ecutprec(20);

  qb.set_nempty (4);

  /* to be implemented later */
  qb.set_fermi_temp (1000);

  qb.randomize_wf();

  //try{qb.load("wf.xml");}catch(...){}

  /* creating output file for total energy and print header information */	
  std::ofstream ofs("Carbon_etotal.dat");
  ofs<<"# total energy   kinetic energy (ekin), local potential (ehart),  non-local potential (enl)"<<std::endl;

  qbox::energies ens = qb.run(0, 50, 10);
  ofs<<ens.etotal()<<" "<<ens.ekin()<<" "<<ens.ehart()<<" "<<ens.enl()<<std::endl;

  qb.save("wf.xml");

  /* creating output file for total energy and print header information */	
  std::ofstream ofs_eig("ch4_eig.dat");
  ofs_eig<<"# eig 0, 1, ...\n";

  /* write eigenvalues to file */
  SlaterDet const* const mysd = qb.wavefunction().sd(0,0);
  for(int i = 0; i < mysd->nst(); i ++ ){
    ofs_eig<< mysd->eig(i)<<" ";
  }
  ofs_eig<<std::endl;

  /* compute density information from qbox (right now only works in PBC) */
  field<double> density=qbox::charge(qb.wavefunction());

  /* add name-value-pair (nvp) to hdf5 */
  hdf5::file("Carbon_density.h5")<<hdf5::nvp<field<double> >("density", density);

  return 0;
}

