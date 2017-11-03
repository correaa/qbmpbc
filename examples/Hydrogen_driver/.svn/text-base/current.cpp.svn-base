#if 0
  echo "compiling..."
  rm -f $0.x
  mpicxx $0 -o $0.x -I../../src.wei -I../../src.wei/include -Wl,-rpath=$HOME/usr/lib -Wall -D_MPBC -DUSE_MPI -DUSE_HDF5 -DSCALAPACK -L../../src.wei -lqb -lfftw -lfftw3  -lboost_system -lboost_filesystem -lboost_mpi -lhdf5 $HOME/usr/lib/libscalapack.a $HOME/usr/lib/libblacsCinit-mpi.a $HOME/usr/lib/libblacs-mpi.a $HOME/usr/lib/libblacsCinit-mpi.a -llapack -lf77blas -latlas -lgfortran -lxerces-c -lmpich -lgsl -lcblas -latlas && echo '--' && echo "done. running ... " && ./$0.x $1 $2 $3 $4 $5 $6 $7 $8 $9
  exit
#endif

#include "driver/qbox.hpp"
#include<boost/numeric/interval/io.hpp> //for cout<<interval
#include<iomanip>
#include<alf/signal.hpp>

using qbox::vector; 
using qbox::cell; 
using hdf5::file; 
using hdf5::nvp; 
using hdf5::make_const_nvp; 
using hdf5::make_nvp;
using namespace boost::filesystem;

int main(int argc, char* argv[]) {

#ifdef _NO_BOOST_MPI
  MPI_Init(&argc,&argv);
#endif

  ofstream qbout("qb_log.xml");
  qbox::session qb(qbout);
  qb.set_nrowmax(1u);

  if(exists("wf.xml")) {

     qb.load("wf.xml");

  } else {

     qb.kpoint_delete(vector(0     , 0 , 0)		);
     qb.kpoint_add   (vector(0.0001, 0., 0), 1.);
     qb.set_nempty   (2);
     qb.set_cell(cell(
                       vector(30,  0,  0),
                       vector( 0, 30,  0),
		       vector( 0,  0, 30)
		));
     std::cout<<"magnetic field "<<qb.magnetic_field()<<std::endl;

     qb.species("hydrogen", "hydrogen_bare_extended.xml");
     vector origin(15,15,15);
     qb.atom("H", "hydrogen", vector( 0., 0, 0)+origin);
     qb.export_xyz("atoms.xyz");
     //qb.set_ecut(40);
     //qb.set_ecutprec(30);
     qb.set_ecut(60);
     qb.set_ecutprec(50);
     qb.set_xc("PBE"); 
     qb.set_wf_dyn("PSD");
     qb.randomize_wf();

  }
  ofstream ofs("etotal.dat");
  ofs<<"# energy "<<endl;
  ofstream ofs_eig("eigv.dat");
  ofs_eig<<"# eig 0, 1, ...\n";

  system_signal::terminal_interrupt::capture();
  clog<<"entering loop ------------------------- press ctrl-c to break"<<endl;

  for ( int i = 0; i < 300 ; i ++ ) {
    qbox::energies ens = qb.run(0, 6, 6);
    ofs<<ens.etotal()<<endl;
    SlaterDet const* const mysd = qb.wavefunction().sd(0,0);
    for(int i = 0; i < mysd->nst(); i ++ ){
       ofs_eig<<std::setw(15) << std::setprecision(15) << mysd->eig(i)<<" ";
    }
    ofs_eig<<endl;
    if(system_signal::terminal_interrupt::flag){
       clog<<"signal handled, result: break the for loop"<<endl;
       break;
    }
  }
  system_signal::terminal_interrupt::release();

  qb.export_xyz("atoms.xyz");
  qb.save("wf.xml");
  file("wf0.h5")<<make_const_nvp("wavefunction", qbox::wavefunction(qb.wavefunction(), 0));
  file("wf1.h5")<<make_const_nvp("wavefunction", qbox::wavefunction(qb.wavefunction(), 1));
  file("density.h5")<<make_const_nvp("density", qbox::density(qb.wavefunction()));
  file("current.h5")<<make_const_nvp("current", qbox::current(qb.wavefunction(),1));
  file("magnetic.h5")<<make_const_nvp("magnetic", qbox::magnetic(qbox::current(qb.wavefunction(),1)));
  std::cout<<"magnetic field "<<qb.magnetic_field()<<std::endl;
  qb.export_commands("commands.qbox");

  return 0;
}

