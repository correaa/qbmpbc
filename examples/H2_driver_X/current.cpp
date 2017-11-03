#if compilation_instructions
mpiCC $0 -o $0.x -Wfatal-errors -Wall \
 -L$HOME/usr/lib$CLUSTER \
 -I../../src.mpbc -I../../src.mpbc/include \
 -D_MPBC -DUSE_MPI -DUSE_HDF5 -DSCALAPACK -D_NO_BOOST_MPI -DUSE_XERCES -DXERCESC_3_0_1 \
 -L../../src.mpbc -lqb -lfftw -lfftw3 -lgfortran -lz \
 -lxerces-c -lhdf5 \
 /g/g91/correaa/usr/libaztec/libscalapack.a /g/g91/correaa/usr/libaztec/libblacsCinit-mpi.a /g/g91/correaa/usr/libaztec/libblacs-mpi.a /g/g91/correaa/usr/libaztec/libblacsCinit-mpi.a -llapack \
 -lboost_system -lboost_filesystem \
 && ./$0.x $1 $2 $3 $4 $5 $6 $7 $8 $9 && rm -f $0.x
exit
#endif

/* usage:
   ./current x/y ecut_min ecut_max ecut_del niter nitscf nite printfreq nempty
*/

#include "driver/qbox.hpp"
#include<boost/numeric/interval/io.hpp> //for cout<<interval
#include<iomanip>
#include<alf/signal.hpp>
#include "boost/units/systems/atomic.hpp"
#include<boost/filesystem/fstream.hpp>
using qbox::cell; 
using hdf5::file; 
using hdf5::nvp; 
using hdf5::make_const_nvp; 
using hdf5::make_nvp;
using std::string;
using boost::lexical_cast;
using namespace boost::filesystem;
int main(int argc, char* argv[]){
#ifdef _NO_BOOST_MPI
  MPI_Init(&argc,&argv);
#else
  boost::mpi::environment env(argc, argv);
#endif

  //qbox::session qb(std::cout);

  double h2_bond = 1.398397338; /*bohr*/
  qbox::vector origin(15,15,15);

  qbox::vector xhat(1,0,0);
  qbox::vector yhat(0,1,0);
  qbox::vector mol_dir;
  string mol_dir_str;
  double ecut_min=10;
  double ecut_max=300;
  double ecut_del=5;
  int niter=50;
  int nitscf=6;
  int nite=6;
  int printfreq=10;
  int nempty=0;

  if (argc>1) {
    if ((string(argv[1])=="x")||(string(argv[1])=="X")) {
       clog<<"mol_dir  = x"<<endl;
       mol_dir_str="x";
       mol_dir = xhat;
    }
    else if ((string(argv[1])=="y")||(string(argv[1])=="Y")) {
       clog<<"mol_dir  = y"<<endl;
       mol_dir_str="y";
       mol_dir = yhat;
    }
    else {
       clog<<"mol_dir not assigned!"<<endl;
       return -1;
    }
    if (argc>2) ecut_min =  boost::lexical_cast<double>(argv[2]);
    if (argc>3) ecut_max =  boost::lexical_cast<double>(argv[3]);
    if (argc>4) ecut_del =  boost::lexical_cast<double>(argv[4]);
    if (argc>5) niter    =  boost::lexical_cast<int>(argv[5]);
    if (argc>6) nitscf   =  boost::lexical_cast<int>(argv[6]);
    if (argc>7) nite     =  boost::lexical_cast<int>(argv[7]);
    if (argc>8) printfreq=  boost::lexical_cast<int>(argv[8]);
    if (argc>9) nempty   =  boost::lexical_cast<int>(argv[9]);
  }
  else {
    clog<<"mol_dir not assigned!"<<endl;
    return -1;
  }

  qbox::vector h1_pos(origin - mol_dir*h2_bond/2);
  qbox::vector h2_pos(origin + mol_dir*h2_bond/2);

  boost::filesystem::ofstream convergence("convergence."+mol_dir_str+".dat");
  convergence<<"# ecut,   Bind_x,   Bind_y,   Bind_z,    Bext,    sigma_x(ppm),   sigma_y,   sigma_z"<<endl;

  for(double ecut=ecut_min; ecut<=ecut_max; ecut+=ecut_del){
      clog<<"starting ecut "<<ecut<<endl;
      boost::filesystem::ofstream qbout("wf/qb_log-"+boost::lexical_cast<std::string>(ecut)+"."+mol_dir_str+ ".xml");
      qbox::session qb(qbout); 
      qb.set_nrowmax(1u); //one

      qb.kpoint_delete(qbox::vector(0     , 0 , 0)	);
      if(!exists("wf/wf.ecut"+lexical_cast<string>(ecut)+"."+mol_dir_str+".xml")) {
         qb.kpoint_add   (qbox::vector(1e-8, 0., 0), 1.);
      }

      qb.set_nempty   (nempty);
      qb.set_cell(cell(
			 qbox::vector(30,  0,  0),
			 qbox::vector( 0, 30,  0),
			 qbox::vector( 0,  0, 30)
      ));
      std::cout<<"magnetic field "<<qb.external_magnetic_field()<<std::endl;
      //qb.species("hydrogen", "hydrogen_bare_extended_withelectron.xml");
      qb.species("hydrogen", "../h2o/hydrogen_pbe.xml");

      qb.atom("H1", "hydrogen", h1_pos); 
      qb.atom("H2", "hydrogen", h2_pos);
      qb.export_xyz("atoms.xyz");

      qb.set_ecut(ecut);
      //qb.set_ecuts(ecut-5);
      qb.set_xc("PBE"); 
      qb.set_wf_dyn("PSD");
      qb.set_ecutprec(ecut/4.0);

      if(exists("wf.ecut"+lexical_cast<string>(ecut)+"."+mol_dir_str+".xml")) {

         qb.load("wf.ecut"+lexical_cast<string>(ecut)+"."+mol_dir_str+".xml");
 
      } else {

         qb.randomize_wf();

      }

      //output monitor files
      boost::filesystem::ofstream ofs_ene("wf/etotal."+boost::lexical_cast<std::string>(ecut)+"."+mol_dir_str+".dat");
      ofs_ene<<"# energy "<<std::setprecision(15)<<endl;
      boost::filesystem::ofstream ofs_eig("wf/eigv."+boost::lexical_cast<std::string>(ecut)+"."+mol_dir_str+".dat");
      ofs_eig<<"# eig 0, 1, ..."<<std::setprecision(15)<<endl;
      boost::filesystem::ofstream ofs_mag("wf/mag."+boost::string_cast(ecut)+"."+mol_dir_str+".dat");
      ofs_mag<<"# magnetic field at point, relative error, factor"<<std::setprecision(15)<<endl;

      clog<<"begining for loop... (press ctrl-c to break)"<<endl;
      system_signal::terminal_interrupt::capture();

      int iter;
      for ( iter = 0; iter <= niter ; iter ++ ){

         qbox::energies ens = qb.run(0, nitscf, nite);
         ofs_ene<<ens.etotal()<<endl;
         SlaterDet const* const mysd = qb.wavefunction().sd(0,0);
         for(int i = 0; i < mysd->nst(); i ++ ) ofs_eig<<mysd->eig(i)<<" "<<endl;

         if(system_signal::terminal_interrupt::flag()){
            clog<<"signal handled, result: break the for loop"<<endl;
            break;
         }

         if ((printfreq == 1) || ( iter % printfreq == 0 && iter > 0 )) {
            clog<<"saving to wf.xml"<<endl;
            qb.save("wf/wf.ecut"+lexical_cast<string>(ecut)+"."+mol_dir_str+".xml");
            clog<<"saving .h5 files"<<endl;
            file("wf/density"+lexical_cast<string>(ecut)+"."+mol_dir_str+".h5")
                       <<make_const_nvp("density", qbox::density(qb.wavefunction()));
            file("wf/current"+lexical_cast<string>(ecut)+"."+mol_dir_str+".h5")
                       <<make_const_nvp("current", qbox::current(qb.wavefunction()));
            file("wf/magnetic"+lexical_cast<string>(ecut)+"."+mol_dir_str+".h5")
                       <<make_const_nvp("magnetic", qbox::magnetic(qbox::current(qb.wavefunction())));
 
            qbox::field<double> den=qbox::density(qb.wavefunction());
            qbox::field<qbox::vector> cur=qbox::current(qb.wavefunction());
            qbox::vector mag_at_h2_dft;
            qbox::field<qbox::vector> mag=qbox::magnetic(cur, h2_pos, &mag_at_h2_dft);
            qbox::vector mag_at_h2 = mag(h2_pos);
            convergence<<ecut<<" "<<mag_at_h2<<" "<<qb.external_magnetic_field()<<" "
                 <<mag_at_h2/qb.external_magnetic_field()/(137.036*137.036)*4.*3.1415926*1e6<<endl;
            convergence<<ecut<<" "<<mag_at_h2_dft<<" "<<qb.external_magnetic_field()<<" "
                 <<mag_at_h2_dft/qb.external_magnetic_field()/(137.036*137.036)*4.*3.1415926*1e6<<" (dft)"<<endl;
         }
      }

      system_signal::terminal_interrupt::release();

      qb.export_history("commands.qbox");
      clog<<"finish ecut "<<ecut<<endl;

   }
   return 0;
}

