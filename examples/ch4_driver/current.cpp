#if 0
  echo "compiling..."
  rm -f $0.x
  time c++ $0 -Wfatal-errors -o $0.x \
	-I../../src.wei -I../../src.wei/include \
	-Wl,-rpath=$HOME/usr/lib -Wall \
	-D_MPBC -DUSE_MPI -DUSE_HDF5 -DSCALAPACK \
	-L../../src.wei -lqb -lfftw -lfftw3 -llapack \
	-lboost_system -lboost_filesystem  -lboost_mpi -lxerces-c \
	-lmpich -lz \
	-lhdf5 $HOME/usr/lib/libscalapack.a $HOME/usr/lib/libblacsCinit-mpi.a $HOME/usr/lib/libblacs-mpi.a $HOME/usr/lib/libblacsCinit-mpi.a  -lmpl  && echo '--' && echo "done. running ... " \
	&& ./$0.x $1 $2 $3 $4 $5 $6 $7 $8 $9 && rm -f $0.x
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
  vector origin(15,15,15);

  vector xhat(1,0,0);
  vector yhat(0,1,0);
  vector mol_dir;

/*
  if (argc>1) {
    if ((string(argv[1])=="x")||(string(argv[1])=="X")) {
       clog<<"mol_dir  = x"<<endl;
       mol_dir = xhat;
    }
    else if ((string(argv[1])=="y")||(string(argv[1])=="Y")) {
       clog<<"mol_dir  = y"<<endl;
       mol_dir = yhat;
    }
    else {
       clog<<"mol_dir not assigned!"<<endl;
       return -1;
    }
  }
  else {
    clog<<"mol_dir not assigned!"<<endl;
    return -1;
  }
*/

  vector c_pos (origin);
  vector h1_pos(origin + vector( 1.2, 1.2, 1.2));
  vector h2_pos(origin + vector( 1.2,-1.2,-1.2));
  vector h3_pos(origin + vector(-1.2, 1.2,-1.2));
  vector h4_pos(origin + vector(-1.2,-1.2, 1.2));

  std::ofstream convergence("convergence.dat");
  convergence<<"# ecut, magnetic at H (a.u.SI)"<<endl;

  for(double ecut=90; ecut<=130; ecut+=10){
      clog<<"starting ecut "<<ecut<<endl;
      ofstream qbout("qb_log-"+boost::lexical_cast<std::string>(ecut)+ ".xml");
      qbox::session qb(qbout); 
      qb.set_nrowmax(1u); //one
      qb.kpoint_delete(vector(0     , 0 , 0)	);
      qb.kpoint_add   (vector(0.0001, 0., 0), 1.);
      qb.set_nempty   (1);
      qb.set_cell(cell(
			 vector(30,  0,  0),
			 vector( 0, 30,  0),
			 vector( 0,  0, 30)
      ));
      std::cout<<"magnetic field "<<qb.magnetic_field()<<std::endl;
      //qb.species("hydrogen", "hydrogen_bare_extended_withelectron.xml");
      qb.species("carbon"  , "../ch4/carbon_pbe.xml"  );
      qb.species("hydrogen", "../ch4/hydrogen_pbe.xml");

      qb.atom("C" , "carbon"  , c_pos );
      qb.atom("H1", "hydrogen", h1_pos); 
      qb.atom("H2", "hydrogen", h2_pos); 
      qb.atom("H3", "hydrogen", h3_pos); 
      qb.atom("H4", "hydrogen", h4_pos); 

      qb.export_xyz("atoms.xyz");

      qb.set_ecut(ecut);
      qb.set_xc("PBE"); 
      qb.set_wf_dyn("PSD");
      qb.set_ecutprec(ecut-5);
      qb.randomize_wf();

      clog<<"calc density..."<<endl;
      qbox::field<vector> cur__LINE__=qbox::current(qb.wavefunction());
      clog<<"..."<<endl;
      qbox::field<double> den__LINE__=qbox::density(qb.wavefunction());
      clog<<"...done"<<endl;
      ofstream ofs("etotal."+boost::lexical_cast<std::string>(ecut)+".dat");
      ofs<<"# energy "<<endl;
      ofstream ofs_eig("eigv."+boost::lexical_cast<std::string>(ecut)+".dat");
      ofs_eig<<"# eig 0, 1, ...\n";

      clog<<"begining for loop... (press ctrl-c to break)"<<endl;
      system_signal::terminal_interrupt::capture();

      int iter;
      for ( iter = 0; iter <= 50 ; iter ++ ){
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

         if ( iter % 10 == 0 && iter > 0 ) {
            clog<<"saving to wf.xml"<<endl;
            qb.save("wf.ecut"+lexical_cast<string>(ecut)+".xml");
            clog<<"saving .h5 files"<<endl;
            file("density"+lexical_cast<string>(ecut)+".h5")<<make_const_nvp("density", qbox::density(qb.wavefunction()));
            file("current"+lexical_cast<string>(ecut)+".h5")<<make_const_nvp("current", qbox::current(qb.wavefunction()));
            file("magnetic"+lexical_cast<string>(ecut)+".h5")<<make_const_nvp("magnetic", qbox::magnetic(qbox::current(qb.wavefunction())));
 
            qbox::field<double> den=qbox::density(qb.wavefunction());
            qbox::field<vector> cur=qbox::current(qb.wavefunction());
            qbox::field<vector> mag=qbox::magnetic(cur);
            convergence<<ecut<<"  C  "<<mag(c_pos )<<" external magnetic field "<<qb.magnetic_field()<<endl;
            convergence<<ecut<<"  H1 "<<mag(h1_pos)<<" external magnetic field "<<qb.magnetic_field()<<endl;
            convergence<<ecut<<"  H2 "<<mag(h2_pos)<<" external magnetic field "<<qb.magnetic_field()<<endl;
            convergence<<ecut<<"  H3 "<<mag(h3_pos)<<" external magnetic field "<<qb.magnetic_field()<<endl;
            convergence<<ecut<<"  H4 "<<mag(h4_pos)<<" external magnetic field "<<qb.magnetic_field()<<endl;
         }
      }

      system_signal::terminal_interrupt::release();

      qb.export_history("commands.qbox");
      clog<<"finish ecut "<<ecut<<endl;

   }
   return 0;
}

