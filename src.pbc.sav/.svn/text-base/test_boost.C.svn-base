////////////////////////////////////////////////////////////////////////////////
//
//
//  test_fftw3_boost.C (from SlaterDet.C)
//
//  /opt/mpich/gnu/bin/mpicxx -I./include -I/home/caiwei/usr/include test_fftw3_boost.C -L/home/caiwei/usr/lib -lfftw3 -lboost_system-gcc34-mt -lboost_filesystem-gcc34-mt -o test_fftw3_boost
//
//
////////////////////////////////////////////////////////////////////////////////

#include <vector>

#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>

//--------------------------------------
//#include "alf/fftw3.hpp"
//--------------------------------------
#include <iostream>
#include "boost/multi_array_utility.hpp"
#include "counted.hpp"

namespace fftw3{
   #include <fftw3.h>
   int const forward=FFTW_FORWARD;
   int const backward=FFTW_BACKWARD;
   using namespace boost;
   class plan : private counted<plan>{
       //fftw3::plan p(in, out, direction, flag);
       plan(fftw_plan const& p) : p_(p) {}
       public:
       // templated constructor will handle N-dimensional array, other cases will be simply overload
           template<class ConstMultiArray, class MultiArray>
	   plan(ConstMultiArray const& in, MultiArray const& out, int direction, unsigned flags=FFTW_ESTIMATE){
	     assert(boost::shape(in)==boost::shape(out));
	     if(in.num_dimensions()>1){
	       for(int i=in.num_dimensions()-1, dense_in=in.strides()[i], dense_out=out.strides()[i]; 
                   i!=-1; dense_in*=in.shape()[i], dense_out*=out.shape()[i], --i){
	           assert(dense_in==in.strides()[i] and dense_out==out.strides()[i]); 
                   // check for dense representation of multidimensional
               }
	     }
             std::cout<<"fftw_plan_many_dft("<<in.num_dimensions()<<","
                 << &shape_int(in)[0]<<","
                 << 1 /*int howmany*/ <<","
                 << (double(*)[2])in.origin() <<","
                 << 0 /*const int *inembed, 0=in.shape()*/ <<","
                 << in.strides()[in.num_dimensions()-1] /*int istride*/ <<"," 
                 << 0 /*int idist ignored for howmany=1*/ <<","
                 << (double(*)[2])out.origin() <<","
                 << 0 /*const int *onembed, 0=out.shape()*/ <<","
                 << out.strides()[out.num_dimensions()-1] /*int ostride*/ <<"," 
                 << 0 /*int odist ignored for howmany=1*/ <<","
                 << direction<<","<< flags<<")\n";

#if 1
             p_ = fftw_plan_many_dft(in.num_dimensions(), &shape_int(in)[0], 1 /*int howmany*/,
                                  (double(*)[2])in.origin() , 0 /*const int *inembed, 0=in.shape()*/,
                                  in.strides()[in.num_dimensions()-1] /*int istride*/, 
                                  0 /*int idist ignored for howmany=1*/,
                                  (double(*)[2])out.origin(), 0 /*const int *onembed, 0=out.shape()*/,
                                  out.strides()[out.num_dimensions()-1] /*int ostride*/, 
                                  0 /*int odist ignored for howmany=1*/,
                                  direction, flags);
#endif

#if 0
             if (in.num_dimensions()==1)
                 p_ = fftw_plan_dft_1d(in.shape()[0], 
                                  (double(*)[2])in.origin(), (double(*)[2])out.origin(), 
                                  direction, flags);
             if (in.num_dimensions()==2)
                 p_ = fftw_plan_dft_2d(in.shape()[0], in.shape()[1], 
                                  (double(*)[2])in.origin(), (double(*)[2])out.origin(), 
                                  direction, flags);
             if (in.num_dimensions()==3)
                 p_ = fftw_plan_dft_3d(in.shape()[0], in.shape()[1], in.shape()[2],
                                  (double(*)[2])in.origin(), (double(*)[2])out.origin(), 
                                  direction, flags);
#endif
           }
           void operator()() const{
                 std::cout<<"fftw_execute(p_)\n";
                 return fftw_execute(p_);
           }
           virtual ~plan(){fftw_destroy_plan(p_);}
        private:
           fftw_plan p_;
   };
   void execute(plan const& p) {p.operator()();}
}


void test_boost()
{
   int i, j, k;
   int np0, np1, np2;
   typedef std::complex<double> complex;

   np0 = 4; np1 = 4; np2 = 4;
   boost::multi_array<complex, 2> A(boost::extents[np0][np1]);
   boost::multi_array<complex, 2> B(boost::extents[np0][np1]);
   boost::multi_array<complex, 2> C(boost::extents[np0][np1]);

   boost::multi_array<complex, 3> R(boost::extents[np0][np1][np2]);
   boost::multi_array<complex, 3> G(boost::extents[np0][np1][np2]);
   boost::multi_array<complex, 3> S(boost::extents[np0][np1][np2]);

   for(i=0; i<R.shape()[0]; i++)
      for(j=0; j<R.shape()[1]; j++)
            A[i][j]=complex(i+j, i-j);

   for(i=0; i<R.shape()[0]; i++)
      for(j=0; j<R.shape()[1]; j++)
         for(k=0; k<R.shape()[2]; k++)
            R[i][j][k]=complex(i+j+k, i-j-k);

   fftw3::plan A_to_B (A, B, FFTW_FORWARD,  FFTW_ESTIMATE);
   fftw3::plan B_to_C (B, C, FFTW_BACKWARD, FFTW_ESTIMATE);
   fftw3::plan R_to_G (R, G, FFTW_FORWARD,  FFTW_ESTIMATE);
   fftw3::plan G_to_S (G, S, FFTW_BACKWARD, FFTW_ESTIMATE);

   std::cout<<"calling execute...\n";
   fftw3::execute( A_to_B );
   fftw3::execute( B_to_C );
   fftw3::execute( R_to_G );
   fftw3::execute( G_to_S );

   for(i=0; i<R.shape()[0]; i++)
      for(j=0; j<R.shape()[1]; j++)
            std::cout<<"A ["<<i<<"]["<<j<<"] = "<<A[i][j]<<"\t"
                     <<"B ["<<i<<"]["<<j<<"] = "<<B[i][j]<<"\t"
                     <<"C ["<<i<<"]["<<j<<"] = "<<C[i][j]<<"\n";

   for(i=0; i<R.shape()[0]; i++)
      for(j=0; j<R.shape()[1]; j++)
         for(k=0; k<1; k++)
            std::cout<<"R ["<<i<<"]["<<j<<"]["<<k<<"] = "<<R[i][j][k]<<"\t"
                     <<"G ["<<i<<"]["<<j<<"]["<<k<<"] = "<<G[i][j][k]<<"\t"
                     <<"S ["<<i<<"]["<<j<<"]["<<k<<"] = "<<S[i][j][k]<<"\n";
}



int main(int argc, char **argv){

   test_boost ();

   return 0;
}


