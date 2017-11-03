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

#include "alf/fftw3.hpp"
#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>

using namespace std;

void test_2()
{
  typedef std::complex<double> complex;
  using boost::multi_array;
  using boost::multi_array_ref;
  using boost::extents; 
  using fftw3::plan;

  int n0, n1;
  int np0, np1, np2;

  np0 = 32; np1 = 32; np2 = 32;
  n0 = np0 * np1 * np2; n1 = n0;

  std::vector<complex> tmp (n0);
  std::vector<complex> ctmp(n1);

  clog<<"sizeof(unsigned)"<<sizeof(unsigned)<<endl;
  clog<<"sizeof(int)"<<sizeof(int)<<endl;
  clog<<"sizeof(ma::size_type)"<<sizeof(multi_array<complex,3>::size_type)<<endl;

  multi_array_ref<complex, 3> wfxyz(&tmp[0], extents[np0][np1][np2]);
  static multi_array<complex, 3> wfXYZ(extents[np0][np1][np2]);

  fftw3::fftw_plan myplan=fftw3::fftw_plan_dft_3d(np0, np1, np2, (double (*)[2])wfxyz.origin(), (double (*)[2])wfXYZ.origin(), FFTW_FORWARD, FFTW_ESTIMATE);
  cout<<"calling fftw3 execute...\n";
  fftw3::fftw_execute(myplan);

  fftw3::fftw_plan myplan_many=fftw3::fftw_plan_many_dft
    (3, &boost::shape_int(wfxyz)[0], 1,
     (double(*)[2])wfxyz.origin(), 0,
     wfxyz.strides()[wfxyz.num_dimensions()-1], 0, 
     (double(*)[2])wfXYZ.origin(), 0, 
     wfXYZ.strides()[wfXYZ.num_dimensions()-1], 0,
     FFTW_FORWARD, FFTW_ESTIMATE);
  cout<<"calling fftw3 many execute...\n";
  fftw3::fftw_execute(myplan_many);


  plan pF_wfxyz2wfXYZ(wfxyz, wfXYZ, FFTW_FORWARD, FFTW_ESTIMATE);

  /* the following line causes segmentation fault */
  std::cout<<"calling execute...\n";
  execute(pF_wfxyz2wfXYZ);

}

int main(int argc, char **argv){

   test_2 ();

   return 0;
}


