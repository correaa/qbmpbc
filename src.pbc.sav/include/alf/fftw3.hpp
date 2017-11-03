// if using with fftw2 include this after fftw.h (from version 2)
#include<complex>
namespace fftw3{ //namespace can also avoid clashing with fftw (fftw2)
#include<fftw3.h>
}
#include <boost/multi_array.hpp>
#include <iostream>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/utility/enable_if.hpp>

//#include "fftw3/allocator.hpp"
#include "boost/multi_array_utility.hpp"
#include "counted.hpp"

using std::clog;
using std::endl;

namespace debug{
  template<class T>
  class lifetime{
  protected:
    lifetime(){
      clog<<"object of type "<<typeid(T).name()<<" created"<<endl;
    }
    ~lifetime(){
      clog<<"object of type "<<typeid(T).name()<<" destroyed"<<endl;
    }
  };
}

namespace fftw3{
  int const forward=FFTW_FORWARD;
  int const backward=FFTW_BACKWARD;
  using namespace boost;
  class plan : private counted<plan>{//fftw3::plan p(in, out, direction, flag);
    plan(fftw_plan const& p) : p_(p) {}
  public:  // templated constructor will handle N-dimensional array, other cases will be simply overload
    template<class ConstMultiArray, class MultiArray>
    plan(ConstMultiArray const& in, MultiArray const& out, int direction, unsigned flags=FFTW_ESTIMATE){
      assert(boost::shape(in)==boost::shape(out));
      //if(in.num_dimensions()>1){
      // check that arrays are dense
	for(int i=in.num_dimensions()-1, dense_in=in.strides()[i], dense_out=out.strides()[i]; i!=-1; dense_in*=in.shape()[i], dense_out*=out.shape()[i], --i){
	  assert(dense_in==in.strides()[i] and dense_out==out.strides()[i]); // check for dense representation of multidimensional
	}
	//}
      switch(in.num_dimensions()){
      /*
	case 1:
	p_ = fftw_plan_dft_1d(in.shape()[0], 
			      (double(*)[2])in.origin(), 
			      (double(*)[2])out.origin(),
			      direction, flags);
      case 2:
	p_ = fftw_plan_dft_2d(in.shape()[0], in.shape()[1],
			      (double(*)[2])in.origin(), 
			      (double(*)[2])out.origin(),
			      direction, flags);
	break;
      case 3:
	p_ = fftw_plan_dft_3d(in.shape()[0], in.shape()[1], in.shape()[2],
			      (double(*)[2])in.origin(), 
			      (double(*)[2])out.origin(),
			      direction, flags);
	break;
      */
      default:
	p_ = fftw_plan_many_dft
	  (in.num_dimensions(), &shape_int(in)[0], 1 /*int howmany*/,
	 (double(*)[2])in.origin() , 0 /*const int *inembed, 0=in.shape()*/,
	 in.strides()[in.num_dimensions()-1] /*int istride*/, 0 /*int idist ignored for howmany=1*/,
	 (double(*)[2])out.origin(), 0 /*const int *onembed, 0=out.shape()*/,
	 out.strides()[out.num_dimensions()-1] /*int ostride*/, 0 /*int odist ignored for howmany=1*/,
	 direction, flags);
	break;
      }
    }
    void operator()() const{return fftw_execute(p_);}
    virtual ~plan(){fftw_destroy_plan(p_);}
  private:
    fftw_plan p_;
  public:
    template<class ConstMultiArray, class MultiArray>
    static plan many(ConstMultiArray const& in, MultiArray const& out, int direction, unsigned flags=FFTW_ESTIMATE){
      assert(in.num_dimensions()>1);
      assert(boost::shape(in)==boost::shape(out));
      for(int i=in.num_dimensions()-1, dense_in=in.strides()[i], dense_out=out.strides()[i]; i!=-1; dense_in*=in.shape()[i], dense_out*=out.shape()[i], --i){
	assert(dense_in==in.strides()[i] and dense_out==out.strides()[i]); // check for dense representation of multidimensional
      }
      return plan
	(fftw_plan_many_dft
	 (in.num_dimensions()-1, (const int *)&(in.shape()[1]), in.shape()[0] /*int howmany*/,
	  (double(*)[2])in.origin() , 0 /*const int *inembed, 0=in.shape()*/,
	  1 /*int istride*/, in.num_elements()/in.shape()[0] /*int idist ignored for howmany=1*/,
	  (double(*)[2])out.origin(), 0 /*const int *onembed, 0=out.shape()*/,
	  1 /*int ostride*/, out.num_elements()/out.shape()[0] /*int odist ignored for howmany=1*/,
	  direction, flags)
	 );
    }
    template<class ConstMultiArray, class MultiArray>
    static plan many(unsigned tdim, ConstMultiArray const& in, MultiArray const& out, int direction, unsigned flags=FFTW_ESTIMATE){
      assert(in.num_dimensions()>1);
      assert(tdim<=in.num_dimensions());
      assert(boost::shape(in)==boost::shape(out));
      int howmany=1;
      for(unsigned i=0; i!=in.num_dimensions()-tdim; ++i)
	howmany*=in.shape()[i];
      for(int i=in.num_dimensions()-tdim, dense_in=in.strides()[i], dense_out=out.strides()[i]; 
	  i!=-1; dense_in*=in.shape()[i], dense_out*=out.shape()[i], --i){
	assert(dense_in==in.strides()[i] and dense_out==out.strides()[i]); // check for dense representation of multidimensional
      }
      return plan
	(fftw_plan_many_dft
	 (tdim, (const int *)&(in.shape()[in.num_dimensions()-tdim]), howmany /*int howmany*/,
	  (double(*)[2])in.origin() , 0 /*const int *inembed, 0=in.shape()*/,
	  1 /*int istride*/, in.num_elements()/howmany /*int idist ignored for howmany=1*/,
	  (double(*)[2])out.origin(), 0 /*const int *onembed, 0=out.shape()*/,
	  1 /*int ostride*/, out.num_elements()/howmany /*int odist ignored for howmany=1*/,
	  direction, flags)
	 );
    }
    static void cleanup(){
      if(get_count()<=0) throw std::logic_error("trying to call fftw3::cleanup() with "+boost::lexical_cast<std::string>(get_count())+" active plans, should be zero");
      fftw_cleanup();
    }
  };
  void execute(plan const& p) {p.operator()();}
}

namespace fftw3{ 
  template<class ConstMultiArray, class MultiArray>
  void dft(ConstMultiArray const& in, MultiArray const& out, int direction){ //need const in out to use with temporaries (marrays views).
    plan p(in, out, direction, FFTW_ESTIMATE);
    execute(p);
  }
  template<class ConstMultiArray, class MultiArray>
  void dft_many(unsigned tdim, ConstMultiArray const& in, MultiArray const& out, int direction){
    plan p=plan::many(tdim, in, out, direction, FFTW_ESTIMATE);
    execute(p);
  }
  template<class ConstMultiArray, class MultiArray>
  void normalized_dft_many(unsigned tdim, ConstMultiArray const& in, MultiArray& out, int direction){
    plan p=plan::many(tdim, in, out, direction, FFTW_ESTIMATE);
    execute(p);
    int howmany=1;
    for(unsigned i=0; i!=in.num_dimensions()-tdim; ++i){
      howmany*=in.shape()[i];
    }
    double normalization=sqrt(out.num_elements()/(double)howmany);
    for(std::complex<double>* i=out.origin(); i!=out.origin()+out.num_elements(); ++i)
      (*i)/=normalization;
  }
  template<class ConstMultiArray, class MultiArray>
  void normalized_dft(ConstMultiArray const& in, MultiArray const& out, int direction){
    dft(in, out, direction);
    for(std::complex<double>* i=out.origin(); i!=out.origin()+out.num_elements(); ++i)
      (*i)/=sqrt((double)out.num_elements());
  }
}
