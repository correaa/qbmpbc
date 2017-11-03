#ifndef FFTW3_HPP_
#define FFTW3_HPP_
/**
 * \brief example
 *   fftw3::plan(A, B, fftw3::forward, fftw3::estimate).execute_normalize();
 */
//if using with fftw2 include this file AFTER fftw.h (from version 2) and not the other way around
#include<complex>  //< before fftw3.h
#include<cstdio> //fix a bug in gcc 4.4
namespace fftw3{   //< this namespace is also to avoid clashing with fftw (i.e. fftw2)
#ifndef FFTW3_H_INCLUDED_FROM_FFTW3_HPP
#define FFTW3_H_INCLUDED_FROM_FFTW3_HPP
#undef FFTW_ESTIMATE
#undef FFTW_MEASURE
#include<fftw3.h>
#endif
}

#include <boost/multi_array.hpp>
#include <iostream>
#include <string>
#include "boost/array_io.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/math/special_functions/fpclassify.hpp> //for boost::math::isnan
//#include "fftw3/allocator.hpp" /* no longer needed, it is optional to the user */

#include "boost/multi_array_utility.hpp"
#include "alf/debug/counted.hpp"

using std::clog;
using std::endl;
using std::complex;

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

  class plan : private debug::counted<plan>{
  private:
    fftw_plan p_;
    int fft_size_;
    complex<double>* out_origin_;
    complex<double>* out_end_;
  public:
    plan(fftw_plan const& p) : p_(p) {}
    void operator()() const{return fftw_execute(p_);}
    void execute() const{
    	fftw_execute(p_);
    	if((boost::math::isnan)(real(*out_origin_)) or (boost::math::isnan)(imag(*out_origin_))) 
    		throw std::range_error("suspect NaN in FFTW output origin, most likely there is a NaN in the input");
    }
    void execute_normalize() const{
      fftw_execute(p_);
      double factor=1./sqrt((double)fft_size_);
      for(complex<double>* i=out_origin_; i!=out_end_; ++i){
		(*i)*=factor;
      }
    }
    virtual ~plan(){fftw_destroy_plan(p_);}
    static void cleanup(){
      if(get_count()<=0)
		throw std::logic_error("trying to call fftw3::cleanup() with "
			  					+boost::lexical_cast<std::string>(get_count())
			       				+" active plans, should be zero");
      fftw_cleanup();
    }
	template<class ConstMultiArrayOfArrays, class MultiArrayOfArrays>
	plan(
		ConstMultiArrayOfArrays const& in, 
		     MultiArrayOfArrays      & out,
		int direction, unsigned flags=FFTW_ESTIMATE
		,typename boost::enable_if_c<is_same<typename ConstMultiArrayOfArrays::element::value_type, std::complex<double> >::value, void*>::type=0
	){
	  BOOST_STATIC_ASSERT(ConstMultiArrayOfArrays::dimensionality==MultiArrayOfArrays::dimensionality);
	  BOOST_STATIC_ASSERT(ConstMultiArrayOfArrays::element::static_size==MultiArrayOfArrays::element::static_size);
	  typedef typename ConstMultiArrayOfArrays::element arr;
	  p_ = fftw_plan_many_dft(
			 in.num_dimensions(), &shape_int(in)[0], /*int howmany*/ arr::static_size ,
	   		 /* in*/ (double(*)[2])(in.origin()->data()), /*const int *inembed*/ 0,
	   		 in.strides()[in.num_dimensions()-1]*arr::static_size /*int istride*/, 
	  		 /*int idist*/ 1 /*ignored if howmany=1*/,
	   		 /*out*/ (double(*)[2])(out.origin()->data()), /*const int *onembed*/ 0,
	   		 out.strides()[out.num_dimensions()-1]*arr::static_size /*int ostride*/, 
	   	     /*int odist*/ 1 /*ignored if howmany=1*/,
	    	 direction, flags
	  );
	  fft_size_=out.num_elements();
	  out_origin_=out.origin()->c_array();
	  out_end_=(out.origin()+out.num_elements())->c_array();
	}
	
    // full templated constructor for N-dimensional array, FFT on every dimension
    template<class ConstMultiArrayComplex, class MultiArrayComplex>
    plan(
	 ConstMultiArrayComplex const& in, MultiArrayComplex& out, 
	 int direction, unsigned flags=FFTW_ESTIMATE,
	 typename enable_if_c<is_same<typename ConstMultiArrayComplex::element, std::complex<double> >::value, void*>::type=0
	){
		BOOST_STATIC_ASSERT((is_same<typename ConstMultiArrayComplex::element, typename MultiArrayComplex::element>::value));
      	if(boost::shape(in)!=boost::shape(out)){
		throw std::runtime_error("multi_arrays for fftw3 transform mismatch size");
		//throw std::runtime_error("multi_arrays for fftw3 transform mismatch size in:"+to_string(boost::shape(in))+" out:"+ to_string(boost::shape(in)));
	  	}
      // check that arrays are dense
      for(int i=in.num_dimensions()-1, dense_in=in.strides()[i], dense_out=out.strides()[i];
          i!=-1; dense_in*=in.shape()[i], dense_out*=out.shape()[i], --i){
		// check for dense representation of multidimensional
		BOOST_ASSERT(dense_in==in.strides()[i] and dense_out==out.strides()[i]); 
      }
      switch (in.num_dimensions()){
	  /*
	  case 1:
	  p_ = fftw_plan_dft_1d(in.shape()[0], 
	  (double(*)[2])in.origin(), 
	  (double(*)[2])out.origin(),
	  direction, flags);
	  break;
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
	   in.strides()[in.num_dimensions()-1] /*int istride*/, 
	   0 /*int idist ignored for howmany=1*/,
	   (double(*)[2])out.origin(), 0 /*const int *onembed, 0=out.shape()*/,
	   out.strides()[out.num_dimensions()-1] /*int ostride*/, 
	   0 /*int odist ignored for howmany=1*/,
	   direction, flags);
		fft_size_=out.num_elements();
		out_origin_=out.origin();
		out_end_=out.origin()+out.num_elements();
		break;
      }
    }
    // special cases of templated constructor for 3d or 2d, 
    // FFT on either first 2 dimensions (for 3d array) or last 1 dimension (for 2d array)
    // if trans_dim = +2, it means transform the first 2 dimensions
    // if trans_dim = -1, it means transform the last  1 dimension

    template<class ConstMultiArray, class MultiArray>
    plan (int trans_dim,
          ConstMultiArray const& in, 
		  MultiArray& out, 
          int direction, unsigned flags=FFTW_ESTIMATE){
      BOOST_ASSERT(boost::shape(in)==boost::shape(out));
      // check that arrays are dense
      for(int i=in.num_dimensions()-1, dense_in=in.strides()[i], dense_out=out.strides()[i];
          i!=-1; dense_in*=in.shape()[i], dense_out*=out.shape()[i], --i){
	// check for dense representation of multidimensional
	BOOST_ASSERT(dense_in==in.strides()[i] and dense_out==out.strides()[i]); 
      }
      switch (in.num_dimensions()){
      case 2: /* 2d array, transform the last one dimension */
	BOOST_ASSERT(trans_dim == -1 );
#if 0 // debug
	clog << "2d array transform last one dimension" << endl;
	clog << "fftw_plan_many_dft( rank=" << abs(trans_dim)
	     << ", n=" << shape_int(in)[1] << ", howmany=" << shape_int(in)[0]
	     << ", in.origin(), inembed=NULL" 
	     << ", istride=" << in.strides()[in.num_dimensions()-1]
	     << ", idist=" << shape_int(in)[1]
	     << ", out.origin(), onembed=NULL"
	     << ", ostride=" << out.strides()[out.num_dimensions()-1]
	     << ", odist=" << shape_int(out)[1] 
	     << ", dir=" << direction << ", flags=" << flags << ")" << endl;
#endif
	p_ = fftw_plan_many_dft(
				abs(trans_dim), &shape_int(in)[1], shape_int(in)[0], /* rank, *n, howmany */
				(double(*)[2])in.origin(), NULL,                     /* in, inembed */
				in.strides()[in.num_dimensions()-1],                 /* istride */ 
				shape_int(in)[1],                                    /* idist */ 
				(double(*)[2])out.origin(), NULL,                    /* out, onembed */
				out.strides()[out.num_dimensions()-1],               /* ostride*/ 
				shape_int(out)[1],                                    /* odist */ 
				direction, flags);
	fft_size_=in.shape()[1];
	out_origin_=out.origin();
	out_end_=out.origin()+out.num_elements();
	break;
      case 3: /* 3d array, transform the first two dimensions */
	BOOST_ASSERT(trans_dim == +2 );
#if 0 // debug
	clog << "3d array transform first two dimensions" << endl;
	clog << "fftw_plan_many_dft( rank=" << abs(trans_dim)
	     << ", n=(" << shape_int(in)[0]<<","<<shape_int(in)[1]<<")"
	     << ", howmany=" << shape_int(in)[2]
	     << ", in.origin(), inembed=NULL" 
	     << ", istride=" << in.strides()[in.num_dimensions()-2]
	     << ", idist=" << 1
	     << ", out.origin(), onembed=NULL"
	     << ", ostride=" << out.strides()[out.num_dimensions()-2]
	     << ", odist=" << 1 
	     << ", dir=" << direction << ", flags=" << flags << ")" << endl;
#endif
	p_ = fftw_plan_many_dft
	  (abs(trans_dim), &shape_int(in)[0], shape_int(in)[2], /* rank, *n, howmany */
	   (double(*)[2])in.origin(), NULL,                     /* in, inembed */
	   in.strides()[in.num_dimensions()-2],                 /* istride */ 
	   1,                                                   /* idist */ 
	   (double(*)[2])out.origin(), NULL,                    /* out, onembed */
	   out.strides()[out.num_dimensions()-2],               /* ostride*/ 
	   1,                                                   /* odist */ 
	   direction, flags);
	fft_size_=in.shape()[0]*in.shape()[1];
	out_origin_=out.origin();
	out_end_=out.origin()+out.num_elements();
	break;
      default:
	BOOST_ASSERT(0);
	break;
      }
    }
  };// class plan
  /* need to disable function execute() to remove error in compilation
   *  function used in two places in the code creates conflicts 
   */
  //void execute(plan const& p){ p.operator()(); }
}//namespace fftw3
#endif
