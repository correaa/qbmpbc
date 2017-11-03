#ifndef _MULTI_ARRAY_UTILITY_HPP_
#define _MULTI_ARRAY_UTILITY_HPP_
/* extentions and utility functions for MultiArray */
#include<boost/multi_array.hpp>
#include "boost/cast.hpp"
//#include "alf/boost/array_io.hpp"

namespace boost{
  /* important:
   * never use shape for construction (it is silently allowed but give erroneous results), use ranges instead (defined below)
   */
  template<class MultiArray>
  boost::array<typename MultiArray::size_type, MultiArray::dimensionality> 
  shape(MultiArray const& ma){
    boost::array<typename MultiArray::size_type, MultiArray::dimensionality> ret;
    std::copy(ma.shape(), ma.shape() + ma.num_dimensions(), ret.begin());
    return ret;
  }
  template<class MultiArray>
  boost::array<int, MultiArray::dimensionality>
  shape_int(MultiArray const& ma){
    boost::array<int, MultiArray::dimensionality> ret;
    ret=shape(ma);
    return ret;
  }
	
//  template<class String, class Size_t, unsigned N>
//  String lexical_cast(boost::array<Size_t,N> const& a){
//    std::ostringstream oss;
//	oss<<a;
//	oss<<"{";
//    std::copy(a.begin(), a.end(), std::ostream_iterator<Size_t>(oss, ","));
//    oss<<"}";
//    return oss.str();
//  };
  
  //std::string lexical_cast(boost::array<unsigned
template<class MultiArray> 
  struct shape_of{
    typedef boost::array<typename MultiArray::size_type, MultiArray::dimensionality> type;
  };
  
  template<class MultiArray>
  typename shape_of<MultiArray>::type multi_array_shape(MultiArray const& M){
    typename shape_of<MultiArray>::type ret; assert(ret.size()==M.num_dimensions());
    std::copy(M.shape(), M.shape() + M.num_dimensions(), ret.begin());
    return ret;
  }
  
  // convert the following three function to multi_array_ranges
  template <typename T, std::size_t NumDims> 
  detail::multi_array::extent_gen<NumDims> ranges(const_multi_array_ref<T, NumDims> const& M){ //this function is adapted from 
    typedef detail::multi_array::extent_gen<NumDims> gen_type;
    gen_type ranges;												 //note my new convention return variables are named after function
    for (int i=0; i != NumDims; ++i){
      typedef typename gen_type::range range_type;
      ranges.ranges_[i] = range_type(M.index_bases()[i],M.index_bases()[i]+M.shape()[i]);
    }
    return ranges;
  }
  template <typename T, std::size_t NumDims> inline
  detail::multi_array::extent_gen<NumDims> ranges(multi_array_ref<T, NumDims> const& M){
    return ranges(static_cast<const_multi_array_ref<T, NumDims> const&>(M));
  }
  template <typename T, std::size_t NumDims> inline
  detail::multi_array::extent_gen<NumDims> ranges(multi_array<T, NumDims> const& M){
    return ranges(static_cast<multi_array_ref<T, NumDims> const&>(M));
  }	
  
  template<typename T, std::size_t NumElem, typename T2> // T2 should be convertible to T 
  boost::array<T, NumElem+1>
  copy_append(boost::array<T, NumElem> const& source, T2 const& t){
    boost::array<T, NumElem+1> ret; 
    std::copy(source.begin(), source.end(), ret.begin());
    ret[NumElem]=t;
    return ret;
  }
  
  
  /*
    template<typename Element>
    struct operator_bracket_dispatch{
    template<class Array>
    static wrap_array_view<MultiArray::value_type> go(Array& ma, int i){
    return wrap_array_view<MultiArray::value_type>(ma[i]);
    }
    static
    Element& go(Array& ma, int i){
    return ma[i];
    }
    };
  */
  
  template<class MultiArray> class wrap_array_view;
  template<typename Element, typename Value>
  struct return_helper{
    typedef wrap_array_view<Value> type;
  };
  template<class MultiArray> class wrap_array_view;
  template<typename Element, typename Value>
  struct return_helper_2_args{
    typedef wrap_array_view<Value> type;
  };
  template<typename Element>
  struct ignore_second_argument_ctr{
    //Element& value_; //this line doesnt work with gcc 4.1
    Element value_;
    //    ignore_second_argument_ctr(Element& val, bool*) : value_(val){}
    ignore_second_argument_ctr(Element val, bool*) : value_(val){}
    //operator Element&(){return value_;}
    operator Element(){return value_;}
  };
  template<typename Element>
  struct return_helper<Element, Element>{
    typedef Element type;
  };
  template<typename Element>
  struct return_helper_2_args<Element, Element>{
    typedef ignore_second_argument_ctr<Element> type;
  };
  template<class MultiArrayRef>
  class wrap_array_view {
    MultiArrayRef ma_;
	typename MultiArrayRef::index* index_bases_;
	bool wrapping_[MultiArrayRef::dimensionality];
  public:
    wrap_array_view(MultiArrayRef const&ma, bool const* wp=0) : ma_(ma), index_bases_(new typename MultiArrayRef::index [MultiArrayRef::dimensionality]){
	  for(unsigned i=0; i!=MultiArrayRef::dimensionality; ++i){
		BOOST_ASSERT(ma_.shape()[i]%2==0);
		BOOST_ASSERT(ma_.index_bases()[i]==0);
		if(wp!=0){
			wrapping_[i]=wp[i];
			if(wp[i]==false){
				index_bases_[i]=0;
			}else{
				index_bases_[i]=-boost::numeric_cast<typename MultiArrayRef::index>(ma_.shape()[i])/2;
			}
		}else{
			wrapping_[i]=true;
			index_bases_[i]=-boost::numeric_cast<typename MultiArrayRef::index>(ma_.shape()[i])/2;
		}
	  }
	};
    typedef typename return_helper<
      typename MultiArrayRef::element&, 
      typename MultiArrayRef::reference>::type return_helper_type;
    return_helper_type operator[](int i){
	  if(wrapping_[0]==true){
	      if(i<ma_.index_bases()[0]) i+=ma_.shape()[0];
	  }
      return typename return_helper_2_args<
      	typename MultiArrayRef::element&, 
      	typename MultiArrayRef::reference>::type(ma_[i], wrapping_+1);
    }
	const typename MultiArrayRef::size_type* shape() const{return ma_.shape();}
	const typename MultiArrayRef::index* index_bases() const{return index_bases_;}
	~wrap_array_view(){
		delete[] index_bases_;
	}
  };
}//namespace boost
#endif
