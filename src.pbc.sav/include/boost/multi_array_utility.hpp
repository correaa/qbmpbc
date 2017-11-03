/* extentions and utility functions for MultiArray */

namespace boost{
  /* important:
   * never use shape for construction (it is silently allowed), use ranges instead (Defined below)
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
	
  template<class String, unsigned N>
  String lexical_cast(boost::array<unsigned,N> const& a){
    std::ostringstream oss("{");
    std::copy(a.begin(), a.end(), std::ostream_iterator<unsigned>(oss, ", "));
    oss<<"}";
    return oss.str();
  };
  
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
  template <typename T, unsigned NumDims> 
  detail::multi_array::extent_gen<NumDims> ranges(const_multi_array_ref<T, NumDims> const& M){ //this function is adapted from 
    typedef detail::multi_array::extent_gen<NumDims> gen_type;
    gen_type ranges;												 //note my new convention return variables are named after function
    for (int i=0; i != NumDims; ++i){
      typedef typename gen_type::range range_type;
      ranges.ranges_[i] = range_type(M.index_bases()[i],M.index_bases()[i]+M.shape()[i]);
    }
    return ranges;
  }
  template <typename T, unsigned NumDims> inline
  detail::multi_array::extent_gen<NumDims> ranges(multi_array_ref<T, NumDims> const& M){
    return ranges(static_cast<const_multi_array_ref<T, NumDims> const&>(M));
  }
  template <typename T, unsigned NumDims> inline
  detail::multi_array::extent_gen<NumDims> ranges(multi_array<T, NumDims> const& M){
    return ranges(static_cast<multi_array_ref<T, NumDims> const&>(M));
  }	
}//namespace boost
