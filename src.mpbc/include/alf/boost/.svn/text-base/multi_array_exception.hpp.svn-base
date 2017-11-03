#ifndef MULTI_ARRAY_EXCEPTION_HPP_
#define MULTI_ARRAY_EXCEPTION_HPP_
#ifndef BOOST_ENABLE_ASSERT_HANDLER
#define RESTORE_BOOST_ENABLE_ASSERT_HANDLER_IN_MA_EXCEPTION
#endif
#define BOOST_ENABLE_ASSERT_HANDLER //enable behavior change
#ifdef BOOST_MULTI_ARRAY_RG071801_HPP
BOOST_STATIC_ASSERT(0 /*multi_array.hpp was included already with different exception handling!*/);
#endif
#include<boost/lexical_cast.hpp>
#include<stdexcept>
#include<string>
#include<boost/assert.hpp>
namespace boost{
  void assertion_failed(char const * expr, char const * function, char const * file, long line){ //define this function
    using namespace std; 
    if(string(expr)=="size_type(idx - index_bases[0]) < extents[0]"){//string may change with the implementation but in any case range_error is-a runtime_error also
      throw 
				range_error(
				//std::string(
					"multi_array index range error (according to "+
					string(file)+":"+
					boost::lexical_cast<string>(line)+
					")"
				);
    }else{
      throw 
				runtime_error(
				//std::string(
					string(file)+":"+
					boost::lexical_cast<string>(line)+
					": "+string(function)+": Assertion `"+string(expr)+"' failed."
				); //other types of errors probably associated with multiarray
    }
  }
}
#include<boost/multi_array.hpp>  //define classes but with new behavior, multi_array shouldn't be included before this
#ifdef RESTORE_BOOST_ENABLE_ASSERT_HANDLER_IN_MA_EXCEPTION
#undef BOOST_ENABLE_ASSERT_HANDLER
#endif
#endif

