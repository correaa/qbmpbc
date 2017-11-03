#ifndef BOOST_TYPE_TRAITS_REMOVE_CONST_REFERENCE
#define BOOST_TYPE_TRAITS_REMOVE_CONST_REFERENCE
#include<boost/type_traits/remove_const.hpp>
#include<boost/type_traits/remove_reference.hpp>
namespace boost{
	template<typename T>
	struct remove_const_reference{
		typedef
			typename boost::remove_const<
				typename boost::remove_reference<
					T
				>::type
			>::type
			type;
	};
}
#endif


