#ifndef BOOST_NUMERIC_MAKE_INTERVAL
#define BOOST_NUMERIC_MAKE_INTERVAL
namespace boost{
namespace numeric{
	template<class T>
	interval<T>
	make_interval(T const& t1, T const& t2){
		return interval<T>(t1, t2);
	}
}
}
#endif

