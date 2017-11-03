#ifndef BOOST_MAKE_ARRAY
#include<boost/array.hpp>

namespace boost{

	template<class T> boost::array<T, 0> make_array(){
		boost::array<T, 0> ret = {{}}; return ret;}
	template<class T> boost::array<T, 1> make_array(T t1){
		boost::array<T, 1> ret = {{t1}}; return ret;}
	template<class T> boost::array<T, 2> make_array(T t1, T t2){
		boost::array<T, 2> ret = {{t1, t2}}; return ret;}
	template<class T> boost::array<T, 3> make_array(T t1, T t2, T t3){
		boost::array<T, 3> ret = {{t1, t2, t3}}; return ret;}
	template<class T> boost::array<T, 4> make_array(T t1, T t2, T t3, T t4){
		boost::array<T, 4> ret = {{t1, t2, t3, t4}}; return ret;}
	template<class T> boost::array<T, 5> make_array(T t1, T t2, T t3, T t4, T t5){
		boost::array<T, 5> ret = {{t1, t2, t3, t4, t5}}; return ret;}
	template<class T> boost::array<T, 6> make_array(T t1, T t2, T t3, T t4, T t5, T t6){
		boost::array<T, 6> ret = {{t1, t2, t3, t4, t5, t6}}; return ret;}
	template<class T> boost::array<T, 6> make_array(T t1, T t2, T t3, T t4, T t5, T t6, T t7){
		boost::array<T, 7> ret = {{t1, t2, t3, t4, t5, t6, t7}}; return ret;}
	template<class T> boost::array<T, 6> make_array(T t1, T t2, T t3, T t4, T t5, T t6, T t7, T t8){
		boost::array<T, 8> ret = {{t1, t2, t3, t4, t5, t6, t7, t8}}; return ret;}
	template<class T> boost::array<T, 6> make_array(T t1, T t2, T t3, T t4, T t5, T t6, T t7, T t8, T t9){
		boost::array<T, 9> ret = {{t1, t2, t3, t4, t5, t6, t7, t8, t9}}; return ret;}
	template<class T> boost::array<T, 6> make_array(T t1, T t2, T t3, T t4, T t5, T t6, T t7, T t8, T t9, T t10){
		boost::array<T,10> ret = {{t1, t2, t3, t4, t5, t6, t7, t8, t9, t10}}; return ret;}

}

#endif
/*
use as 

 f(boost::make_array(1.0, 2.0, 3.0)); //passes boost::array<double, 3>
 using boost::make_array;
 f(make_array(1,2,3)); //passes boost::array<int, 3>
*/

