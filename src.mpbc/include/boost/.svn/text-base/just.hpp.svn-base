#if 0
ln -sf $0 .$0.cpp && c++ -Wall `#-Wfatal-errors` -I$HOME/usr/include -I$HOME/prj .$0.cpp -Wl,-rpath=$HOME/usr/lib -L$HOME/usr/lib -lboost_filesystem -lboost_system -D_TEST_BOOST_JUST_HPP -o ./.$0.x && ./.$0.x $1 $2 $3 $4 $5 $6 $7 $8 $9
rm -f .$0.x .$0.cpp
exit
#endif
#ifndef BOOST_JUST_HPP
#define BOOST_JUST_HPP
#include<boost/ref.hpp>
namespace boost{
template<class T>
struct just /*: T*/{
	typedef T type;
};
template<class T>
just<T> _(T const& t){
	return just<T>(t);
}

template<class T>
struct just<T&> /*:boost::reference_wrapper<T> */{
	typedef boost::reference_wrapper<T> type;
};
template<>
struct just<double>{
	double impl_;
	typedef just<double> type;
	just(double const& d) : impl_(d){}
	operator double const&() const{return impl_;}
	double& operator+=(double const& d){return impl_+=d;}
};
}
#endif

#define BOOST_INHERIT_UNARY_CONSTRUCTOR(MyclasS, MybaseclasS) \
	template<typename A> MyclasS(A& arg) : MybaseclasS(arg) {}

#ifdef _TEST_BOOST_JUST_HPP
struct number : boost::just<double&>::type{
	BOOST_INHERIT_UNARY_CONSTRUCTOR(number, boost::just<double&>::type)
};

#include<iostream>
int main(){
	double d=5;
	number n(d);
	n+=4;
	std::cout << n << std::endl;
}
#endif

