#if 0
ln -sf $0 .$0.cpp && c++ -Wall `#-Wfatal-errors` -I$HOME/usr/include -I$HOME/prj .$0.cpp -Wl,-rpath=$HOME/usr/lib -L$HOME/usr/lib -lboost_filesystem -lboost_system -D_NAMED_HPP_TEST -o ./.$0.x && ./.$0.x $1 $2 $3 $4 $5 $6 $7 $8 $9
rm -f .$0.x .$0.cpp
exit
#endif
#ifndef BOOST_NAMED_HPP
#define BOOST_NAMED_HPP
#include<string>
#include<boost/lexical_cast.hpp>
#include "boost/just.hpp"

namespace boost{
using std::string;
template<class T>
struct named : virtual just<T>::type{
	named(string name, T const& t) : just<T>::type(t), name_(name){}
	named(T const& t, string name) : just<T>::type(t), name_(name){}
	named(T const& t)              : just<T>::type(t), name_(boost::lexical_cast<string>(this)){}
	named();
	operator named<T const&>() const{
		return named<T const&>(
			name_, 
			(typename just<T>::type const&)(*this)
		);
	}
	string const& name() const{return name_;}
	virtual ~named(){}
	protected:
	string name_;
};
template<class T>
named<T>::named() : just<T>::type(T()), name_(boost::lexical_cast<string>(this)){}

template<class T>
named<T> make_named(T& t, string name){return named<T>(t, name);}
}
#endif

#ifdef _NAMED_HPP_TEST
#include<iostream>
using std::cout; using std::endl;
struct dummy{
	int i_;
	dummy(int i) : i_(i){}
};
int main(){
	dummy d(3);
	boost::named<dummy> a(d, "dummy");
	double b = 4.;
	boost::named<double const&> bn(b, "bn");
	b += 1;
	cout << boost::named<double>(b, "kk") << endl;
	cout << boost::make_named(b, "dd") <<endl;
}
#endif

