#ifdef compile_instructions
ln -sf $0 .$0.cpp && c++ -std=c++0x .$0.cpp `#-Wfatal-errors` -L$HOME/lib `pkg-config --libs gsl` -D_TEST_DERIVATIVE_PARTIAL_HPP -o ./.$0.x && ./.$0.x $@
rm -f .$0.x .$0.cpp
exit
#endif
// support for partial derivative in named functions
#include "../../gsl/derivative.hpp"
#include <boost/static_assert.hpp>
#include <boost/fusion/container/map.hpp>
#include <boost/fusion/include/at_key.hpp>
#include <boost/fusion/sequence/intrinsic/value_at_key.hpp> 
#include <boost/fusion/container/generation/make_map.hpp>
#include <boost/fusion/container/generation/make_set.hpp>
//#include <boost/fusion/container/generation/make_vector.hpp>
//#include <boost/fusion/algorithm/transformation/push_front.hpp>
#include <boost/fusion/functional/adapter/unfused.hpp>
#include <boost/fusion/functional/adapter/fused.hpp>
#include <boost/fusion/functional/invocation/invoke_function_object.hpp> 
struct p{};
struct q{};

namespace gsl{
namespace derivative{

using namespace boost::fusion;
template<
	class FusionMap, // e.g. fusion::map<T1,T2,...Tn>
	class Functor,   // double(FreeFunction)(FusionMap const&),
	class FreeArgumentMapKey
> 
struct bind_free_at{
	mutable FusionMap m;
	Functor f; 
	bind_free_at(Functor const& f, FusionMap const& fm) :  m(fm), f(f){}
	//bind_free_at(bind_free_at const& other) :  FusionMap((FusionMap const&)other), f(other.f){}
	double operator()(typename result_of::value_at_key<FusionMap, FreeArgumentMapKey>::type const& free_value) const{
		at_key<FreeArgumentMapKey>(m) = free_value;
		return f(m);
	}
};

template<class ParameterKey, class Functor, class Args=typename Functor::mapped_arguments> //, class FusionMap>
struct d_{
	Functor f;
	d_(Functor const& f) : f(f){}
	template <class Seq>
	struct result{typedef double type;};
	typedef typename Functor::mapped_arguments mapped_arguments;
	double operator()(mapped_arguments args) const{
		bind_free_at<Args, Functor, ParameterKey> bf(f, args);
		double x = at_key<ParameterKey>(args);
		return gsl::derivative::central(bf, x);
	}
	template<class... DomainTypeS>
	double operator()(DomainTypeS const& ... domainS) const{
		//BOOST_STATIC_ASSERT((sizeof...(DomainTypeS))>1);
		return (*this)(mapped_arguments(domainS ... ));
	}
};

template<class ParameterKey, class Functor, class Args=typename Functor::mapped_arguments>
d_<ParameterKey, Functor, Args> partial(Functor const& f){
	return d_<ParameterKey, Functor, Args>(f);
};

}
}


#ifdef _TEST_DERIVATIVE_PARTIAL_HPP
#include<iostream>
using std::clog;
using std::endl;

using namespace boost::fusion;


struct H{
	double operator()(double p, double q) const{
		return p*p + q*q;
	}
	typedef double result_type;
	typedef result_of::make_map<
			p     , q     , 
			double, double
		>::type mapped_arguments;
	result_type operator()(mapped_arguments const& args) const{
			//return invoke_function_object(*this, args);
			return operator()(at_key<p>(args), at_key<q>(args));
	}
};

int main (void){
	double uno=1., dos =2.;
	H h;
	clog << h(1.,2.) << endl;
	clog << h(make_map<p,q>(1.,2)) << endl;
	clog << gsl::derivative::partial<q>(h)(make_map<p,q>(uno,dos)) << endl;
	//clog << gsl::derivative::partial<q>(h)(uno,dos) << endl;

	using namespace gsl::derivative;
	clog << partial<q>(partial<q>(h))(1.,2) << endl;

	using namespace gsl::derivative;
	gsl::derivative::d_<q, H> const dd(gsl::derivative::partial<q>(h));
	//gsl::derivative::d_impl<q, H> dd(gsl::derivative::partial<q>(h));
	//clog << gsl::derivative::partial<q>(h)(make_map<p,q>(uno, dos)) << endl;
	//clog << gsl::derivative::partial<q>(gsl::derivative::partial<q>(h))(make_map<p,q>(1.,2)) << endl;
	//gsl::derivative::partial<q>(gsl::derivative::partial<q>(h))(make_map<p,q>(uno, dos));
	return 0;
}

#endif


