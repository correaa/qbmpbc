#if 0
	cp $0 $0.cpp
	c++ \
	-Wall -std=c++0x -I$HOME/usr/include $0.cpp -L$HOME/usr/lib -Wl,-rpath=$HOME/usr/lib -llapack -lptf77blas -latlas -lgfortran -lpthread -lboost_system -o $0.cpp.x -D_TEST_LINEAR_SPACE_HPP && ./$0.cpp.x
	rm -f $0.cpp
	exit
#endif

#include<boost/multi_array.hpp>
#include<boost/operators.hpp>
#include<boost/type_traits/function_traits.hpp> 
#include<iostream>
using std::clog; using std::endl;
// useful for std::vector and boost::array

template<class T>
struct elementwise_addition{
	static T apply(T const& t1, T const& t2){
		BOOST_ASSERT(t1.size()==t2.size());
		T ret(t1);
		for(unsigned i=0; i!=ret.size(); ++i){
			ret[i]+=t2[i];
		}
		return ret;
	}
	T operator()(T const& a1, T const& a2) const{
		BOOST_ASSERT(a1.size()==a2.size());
		T ret(a1);
		for(unsigned i=0; i!=ret.size(); ++i){
			ret[i]+=a2[i];
		}
		return ret;
	}
};

template<typename T, unsigned D>
struct elementwise_addition<boost::multi_array<T, D> >{
	typedef boost::multi_array<T, D> args;
	static boost::multi_array<T, D> apply(args& a1, args const& a2){
		args ret(a1);
		for(unsigned i=0; i!=a1.num_elements(); ++i){
			ret.origin()[i]+=a2.origin()[i];
		}
		return ret;
	}
};

template<typename S, class T>
struct elementwise_multiplication{
	static T apply(S const& f, T& a){
		T ret(a);
		for(unsigned i=0; i!=ret.size(); ++i){
			ret[i]*=f;
		}
		return ret;
	}
	T operator()(S const& f, T& a) const{
		T ret(a);
		for(unsigned i=0; i!=ret.size(); ++i){
			ret[i]*=f;
		}
		return ret;
	}
};

template<class T>
struct elementwise_inner_product{
	static double conj(double const& d){return d;}
	static typename T::value_type apply(T const& t1, T const& t2){
		//assert(t1.size()==t2.size());
		typename T::value_type ret(0);
		for(unsigned i=0; i!=t1.size(); ++i){
			ret+=(t1[i]*conj(t2[i]));
		}
		return ret;
	}
};

template<class T, unsigned D>
struct elementwise_inner_product<boost::multi_array<T, D> >{	
	static T apply(boost::multi_array<T, D> const& t1, boost::multi_array<T, D> const& t2){
		T ret(0);
		for(unsigned i=0; i!=t1.num_elements(); ++i){
			ret+=t1.origin()[i]*t2.origin()[i];
		}
		return ret;
	}
};

template<typename S, class T, unsigned D>
struct elementwise_multiplication<S, boost::multi_array<T, D> >{
	static boost::multi_array<T, D> apply(S const& f, boost::multi_array<T, D>& a){
		boost::multi_array<T, D> ret(a);
		for(unsigned i=0; i!=ret.num_elements(); ++i){
			ret.origin()[i]*=f;
		}
		return ret;
	}
};

/*
template <class _Tp1, class _Tp2, class _Tp3>
struct multiplies : public std::binary_function<_Tp1, _Tp2, _Tp3>{
       _Tp3
      operator()(const _Tp1& __x, const _Tp2& __y) const
       { return __x * __y; }
}; 
*/
// generalization of std::multiplies<T>
namespace stdx{
template <class T1, class T2, class R>
struct multiplies : public boost::function_traits<R(T1 const&,T2 const&)> {//public std::binary_function<_Tp1, _Tp2, _TpRes>{
	typedef boost::function_traits<R(T1 const&,T2 const&)> base_; 
	typedef typename base_::result_type result_type;
	typedef typename base_::first_argument_type first_argument_type;
	typedef typename base_::second_argument_type second_argument_type;
	result_type operator()(first_argument_type t1, second_argument_type t2) const{
		return t1 * t2;
	}
};
}

template<class T> struct elementwise_scalar{typename T::value_type typedef type;};

template<typename T, unsigned D>
struct elementwise_scalar<boost::multi_array<T, D> >{typedef T type;};

template<
  class Element,
  class Commutative_Associative_Addition=elementwise_addition<Element>,
  class Field_Scalar=typename elementwise_scalar<Element>::type,
  class Associative_ScalarMultiplication=elementwise_multiplication<Field_Scalar, Element> 
>
class linear_space : public Element,
  boost::additive      <linear_space<Element, Commutative_Associative_Addition, Field_Scalar, Associative_ScalarMultiplication>,//EBCO
  boost::multiplicative<linear_space<Element, Commutative_Associative_Addition, Field_Scalar, Associative_ScalarMultiplication>, Field_Scalar >//EBCO
  >
{
  public:
  typedef Element                          element;
  typedef Field_Scalar                     scalar;
  typedef Commutative_Associative_Addition addition;
  typedef Associative_ScalarMultiplication scalar_multiplication;

	static addition addition_instance;
	static scalar_multiplication scalar_multiplication_instance;

  linear_space(element const& e) : element(e){}
  friend linear_space& operator+=(linear_space& self, linear_space const& other){	//static const addition add=addition();
		self = addition_instance(self, other);			//addition::apply(self, other ); 
		return self;
  }
  friend linear_space& operator*=(linear_space& self, scalar const& s){
		(element&)self=scalar_multiplication_instance(s,(element&)self); //scalar_multiplication::apply(s,(element&)self);
		return self;
  }
  friend linear_space& operator-=(linear_space& self, linear_space const& other){self += (-scalar(1))*other ; return self;}
  friend linear_space& operator/=(linear_space& self, scalar const& s){self *= (scalar(1)/s)       ; return self;}

  friend linear_space  operator-(linear_space const& v){linear_space ret(v); ret*=-scalar(1); return ret;} 
  friend linear_space const& operator+(linear_space const& v){return v;                              }
};

//static definitions
template<
  class Element,
  class Commutative_Associative_Addition,
  class Field_Scalar,
  class Associative_ScalarMultiplication 
> typename linear_space<Element, Commutative_Associative_Addition, Field_Scalar, Associative_ScalarMultiplication>::addition
linear_space<Element, Commutative_Associative_Addition, Field_Scalar, Associative_ScalarMultiplication>::addition_instance;

template<
  class Element,
  class Commutative_Associative_Addition,
  class Field_Scalar,
  class Associative_ScalarMultiplication 
> typename linear_space<Element, Commutative_Associative_Addition, Field_Scalar, Associative_ScalarMultiplication>::scalar_multiplication
linear_space<Element, Commutative_Associative_Addition, Field_Scalar, Associative_ScalarMultiplication>::scalar_multiplication_instance;


//templace<class LinearSpaceElement>
//class one_dimensional_subspace : protected LinearSpaceElement{
//	one_dimensional_subspace(LinearSpaceElement const& e) : LinearSpaceElement(e){}
//};

template<
	class LinearSpaceElement,
	class InnerProduct=elementwise_inner_product<typename LinearSpaceElement::element> >
class inner_product_space : public LinearSpaceElement{
	public:
	typedef	LinearSpaceElement linear_space;
	typedef InnerProduct inner_product;
	inner_product_space(linear_space const& e) : LinearSpaceElement(e){}
	friend typename linear_space::scalar operator,(inner_product_space const& self, inner_product_space const& other){
		return inner_product::apply(
			self, 
			other
		);
	}
	friend typename linear_space::scalar inner_product(inner_product_space const& self, inner_product_space const& other){
		return inner_product::apply(self, other);
	}
};

#ifdef _TEST_LINEAR_SPACE_HPP

#include<boost/array.hpp>
#include "alf/boost/array_io.hpp"
namespace cartesian{
typedef linear_space<boost::array<double, 3> > vector;
}

typedef linear_space<boost::multi_array<std::complex<double>, 2> > gl;
typedef inner_product_space<linear_space<boost::multi_array<std::complex<double>, 2> >  > glp;
#include <iostream>
//#include <pstade/oven/initial_values.hpp>
#include<boost/smart_ptr.hpp>
#include<alf/mathematica/archive.hpp>
using std::cout; using std::clog; using std::endl;
using namespace boost;

//using cartesian::vector;

int main(){
	typedef array<double, 3> double3;
	clog<<sizeof(array<double, 3>)<<endl;
	clog<<sizeof(cartesian::vector)<<endl;
	cartesian::vector v((double3){{1,0,1}});
	clog<<v/2.<<endl;
	inner_product_space<cartesian::vector> vp(((double3){{0,1,0}}));
	inner_product_space<cartesian::vector> vp2(v);
	clog<< (vp , vp2) <<endl;
	
	
	boost::multi_array<std::complex<double>, 2> a(boost::extents[2][2]);
	boost::multi_array<std::complex<double>, 2> b(boost::extents[2][2]);
	for(int i=0; i!=2; ++i)
		for(int j=0; j!=2; ++j){
			a[i][j]=i; b[i][j]=std::complex<double>(j,1);
		}
	gl gl1(a);
	gl gl2(b);
	gl gl3 = gl1 + gl2;
	mathematica::oarchive oa(cout);
	oa<<gl3;
	auto d=3;
	glp glp1(gl1);
	glp glp2(gl2);
	std::complex<double> c=inner_product(glp1,glp1);
	glp1/=sqrt(c);
	cout<<(glp1, glp1)<<endl;
    return 0;
}
#endif

