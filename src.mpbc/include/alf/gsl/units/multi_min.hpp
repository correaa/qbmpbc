#ifdef compile_instructions
ln -sf $0 .$0.cpp && c++ -std=c++0x -I$HOME/usr/include -I$HOME/prj .$0.cpp -Wall `#-Wfatal-errors` -L$HOME/lib `pkg-config --libs gsl` -D_TEST_GSL_MULTI_MIN_BOOST_UNITS_HPP -o ./.$0.x && ./.$0.x $@
rm -f .$0.x .$0.cpp
exit
#endif
#ifndef GSL_MULTI_MIN_BOOST_UNITS_HPP
#define GSL_MULTI_MIN_BOOST_UNITS_HPP
#include "../../gsl/multi_min.hpp"
#include<boost/units/quantity.hpp>
#include<boost/fusion/algorithm/transformation/transform.hpp>
#include<boost/fusion/adapted/boost_array.hpp>
#include<boost/fusion/algorithm/iteration/for_each.hpp>
#include<boost/fusion/container.hpp>
#include<boost/fusion/view/zip_view.hpp>

namespace gsl{
namespace multi_min{
namespace units{
using namespace boost::units;
using namespace boost::fusion;


struct from_value{ // binary fusion transformation
	template <typename SignatureT> struct result;
	template <typename Q>
	struct result<from_value(double const&, Q const&)>{
		typedef Q type;
	};
	template<typename Q>
	Q operator()(double const& value, Q const&) const{
		return Q::from_value(value);
	}
};

struct set_first_to_second{
	template<typename T>
	void operator()(T t) const{
		deref(boost::fusion::begin(t)) = deref(boost::fusion::next(boost::fusion::begin(t))).value();   
	}
};

template<class FusionVectorOfQuantities>
boost::array<double, 
	FusionVectorOfQuantities::size::value
> to_value(FusionVectorOfQuantities& fv){
	typedef boost::array<
		double, 
		FusionVectorOfQuantities::size::value
	> array; 
	array ret;
	typedef boost::fusion::vector<array&, FusionVectorOfQuantities&> zip_type;
	for_each(
		boost::fusion::zip_view<zip_type>(zip_type(ret, fv)), set_first_to_second()
	);
	return ret;
}

template<class Function>
struct adimensional{
	Function const& f_;
	adimensional(Function const& f) : f_(f){}
	double operator()(boost::array<double, Function::domain_type::size::value> const& x) const{
		return f_(boost::fusion::transform(x, typename Function::domain_type(), from_value())).value();
	}
};

template<class Functor>
class minimizer :
	public boost::iterator_facade<
		minimizer<Functor>, 
		std::pair<typename Functor::domain_type, typename Functor::result_type > const, 
		boost::forward_traversal_tag, 
		std::pair<typename Functor::domain_type, typename Functor::result_type > const /*noref*/ 
	>
{
	adimensional<Functor> af_;
	gsl::multi_min::minimizer<Functor::domain_type::size::value> impl_;
	public:
	minimizer(
		Functor const& f, 
		typename Functor::domain_type const& starting_point, 
		typename Functor::domain_type const& step_size
	) : af_(f), impl_(af_, to_value(starting_point), to_value(step_size) ){
	}
		// thin gsl interface http://www.gnu.org/software/gsl/manual/html_node/Multimin-Iteration.html
		void iterate(){impl_.iterate();}
		typename Functor::domain_type /*const&*/ x() const{ return
			boost::fusion::transform(impl_.x(), typename Functor::domain_type(), from_value());
	 		//impl_.x(); 
		}
		typename Functor::result_type minimum() const{return Functor::result_type::from_value(impl_.minimum());}
		typename Functor::domain_type /*const&*/ gradient() const{
			return 
				boost::fusion::transform(
					impl_.gradient(), 
					typename Functor::domain_type(), 
					from_value()
				)
			;
		}
		// there is no characteristic step_size because all variables can have different dimension
		// double characteristic_step_size() /*size()*/ const{return gsl_multimin_fminimizer_size(pimpl_);} // Minimizer specific characteristic size for the minimizer s, name could be confusing characteristic_size()? step_size()? 

		// iterator interface
		protected:
		std::pair<typename Functor::domain_type, typename Functor::result_type > dereference() const{
				return std::pair<typename Functor::domain_type, typename Functor::result_type >(x(), minimum());
		}
		//bool equal(minimizer const& other);?
		void increment(){this->iterate();}
    friend class boost::iterator_core_access;
	};

}}}
#endif

#ifdef _TEST_GSL_MULTI_MIN_BOOST_UNITS_HPP
#include<boost/assign.hpp>
#include<boost/units/systems/si.hpp>
#include<boost/units/cmath.hpp>
#include"boost/make_array.hpp"
#include <boost/fusion/container/generation/make_vector.hpp>
#include<fstream>

using namespace boost::units;
struct f{
	quantity<si::dimensionless> const& a;
	quantity<si::area> const& b;
	quantity<si::area> operator()(
		quantity<si::length> const& l
	){
		return b + a*pow<2>(l);
	}
};

using namespace boost::fusion;

struct residual_squares{
	typedef std::vector< std::pair<quantity<si::length>, quantity<si::area> > > data_type;
	data_type const& data;
	residual_squares(
		data_type const& data
	) : data(data) {}
	//quantity<si::area> 
	typedef boost::fusion::vector<quantity<si::dimensionless>, quantity<si::area> > domain_type;
	typedef power_typeof_helper<quantity<si::area>, static_rational<2> >::type result_type;
	result_type operator()(domain_type const& x) const{
		quantity<si::dimensionless> const& a = at_c<0>(x);
		quantity<si::area> const& b = at_c<1>(x);
		power_typeof_helper<quantity<si::area>, static_rational<2> >::type ret(0);
		f f_={a, b};
		for(data_type::const_iterator it=data.begin(); it!=data.end(); ++it){
			ret = ret + pow<2>(  it->second - f_(it->first) );
		}
		return ret;
	}
};

struct a_class{
	double p0, p1, p2, p3, p4;
	a_class(double p0, double p1, double p2, double p3, double p4) : 
		p0(p0), p1(p1), p2(p2), p3(p3), p4(p4){}
	double operator()(boost::array<double, 2> const& d) const{
		double x = d[0];
		double y = d[1];
		return 
			p2 * (x - p0) * (x - p0) +
			p3 * (y - p1) * (y - p1) + p4
		;
	}
};

struct orthoromb_area{
	quantity<si::volume> v;
	orthoromb_area(quantity<si::volume> const& v) : v(v){}
	typedef boost::fusion::vector<quantity<si::length>, quantity<si::length> > domain_type;
	typedef quantity<si::area> result_type;
	quantity<si::area> operator()(domain_type const& args) const{
		quantity<si::length> const& a = at_c<0>(args); 
		quantity<si::length> const& b = at_c<1>(args);
		assert(a>0.*si::meter and b>0.*si::meter);
		quantity<si::length> const c = v/a/b;
		return a*b + b*c + a*c;
	}
};

using std::cout; using std::clog;
int main(){
	typedef std::vector<std::pair<quantity<si::length>, quantity<si::area> > > data_type;
	data_type data = 
		boost::assign::map_list_of //<std::pair<quantity<si::length>, quantity<si::area> > >
		(1.*si::meter,   2.2*pow<2>(si::meter))
		(3.*si::meter,  10.1*pow<2>(si::meter))
		(2.*si::meter,   5.2*pow<2>(si::meter))
		(4.*si::meter,  17.1*pow<2>(si::meter))
	;
	std::sort(data.begin(), data.end());
	{
		std::ofstream ofs("multi_min.dat");
		for(data_type::const_iterator it = data.begin(); it!=data.end(); ++it){
			ofs << it->first << " " << it->second << std::endl;
		}
	}
	residual_squares r(data);
	gsl::multi_min::units::adimensional<residual_squares> ar(r);
	a_class a(1.0,2.0, 10.0, 20.0, 30.0);
	gsl::multi_min::minimizer<2> mnr(
		ar, 
		boost::make_array(10.,10.),
		boost::make_array(1.,1.)
	);
	for(unsigned i=0; i!=100; ++i){
		++ mnr;
		cout << at_c<0>(mnr -> first) << " " << at_c<1>(mnr -> first) <<" " << mnr->second << endl;
	}
	{
			boost::fusion::vector<quantity<si::length>, quantity<si::dimensionless> > starting_point(10.*si::meter, 10.);
			boost::fusion::vector<quantity<si::length>, quantity<si::dimensionless> > step_size(1.*si::meter, 1.);

		gsl::multi_min::units::minimizer<
			residual_squares
		> mnr(
			r,
			residual_squares::domain_type(10. , 10.*pow<2>(si::meter)),
			residual_squares::domain_type(10. , 10.*pow<2>(si::meter))
		);
		for(unsigned i=0; i!=100; ++i){
			++ mnr;
			* mnr;
			cout << at_c<0>(mnr -> first) << " " << at_c<1>(mnr -> first) <<" " << mnr->second << endl;
		}
	}

	orthoromb_area oa(1.*pow<3>(si::meter));
	clog << oa(boost::fusion::make_vector(1.*si::meter, 2.*si::meter) ) << endl;
	gsl::multi_min::units::minimizer<orthoromb_area> mnru(
		oa, 
		boost::fusion::make_vector(1.*si::meter, 2.*si::meter), 
		boost::fusion::make_vector(0.1*si::meter, 0.1*si::meter)
	);
	for(unsigned i=0; i!=50; ++i){
		++ mnru;
		cout << at_c<0>(mnru -> first) << " " << at_c<1>(mnru -> first) <<" " << mnru ->second << endl;
	}

	return 0;
}
#endif

