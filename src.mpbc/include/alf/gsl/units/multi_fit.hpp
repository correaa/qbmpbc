#ifdef compile_instructions
ln -sf $0 .$0.cpp && c++ -std=c++0x -Wl,-rpath=$HOME/usr/lib$CLUSTER -I$HOME/usr/include -I$HOME/prj .$0.cpp -Wall `#-Wfatal-errors` -L$HOME/usr/lib `pkg-config --libs gsl`-lboost_regex -lboost_system -lboost_filesystem -laspell -D_TEST_GSL_FIT_LINEAR_BOOST_UNITS_HPP -o ./.$0.x && ./.$0.x $@
rm -f .$0.x .$0.cpp
exit
#endif
#ifndef GSL_FIT_NONLINEAR_BOOST_UNITS_HPP
#define GSL_FIT_NONLINEAR_BOOST_UNITS_HPP

#include "../../gsl/multi_fit.hpp"
#include<boost/units/quantity.hpp>
#include<boost/units/systems/si.hpp>
#include<boost/units/systems/si/codata_constants.hpp>
#include<boost/units/systems/si/io.hpp>
#include<boost/fusion/container.hpp>
#include<boost/fusion/include/container.hpp>
#include<boost/tuple/tuple.hpp>
#include<fstream>
#include<boost/fusion/algorithm/transformation/transform.hpp>
#include<boost/fusion/include/transform.hpp>
#include<boost/version.hpp>
#include<boost/fusion/adapted/boost_array.hpp> //needs boost >1.44, commenting will create problems later when array is not found to be adapted
#include<boost/fusion/algorithm/iteration/for_each.hpp>
#include<boost/fusion/include/for_each.hpp>
#include<boost/fusion/view/zip_view.hpp>
#include<boost/fusion/include/zip_view.hpp>
#include"boost/type_traits/remove_const_reference.hpp"
#include <boost/fusion/include/adapt_struct.hpp>
#include "boost/units/systems/atomic/si_conversion.hpp"
namespace gsl{
namespace fit{
namespace units{
using namespace boost::units;
using namespace boost::fusion;

struct set_first_to_second{
	template<typename T>
	void operator()(T t) const{
		deref(boost::fusion::begin(t)) = deref(boost::fusion::next(boost::fusion::begin(t))).value();   
	}
};
template<class FusionSeqOfQuantities>
boost::array<
	double, 
	result_of::size<FusionSeqOfQuantities>::type::value //FusionVectorOfQuantities::size::value
> to_value(FusionSeqOfQuantities& fv){
	typedef boost::array<
		double, 
		result_of::size<FusionSeqOfQuantities>::type::value //FusionVectorOfQuantities::size::value
	> array; 
	array ret;
	typedef boost::fusion::vector<array&, FusionSeqOfQuantities&> zip_type;
	for_each(
		boost::fusion::zip_view<zip_type>(zip_type(ret, fv)), set_first_to_second()
	);
	return ret;
}

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

template<class QNum>
struct divide{
	template <typename SignatureT> struct result;
	template<class Q> struct 
		result<divide<QNum>(Q)> : 
		divide_typeof_helper<
			QNum,
			typename boost::remove_const_reference<Q>::type
		>{}
	;
};

template<class Model>
struct adimensional : Model{
	enum {static_size = result_of::size<typename Model::parameters_type>::type::value}; 	//enum {static_size = Model::parameters_type::size::value}; //enum {static_size = Model::size::value};
	typedef boost::array<double, static_size> array;
	//adimensional(Model const& m) : Model(m){}
	adimensional(array const& a) 
	:	Model(boost::fusion::transform(
			a, 
			//Model(),
			typename Model::parameters_type(), 
			from_value()
		))
	{}
	double operator()(double const& d) const{
		return Model::operator()(Model::domain_type::from_value(d)).value();
	}
	array da(double const& d) const{
		typename Model::parameters_gradient_type v = Model::da(Model::domain_type::from_value(d));
		array ret = to_value(v);
		return ret;
	}
};

// apparently parameters type and result type have to be in a consistent unit system otherwise nasty things happen!!! to fix, not easy to check via templates
template<class UnitsModel> class 
	solver
		: 
		public boost::iterator_facade<
			solver<UnitsModel>, 
			UnitsModel ,//const, 
			boost::forward_traversal_tag, 
			UnitsModel //const /*noref*/ 
		>
	{
	gsl::fit::solver<gsl::fit::units::adimensional<UnitsModel> > impl_;
	public:
	solver(
		std::vector<std::pair<typename UnitsModel::domain_type, typename UnitsModel::result_type> > const& data,
		UnitsModel const& init
		//typename UnitsModel::parameters_type const& parameters
	) : 
		impl_(
				reinterpret_cast<std::vector<std::pair<double, double> > const&>(data),
				to_value(init.parameters)
		)
	{
	}
	//typedef UnitsModel model_type;
	//#if 0
	UnitsModel /*to be as_vector*/ position() const{
		return transform(
				impl_.position(), 
				typename UnitsModel::parameters_type(), 
				from_value()
		);
	}
	private:
	UnitsModel dereference() const{ // iterator_facade implements operator*
		return UnitsModel( 
			transform(
				impl_.position(), 
				typename UnitsModel::parameters_type(), 
				from_value()
			)
		);
	}
	void increment(){impl_.iterate();} //iterator_facade implements operator++
	friend class boost::iterator_core_access;
	//#endif
};

template<class Domain, class Result, class Parameters>
struct model_facade{
	typedef Domain     domain_type;
	typedef Result     result_type;
	typedef Parameters parameters_type; //must be fusion vector
	parameters_type parameters;
	typedef 
		typename boost::fusion::result_of::as_vector<
			typename boost::fusion::result_of::transform<
				parameters_type, 
				typename gsl::fit::units::divide<result_type> 
			>::type
		>::type
	parameters_gradient_type;
	virtual result_type operator()(domain_type const& d) const = 0;
	virtual parameters_gradient_type da(domain_type const& p) const = 0;
};


}}}
#endif

#ifdef _TEST_GSL_FIT_LINEAR_BOOST_UNITS_HPP

#include "alf/latex.hpp"
#include "alf/latex/vector_io.hpp"
#include<boost/assign.hpp>
#include<boost/units/systems/si.hpp>

using namespace boost::units;
using namespace boost::fusion;

struct dulong{
	typedef quantity<atomic::temperature> domain_type;
	typedef quantity<atomic::energy> result_type;
	typedef boost::fusion::vector<
		result_type, 
		quantity<atomic::heat_capacity>
	> parameters_type;
	parameters_type parameters;
	typedef 
		boost::fusion::result_of::as_vector<
			boost::fusion::result_of::transform<
				parameters_type, 
				gsl::fit::units::divide<result_type> 
			>::type
		>::type
	parameters_gradient_type;

	result_type const& a;
	quantity<atomic::heat_capacity> const& b;
	dulong(
		parameters_type const& p
	) : parameters(p), a(boost::fusion::at_c<0>(parameters)), b(boost::fusion::at_c<1>(parameters)){}
	result_type operator()(domain_type const& t) const{
		using namespace boost::units::auto_conversion_operators::to_right;
		return result_type(a + b*t);
	}
	parameters_gradient_type da(domain_type const& t) const{
		return 
			parameters_gradient_type(
				1.,
				boost::fusion::result_of::at_c<parameters_gradient_type, 1>::type(t)
			)
		;
	}
	friend latex::ostream& operator<<(latex::ostream& os, dulong const& ecm){
		using latex::math;
		return os << (
			latex::math("(a + b \\,\\#1)\\& \\quad \\text{with} \\quad a")=ecm.a, math("b")=ecm.b
		);
	}
};


struct exponential_convergence : gsl::fit::units::model_facade<
		quantity<atomic::time>,
		quantity<atomic::energy>,
		boost::fusion::vector<
			quantity<atomic::energy>, 
			quantity<atomic::energy>, 
			quantity<atomic::time>
		>		
	>{
	parameters_type parameters;
	quantity<atomic::energy> const& a     ;
	quantity<atomic::energy> const& b     ;
	quantity<atomic::time>   const& lambda;

	exponential_convergence(
		parameters_type const& p
	) : parameters(p), a(boost::fusion::at_c<0>(parameters)), b(boost::fusion::at_c<1>(parameters)), lambda(boost::fusion::at_c<2>(parameters)){}

	result_type operator()(domain_type const& t) const{
		return a + b*exp(-t/lambda);
	}
	parameters_gradient_type da(domain_type const& t) const{
		return 
			parameters_gradient_type(
				1.,
				exp(-t/lambda),
				b*exp(-t/lambda)*t/pow<2>(lambda)
			)
		;
	}
	friend latex::ostream& operator<<(latex::ostream& os, exponential_convergence const& ecm){
		using latex::math;
		return os << (
			latex::math("(a + b e^{-\\#1/\\lambda})\\& \\quad \\text{with} \\quad a")=ecm.a, math("b")=ecm.b, math("\\lambda")=ecm.lambda
		);
	}
};

//BOOST_FUSION_ADAPT_STRUCT(
//	exp_conv_model,
//	(quantity<atomic::energy>, a     )
//	(quantity<atomic::energy>, b     )
//	(quantity<atomic::time>  , lambda)
//)

/*
namespace boost { namespace fusion { namespace traits {
    template<>
    struct tag_of<exp_conv_model>
    {
        typedef exp_conv_model_tag type;
    };
}}}*/


/*
using namespace boost::units;
using namespace boost::fusion;

struct vinet{
	typedef boost::mpl::size_t<3> size;
	quantity<si::volume       > V0;
	quantity<si::pressure     > B0;
	quantity<si::dimensionless> B1;
	vinet(
		quantity<si::volume> const&   V0, 
		quantity<si::pressure> const& B0, 
		quantity<si::dimensionless> const& B1
	) : V0(V0), B0(B0), B1(B1){}
	struct result{typedef quantity<si::pressure> type;};
	typedef 
		quantity<si::volume> 
		domain_type
	;
	typedef 
		vector<
			quantity<si::volume>, 
			quantity<si::pressure>, 
			quantity<si::dimensionless>
		>
		parameters_type
	;
	typedef 
		result_of::as_vector<
			result_of::transform<
				parameters_type, 
				gsl::fit::units::divide<result::type> 
			>::type
		>::type
		gradient_type
	;
	vinet(parameters_type const& parameters) : 
		V0(at_c<0>(parameters)), 
		B0(at_c<1>(parameters)), 
		B1(at_c<2>(parameters)){}
	quantity<si::pressure> operator()(quantity<si::volume> const& V) const{
		double X = pow<static_rational<1,3> >(V/V0);
		return 3.*B0*(1.-X)*pow<-2>(X)*exp(3./2.*(B1-1.)*(1-X)); // Jeanloz88
	}
	gradient_type da(quantity<si::volume> const& V) const{
		double X = pow<static_rational<1,3> >(V/V0);
		return 
			gradient_type(
				(B0*((1.5 - 1.5*B1)*V + V0*X*(2.-2.5*X+1.5*B1*X)))/(exp(1.5*(-1. + B1)*(-1. + X))*V*V0),
				(3. - 3.*X)/(exp(1.5*(-1. + B1)*(-1. + X))*pow<2>(X)),
				(B0*(4.5*V + V0*(4.5 - 9.*X)*X))/(exp(1.5*(-1. + B1)*(-1. + X))*V)
			)
		;
	}
};
*/

int main(){
	using namespace boost::units;
	typedef 
		std::vector<
			std::pair<
				quantity<atomic::time>, 
				quantity<atomic::energy> 
			> 
		> 
		data_type
	;
	data_type data; data.reserve(10);
	for(unsigned i = 0.; i != data.capacity(); ++i){
		quantity<atomic::time> t = i*1.*atomic::time_unit;
		data.push_back(std::make_pair(t, 1.23*atomic::hartree + 4.56*atomic::hartree*exp(-t/(0.789*atomic::time_unit))));
	}
	/*
	data_type data = boost::assign::list_of<data_type::value_type>
		(1.0*atomic::time_unit,  1. *atomic::hartree)
		(2.1*atomic::time_unit,  4. *atomic::hartree)
		(3.2*atomic::time_unit,  5. *atomic::hartree)
		(4.3*atomic::time_unit,  6.5*atomic::hartree)
		(6.3*atomic::time_unit,  7. *atomic::hartree)
		(7.1*atomic::time_unit,  7.1*atomic::hartree)
		(7.2*atomic::time_unit,  7.1*atomic::hartree)
		(7.3*atomic::time_unit,  7.2*atomic::hartree)
	;*/
	latex::ostream lout("multi_fit.pdf");
	lout << latex::section("Several Ways to Present Data");
	lout << latex::subsection("Plain Text");
	lout << data;
	lout << latex::subsection("Math Mode");
	latex::equation e(""); e << data;
	lout << e;
	lout << latex::subsection("Tabular");
	lout << "\\begin{tabular}{cc}\n";
	lout << "\\hline\\hline";
	lout << data_type::value_type::first_type::unit_type() << "&" << data_type::value_type::second_type::unit_type() << "\\\\\n";
	lout << "\\hline";
	for(data_type::const_iterator it = data.begin(); it != data.end(); ++it){
		lout << latex::math(it->first.value()) << "&" << latex::math(it->second.value()) << "\\\\\n";
	}
	lout << "\\hline\\hline";
	lout << "\\end{tabular}\n";
	lout << latex::newline;

	lout << latex::subsection("Graphically");
	latex::tikz::picture p; p << data;
	lout << p;

	gsl::fit::units::solver<exponential_convergence> s(
		data, 
		exponential_convergence(boost::fusion::make_vector(1.23*atomic::hartree, 4.56*atomic::hartree, 0.989*atomic::time_unit)) 
	);

	latex::pgfplots::axis::units<atomic::time, atomic::energy> ax;
	latex::pgfplots::axis::units<atomic::time, atomic::energy>::coordinates_type coords; coords(data);
	ax << coords;

	latex::pgfplots::axis::units<atomic::time, atomic::energy>::coordinates_type coords_model("no markers"); 
	for(quantity<atomic::time> t= 0; t < 8.*atomic::time_unit; t += 0.1*atomic::time_unit){
		coords_model << latex::pgfplots::axis::units<atomic::time, atomic::energy>::coordinates_type::pair_type(t, (*s)(t) );
	}
	ax << coords_model;


	latex::tikz::picture tkzp; tkzp << ax;
	lout << tkzp;
	lout << latex::newline;

	//ed = boost::fusion::make_vector(7.*atomic::hartree, -5.*atomic::hartree, 4.*atomic::time_unit); 
	for(unsigned i = 0; i != 10; ++i){
		lout << *s << latex::newline;
		++s;
	}
	{
		latex::pgfplots::axis::units<atomic::time, atomic::energy> ax;
		latex::pgfplots::axis::units<atomic::time, atomic::energy>::coordinates_type coords; coords(data);
		ax << coords;

		latex::pgfplots::axis::units<atomic::time, atomic::energy>::coordinates_type coords_model("no markers"); 
		for(quantity<atomic::time> t= 0; t < 8.*atomic::time_unit; t += 0.1*atomic::time_unit){
			coords_model << latex::pgfplots::axis::units<atomic::time, atomic::energy>::coordinates_type::pair_type(t, (*s)(t) );
		}
		ax << coords_model;
		latex::tikz::picture tkzp; tkzp << ax;
		lout << tkzp;
	}

	/*
	gsl::fit::units::solver<vinet> s_u(
		data, 
		vinet::parameters_type(6.*pow<3>(si::meter), 10.*si::pascal, 1.2)
	);*/
	/*
	for(unsigned i = 0; i!=10; ++i){
		++s_u;
	}
	clog << (*s_u)(6.1*pow<3>(si::meter)) << endl;

	std::vector<
		std::pair<double, double> 
	> const& adim_data = (std::vector<std::pair<double, double> > const&)data;

	gsl::fit::solver<gsl::fit::units::adimensional<vinet> > s(adim_data, (boost::array<double, 3>){{6., 10., 1.2}});
	for(unsigned i = 0; i!=10; ++i){
		clog << s.position()[1] << endl; 
		s.iterate();
	}
	clog << "fitting done" <<endl;
	{std::ofstream ofs("fit.dat");{
		for(
			quantity<si::volume> v = 0.5*pow<3>(si::meter); 
			v<5.0*pow<3>(si::meter); 
			v+=0.01*pow<3>(si::meter)){
			ofs<<v << " " << (*s)(v.value()) << endl;
		}
	}}
	boost::fusion::vector< quantity<si::length>, quantity<si::mass> > v;
	clog << at_c<0>(v) << " " << at_c<1>(v) << endl;
*/
	{
		std::vector<
			std::pair<
				quantity<atomic::temperature>, 
				quantity<atomic::energy> 
			> 
		> data; data.reserve(10);
		for(unsigned i = 0.; i != data.capacity(); ++i){
			quantity<si::temperature> t = i*100.*si::kelvin;
			using namespace boost::units::auto_conversion_operators::to_left;
			quantity<atomic::temperature> ta(t);
			quantity<atomic::energy> en(
				1.23*atomic::hartree + 4.56*atomic::k_B*ta
		);
			data.push_back(
				std::pair<
					quantity<atomic::temperature>, 
					quantity<atomic::energy> 
				>(
					t,
					en
			));
		}
		latex::pgfplots::axis::units<atomic::temperature, atomic::energy> ax; 
		ax << data;
		typedef latex::pgfplots::axis::units<atomic::temperature, atomic::energy>::coordinates_type coordinates;
		dulong d(boost::fusion::make_vector(1.23*atomic::hartree, 3.*atomic::boltzman_constant));
		{
			coordinates coords("no markers");
			for(quantity<atomic::temperature> t = 0; t < quantity<atomic::temperature>(1000.*si::kelvin) ; t+= quantity<atomic::temperature>(10.*si::kelvin)){
				coords << coordinates::pair_type(t, d(t) );
			}
			ax << coords;
		}
		gsl::fit::units::solver<dulong> sol(data , d);
		for(unsigned i=0; i<10; ++i) ++sol;
		{
			coordinates coords("no markers");
			for(quantity<atomic::temperature> t = 0; t < quantity<atomic::temperature>(1000.*si::kelvin) ; t+= quantity<atomic::temperature>(10.*si::kelvin)){
				coords << coordinates::pair_type(t, (*sol)(t) );
			}
			ax << coords;
		}
		latex::tikz::picture p; 
		p << ax; //data;
		lout << latex::newline << p;
	}
	return 0;
}
#endif
// vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent:

