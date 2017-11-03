#if 0
<<<<<<< .mine
ln -sf $0 $0.cpp && c++ -Wall `#-Wfatal-errors` $0.cpp -L$HOME/lib -I$HOME/usr/include -I$HOME/prj `pkg-config --libs gsl` -lboost_regex -lboost_system -lhunspell -D_TEST_INTERPOLATION_HPP -o ./$0.x && ./$0.x $@; 
rm -f $0.cpp 0.x ; exit
=======
ln -sf $0 $0.cpp && c++ -Wall `#-Wfatal-errors` $0.cpp -L$HOME/lib -I$HOME/usr/include -I$HOME/prj `pkg-config --libs gsl` -lboost_regex -lboost_system -D_TEST_INTERPOLATION_HPP -o ./$0.x && ./$0.x $@; 
rm -f $0.cpp 0.x ; exit
>>>>>>> .r262
#endif
#ifndef INTERPOLATION_HPP
#define INTERPOLATION_HPP

#include<gsl/gsl_interp.h> //doesn't keep data
#include<gsl/gsl_spline.h> //does    keep data
#include<map>
#include<vector>
#include<string>
#include<cassert>
#include<boost/numeric/interval.hpp>
#include<iostream>
#include<boost/lambda/lambda.hpp>
#include<boost/lambda/bind.hpp>
#include<boost/spirit/home/phoenix.hpp>
#include<boost/function.hpp>
namespace gsl{
	namespace interpolation{
	typedef boost::numeric::interval<double> interval;
	class data{
		public:
		data(std::map<double, double> const& data) : 
			xa(data.size()), ya(data.size()){
			{unsigned idx=0; 
			for(std::map<double, double>::const_iterator it=data.begin(); it!=data.end(); ++it, ++idx){
				xa[idx]=it->first;
				ya[idx]=it->second;
			}}
		}
		std::pair<double, double> operator[](unsigned i) const{return std::pair<double, double>(xa[i], ya[i]);}
		size_t size() const{return xa.size();}
		protected:
		// since gsl_interp doesn't keep the data I have to keep it here,
		// moreover its has to be stored as linear array
		std::vector<double> xa;
		std::vector<double> ya;
	};
	typedef const gsl_interp_type* type;
	type const linear = gsl_interp_linear;
	type const polynomial = gsl_interp_polynomial;
	type const cspline = gsl_interp_cspline;
	type const cspline_periodic = gsl_interp_cspline_periodic;
	type const akima = gsl_interp_akima;
	type const akima_periodic = gsl_interp_akima_periodic;

	//same as interpolation but natively keeps a copy of the data
	using namespace boost::lambda;
	using namespace boost::tuples;
	/*template<class Derived>
	struct self_binder{
		//typedef spline Derived;
		template<class LambdaExp>
		lambda_functor<lambda_functor_base<
			action<3, 
				function_action<3> 
			>, 
			tuple<
				double (Derived::* const)(const double&)const, 
				const Derived, //this class
				lambda_functor<LambdaExp> const
			> 
		> >
		operator()(lambda_functor<LambdaExp> const& exp) const{
			return bind(
				static_cast<double(Derived::*)(double const&) const>(&Derived::operator()),
				static_cast<Derived const&>(*this),
				exp
			);
		}
	};*/
	class spline : 
		//public self_binder<spline>, 
		public data //,
		//public boost::phoenix::function<gsl::interpolation::spline> // reference to reference problem in gcc 4.1
	{
		//gsl_spline* const pimpl_;
		gsl_interp* const pimpl_; //using interpolation to make metainfo extraction easier
		public:
		spline(spline const& other) : 
			data(other), 
			//boost::phoenix::function<gsl::interpolation::spline>(*this), 
			pimpl_(gsl_interp_alloc(cspline, other.size() ) ) {
			gsl_interp_init(pimpl_, &xa[0],&ya[0],other.size());
		}
		spline(type const& t, std::map<double, double> const& d) : 
			data(d),
			//boost::phoenix::function<gsl::interpolation::spline>(*this), 
			pimpl_(gsl_interp_alloc(t, d.size())){
			assert(d.size()>gsl_interp_min_size(pimpl_));
			gsl_interp_init(pimpl_,&/*idata.*/xa[0], &/*idata.*/ya[0], d.size());
		}
		spline(std::map<double, double> const& d) :
			data(d),
			//boost::phoenix::function<gsl::interpolation::spline>(*this),
			pimpl_(gsl_interp_alloc(cspline, d.size()))
		{
			assert(d.size()>gsl_interp_min_size(pimpl_));
			gsl_interp_init(pimpl_,&/*idata.*/xa[0], &/*idata.*/ya[0], d.size());
		}
		~spline(){gsl_interp_free(pimpl_);}
		std::string name() const{return gsl_interp_name(pimpl_);}

		//template <typename Arg> typename result<Arg>::type;
	   template <typename Arg> struct result{
	   	typedef double type; //can be generalized
	   };
		//boost::phoenix::actor<boost::phoenix::composite<boost::phoenix::detail::function_eval<1>, boost::fusion::vector<boost::phoenix::value<gsl::interpolation::spline>, boost::phoenix::argument<0>, boost::fusion::void_, boost::fusion::void_, boost::fusion::void_, boost::fusion::void_, boost::fusion::void_, boost::fusion::void_, boost::fusion::void_, boost::fusion::void_> > >
		//boost::phoenix::actor<boost::phoenix::composite<boost::phoenix::detail::function_eval<1>, boost::fusion::vector<boost::phoenix::value<const gsl::interpolation::spline>, boost::phoenix::argument<0>, boost::fusion::void_, boost::fusion::void_, boost::fusion::void_, boost::fusion::void_, boost::fusion::void_, boost::fusion::void_, boost::fusion::void_, boost::fusion::void_> > >

		//using boost::phoenix::function<gsl::interpolation::spline>::operator();
		double operator()(double const& x) const{
			return evaluate(x);
		}
		typedef spline this_type;
		double evaluate(double const& x) const{
			//accelerator a();
			return gsl_interp_eval(pimpl_, &xa[0], &ya[0],x, NULL);
		}
		double evaluate_derivative(double const& x) const{
			return gsl_interp_eval_deriv(pimpl_, &xa[0], &ya[0], x, NULL);
		}
		double evaluate_integral(interval const& limits) const{
			return gsl_interp_eval_integ(pimpl_, &xa[0], &ya[0], lower(limits), upper(limits), NULL); 
		}
		interval domain() const{
			return interval(xa.front(), xa.back());
		}
	};
	class accelerator{
		gsl_interp_accel* const pimpl_;
		accelerator() : pimpl_(gsl_interp_accel_alloc()){}
		~accelerator(){gsl_interp_accel_free(pimpl_);}
		void reset(){gsl_interp_accel_reset(pimpl_);}
	};
	}//ns interpolation
}//ns gsl

#ifndef GSL_INTERPOLATION_LATEX
#include "alf/latex.hpp"
latex::pgfplots::axis& operator<<(
	latex::pgfplots::axis& ax, 
	gsl::interpolation::spline const& s
){
	latex::pgfplots::coordinates c_data("only marks");
	for(unsigned idx = 0; idx!=s.size(); ++idx){
		c_data<<s[idx];
	}
	latex::pgfplots::coordinates c_interp("no marks");
	for(double x=lower(s.domain()); 
		x<upper(s.domain()); 
		x+=width(s.domain())/100.){
		c_interp << std::pair<double, double>(x, s(x));
	}
	return ax << c_interp << c_data;
}
#endif

#ifdef _TEST_INTERPOLATION_HPP
#include<iostream>
#include<boost/numeric/interval/io.hpp>
int main(){
	using namespace std;
	std::map<double, double> data;
	data[1]=3;
	data[2]=5;
	data[3]=6;
	data[6]=8;

	gsl::interpolation::spline p(data);
	std::cout<<p.domain()<<std::endl;
	cout<< p(1.5) <<endl;

	gsl::interpolation::spline pcopy(p);
	cout << pcopy(1.5) << endl;
	double dos = 2.;
//	double resp = p(boost::phoenix::arg_names::arg1)(dos);
//	double res =  (p(boost::lambda::_1+1.)+boost::lambda::_1)(1.5);
//	cout << res << " " << p(1.5+1.) + 1.5 << endl;
	
	{
		using namespace boost::phoenix; using namespace arg_names;
		function<gsl::interpolation::spline> php(p);
		cout << "cacac " << (php(arg1)*arg1)(dos) << endl; 
	   cout << "sss " << endl;
	}
	return 0;

}
#endif

#endif
/* vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent: */

