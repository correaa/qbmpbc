#if 0
	cp $0 $0.cpp && c++ $0.cpp -L$HOME/lib -lgsl -lcblas -latlas -D_TEST_MONTE_CARLO_HPP -o ./$0.cpp.x && ./$0.cpp.x $1 $2 $3 $4 $5 $6 $7 $8 $9
	rm -f $0.cpp.x $0.cpp
	exit
#endif
#ifndef MONTE_CARLO_HPP
#define MONTE_CARLO_HPP

#include <cstdlib>
//#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include<boost/function.hpp>
#include<boost/array.hpp>

#include<boost/numeric/interval.hpp>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp> 
#include <boost/accumulators/statistics/weighted_variance.hpp>

#include <boost/operators.hpp> //dereferenceable
#include <limits>
#include <iostream>
/// Wrapper to GSL library http://www.gnu.org/software/gsl/
namespace gsl{
/// Monte Carlo (scalar) integration in arbitrary dimensions http://www.gnu.org/software/gsl/manual/html_node/Monte-Carlo-Integration.html
namespace monte_carlo{
	using std::clog; using std::endl; 
	using boost::array;
	namespace bn = boost::numeric;
	typedef bn::interval<double> interval; //< Defines an 1D integration interval
	
	/// Stores the result (estimate and error) of a Monte Carlo integration
	struct result : std::pair<double,double>{
		double const& estimate() const{return first;}
		double const& error() const{return second;}
		bn::interval<double> interval() const{return bn::interval<double>(this->estimate()-error(),estimate()+error());} //make it operator?
		result(std::pair<double, double> const& p) : std::pair<double, double>(p){}
		friend std::ostream& operator<<(std::ostream& os, result const& r){return os<<r.estimate()<<" \u00B1 "<<r.error();} // \ub00B1 is unicode +/-
	};
	/// Function inteface for integration, it captures a suitable function object, lambda expression, or free function
	template<size_t NDim>
	struct function : gsl_monte_function, boost::function<double(array<double, NDim> const&)>{
		private:
		void init(){
			gsl_monte_function::f = &(function<NDim>::f_free_);
			gsl_monte_function::dim = NDim;
			gsl_monte_function::params = this;
		}
		public:
		template<class T>
		function(T f_object) : boost::function<double(array<double, NDim> const&)>(f_object){init();}
		function(function const& f) : 
			gsl_monte_function((gsl_monte_function const&)f), 
			boost::function<
				double(array<double, NDim> const&)
			>((boost::function<double(array<double, NDim> const&)> const&)f){
			init();
		}
		function(
			boost::function<double(array<double, NDim> const&)> const& f_object
		) : boost::function<double(array<double, NDim> const&)>(f_object){
			init();
		}
		private:
		static double f_free_(double x[], size_t, void* self){
			return ((function<NDim>*)self)->operator()(*(boost::array<double, NDim>*)x);
		}
	};
	/// Limits is a set of intervals that defines a multidimensional rectangle (cartesian product) for integration
	template<size_t NDim> 
	struct limits : array<bn::interval<double>, NDim>{
		typedef array<bn::interval<double>, NDim> base;
		limits(base const& arr) : base(arr){}
	};
	template<size_t NDim> class integral; //< Represents an indefinite integral, when associated with a domain (limits) it can be numerically evaluated
	template<size_t NDim> class integrator; //< Integrator represents and integration method with its concrete options
	using namespace boost::accumulators;
	// why <3>?
	class iterator : public boost::dereferenceable<iterator, result const*>{ //< Lazy forward iterator that will dereference succesive answers
		accumulator_set<double, stats<tag::weighted_mean(immediate), tag::sum_of_weights>, double> acc_mean;
		integrator<3> const* ir_; // change to smart pointer?
		integral<3> const* id_; // change to smart pointer?
		iterator(integrator<3> const& ir, integral<3> const& id) : ir_(&ir), id_(&id){}
		friend class integrator<3>;
		public:
		iterator& operator+=(size_t iters);
		result operator*() const{
			return result(std::pair<double,double>(
				boost::accumulators::weighted_mean(acc_mean),
				std::sqrt(1./boost::accumulators::sum_of_weights(acc_mean))
			));
		}
	};
	/// Abstract Monte Carlo integrator, an integrator has an associated method, its paramaters (if any) and a random number generator.
	template<size_t NDim>
	class integrator{
		protected:
		gsl_rng* rP_;
		integrator(){ // if user were able to choose the random source it would be implmented here
			const gsl_rng_type *T;
			gsl_rng_env_setup();
			T = gsl_rng_default;
			rP_ = gsl_rng_alloc(T);
		}
		virtual result operator()(integral<NDim> const&, size_t chunk) const=0;
		friend class iterator;
		public:
		iterator operator()(integral<NDim> const& il) const{return iterator(*this, il);}
		virtual ~integrator(){gsl_rng_free(rP_);}
	};
	/// The integrator performs a number of evaluations of the function to improve the result estimation.
	iterator& iterator::operator+=(size_t iters){
		gsl::monte_carlo::result partial = ir_->operator()(*id_, iters);
		acc_mean(partial.estimate(), weight=1./(partial.error()*partial.error()));
		return *this;
	}
	/// Represents a definite integral in a concrete rectangular domain, and integrator applied to it can provide estimates
	template<size_t NDim>
	class integral : std::pair<function<NDim>, limits<NDim> >{
		typedef std::pair<function<NDim>, limits<NDim> > base;
		public:
		integral(function<NDim> const& f, limits<NDim> const& l) : base(f,l){}
		function<NDim> const& integrand() const{return this->first;}
		limits<NDim> const& domain() const{return this->second;}
		iterator operator|(integrator<NDim> const& ir) const{return ir(*this);}
		//iterator operator|() const{return vegas::integrator<NDim>()(*this);}
	};
	/// Plain Monte Carlo integration http://www.gnu.org/software/gsl/manual/html_node/PLAIN-Monte-Carlo.html
	namespace plain{
		struct workspace{
			gsl_monte_plain_state* sp_;
			workspace(size_t dim) : sp_(gsl_monte_plain_alloc(dim)){reset();}
			void reset(){gsl_monte_plain_init(sp_);}
			~workspace(){gsl_monte_plain_free(sp_);}
		};
		/// Plain Monte Carlo integrator, no options (except for static dimensionality)
		template<size_t NDim>
		class integrator : workspace, public monte_carlo::integrator<NDim>{
			public:
			integrator() : workspace(NDim){}
			friend class iterator;
			using monte_carlo::integrator<NDim>::operator();
			protected:
			monte_carlo::result operator()(integral<NDim> const& i, size_t calls) const{
				double result_;
				double abserr_;
				{
					array<double, NDim> xl_;
					array<double, NDim> xu_;
					for(size_t j=0; j!=NDim; ++j){
						xl_[j] = lower(i.domain()[j]);
						xu_[j] = upper(i.domain()[j]);
					}
					gsl_monte_plain_integrate(
						const_cast<gsl_monte_function*>((const gsl_monte_function*)&i.integrand()),       //gsl_monte_function * f, 
						&xl_[0],   //const double xl[],
						&xu_[0],   //const double xu[], 
						NDim,      //size_t dim, 
						calls,     //size_t calls, 
						this->rP_, //gsl_rng * r,
						this->sp_, //gsl_monte_miser_state * s,
						&result_,  //double * result, 
						&abserr_   //double * abserr)
					);
				}
				return result(std::pair<double, double>(result_, abserr_) );
			}
		};
	}
	/// Miser Monte Carlo integration http://www.gnu.org/software/gsl/manual/html_node/MISER.html
	namespace miser{
		/// Miser parameters (see http://www.gnu.org/software/gsl/manual/html_node/MISER.html)
		/** Miser parameters can be chained parameters().estimate_fraction(0.15).variance_scaling().estimate_min_calls_per_bisection(100). It is initialized with default (by GSL manual) values, to restore default (or set recommend) value call the option with no value. */
		template<size_t NDim>
		struct parameters : gsl_monte_miser_params{
			parameters(){
				//default values go here
				estimate_frac = 0.1;
				min_calls = 16 * NDim;
				min_calls_per_bisection = 32 * min_calls;
				alpha = 2.;
				dither = 0.;
			}
			//recommended values go here if different from default values
			parameters& estimate_fraction(double f=0.10){estimate_frac=f; return *this;}
			parameters& estimate_min_calls(size_t mc = 16*NDim){min_calls=mc; return *this;}
			parameters& estimate_min_calls_per_bisection(size_t mc){min_calls_per_bisection=mc; return *this;}
			parameters& estimate_min_calls_per_bisection(){min_calls_per_bisection=32*min_calls; return *this;}
			parameters& variance_scaling(double a = 2.){alpha=a; return *this;}
			parameters& bisection_dither(double d=0.1){dither=d; return *this;}
		};
		struct workspace{
			gsl_monte_miser_state* sp_;
			workspace(size_t dim) : sp_(gsl_monte_miser_alloc(dim)){reset();}
			void reset(){gsl_monte_miser_init(sp_);}
			~workspace(){gsl_monte_miser_free(sp_);}
		};
		template<size_t NDim>
		class integrator : workspace, public monte_carlo::integrator<NDim>{
			public:
			integrator(parameters<NDim> const& p = parameters<NDim>()) : workspace(NDim){
				gsl_monte_miser_params_set(this->sp_, &p);
			}
			friend class iterator;
			protected:
			monte_carlo::result operator()(integral<NDim> const& i, size_t calls) const{
				double result_;
				double abserr_;
				{
					//monte_carlo::function<NDim> f_(i.integrand());
					array<double, NDim> xl_;
					array<double, NDim> xu_;
					for(size_t j=0; j!=NDim; ++j){
						xl_[j] = lower(i.domain()[j]);
						xu_[j] = upper(i.domain()[j]);
					}
					gsl_monte_miser_integrate(
						const_cast<gsl_monte_function*>((const gsl_monte_function*)&i.integrand()),       //gsl_monte_function * f, 
						&xl_[0],   //const double xl[],
						&xu_[0],   //const double xu[], 
						NDim,      //size_t dim, 
						calls,     //size_t calls, 
						this->rP_,    //gsl_rng * r,
						this->sp_,    //gsl_monte_miser_state * s,
						&result_,  //double * result, 
						&abserr_   //double * abserr)
					);
				}
				return result(std::pair<double, double>(result_, abserr_) );
			}
		};
	}
	/// Vegas Monte Carlo integration http://www.gnu.org/software/gsl/manual/html_node/VEGAS.html
	namespace vegas{
		/// Vegas parameters (see http://www.gnu.org/software/gsl/manual/html_node/VEGAS.html)
		/** Vegas parameters can be chained parameters().rebinning_stifness(1.5).iteratotions_per_call().estimate_min_calls_per_bisection(100). It is initialized with default (by GSL manual) values, to restore default (or set recommend) value call the option with no value. */
		template<size_t NDim>
		struct parameters : gsl_monte_vegas_params{
			parameters(){
				//default values go here
				alpha = 1.5;
				iterations = 5u;
				stage = 0;
				mode = GSL_VEGAS_MODE_IMPORTANCE;
				verbose = -1;
				ostream = NULL;
			}
			parameters& rebinning_stifness(double a=1.5){alpha=a; return *this;}
			parameters& iterations_per_call(size_t ipc=5){iterations=ipc; return *this;}
		};
		struct workspace{
			gsl_monte_vegas_state* sp_;
			workspace(size_t dim) : sp_(gsl_monte_vegas_alloc(dim)){reset();}
			void reset(){gsl_monte_vegas_init(sp_);}
			~workspace(){gsl_monte_vegas_free(sp_);}
		};
		template<size_t NDim>
		class integrator : workspace, public monte_carlo::integrator<NDim>{
			public:
			integrator() : workspace(NDim){}
			integrator(parameters<NDim> const& p /*=parameters<NDim>()*/) : workspace(NDim){
				gsl_monte_vegas_params_set(this->sp_, &p);
			}
			friend class iterator;
			using monte_carlo::integrator<NDim>::operator();
			protected:
			monte_carlo::result operator()(integral<NDim> const& i, size_t calls) const{
				double result_;
				double abserr_;
				{
					array<double, NDim> xl_;
					array<double, NDim> xu_;
					for(size_t j=0; j!=NDim; ++j){
						xl_[j] = lower(i.domain()[j]);
						xu_[j] = upper(i.domain()[j]);
					}
					gsl_monte_vegas_integrate(
						const_cast<gsl_monte_function*>((const gsl_monte_function*)&i.integrand()),       //gsl_monte_function * f, 
						&xl_[0],   //const double xl[],
						&xu_[0],   //const double xu[], 
						NDim,      //size_t dim, 
						calls,     //size_t calls, 
						this->rP_,    //gsl_rng * r,
						this->sp_,    //gsl_monte_miser_state * s,
						&result_,  //double * result, 
						&abserr_   //double * abserr)
					);
				}
				return result(std::pair<double, double>(result_, abserr_) );
			}
		};
	}
}}
#endif

#ifdef _TEST_MONTE_CARLO_HPP
#include<boost/numeric/interval/io.hpp> //for cout<<interval
#include<list>
double g(boost::array<double,3> const& k){
 //double A = 1.0 / (M_PI * M_PI * M_PI);
 //return A / (1.0 - cos (k[0]) * cos (k[1]) * cos (k[2]));
 return k[0]*k[1]*k[2];
}
double gnc(boost::array<double,3>& k){
 return k[0]*k[1]*k[2];
}
struct gp{
	double scale;
	gp(double s=1) : scale(s){}
	double operator()(boost::array<double,3> const& k) const{
	 return k[0]*k[1]*k[2];
	}
};
int main(){
	using boost::numeric::interval;
	boost::array<interval<double>,3> d={{
		interval<double>(0.,M_PI), 
		interval<double>(0.,M_PI), 
		interval<double>(0.,M_PI)
	}};
	/*
	{ // original gsl interface
		gsl::monte_carlo::integrator<3> i;
		//boost::function<double(boost::array<double,3> const&)> gf = gp(1.);
	  gsl::monte_carlo::result result = i.integrate(gp(1.), d, 50000);
		std::cout<<"simple total "<<result.estimate()<<" +/- "<<result.error()<<std::endl;
	}
	{ //smart interface
		gsl::monte_carlo::result smart = gsl::monte_carlo::integral<3>(&g, d, 0.626,0,50000);
		std::cout<<"smart total "<<smart.estimate()<<" +/- "<<smart.error()<<std::endl;
	}*/
	boost::function<double(boost::array<double,3> const&)> bf(&g);
	gsl::monte_carlo::function<3> mcf(bf);
//	double mfcat = mcf((boost::array<double,3>){{1.,2.,3.}}); clog<<mfcat<<endl;
	using namespace gsl::monte_carlo;
//	gsl::monte_carlo::integral<3> ii(mcf, gsl::monte_carlo::limits<3>(d) );
	//miser::integrator<3> mir;
//	boost::numeric::interval<double> inter = miser::integrator<3>(miser::parameters<3>().bisection_dither(0.1))(integral<3>(&g, d) , (size_t)500000);
	//boost::numeric::interval<double> inter = vegas::integrator<3>()(integral<3>(&g, d) , (size_t)500000);

	//gsl::monte_carlo::result i = (integral<3>(&g, d)|vegas::integrator<3>())[1000];
	boost::numeric::interval<double> inter;
	clog
		<<*(integral<3>(gp(),d)|plain::integrator<3>()+=50000)<<endl
		<<*(integral<3>(gp(),d)|miser::integrator<3>()+=50000)<<endl
		<<*(integral<3>(gp(),d)|vegas::integrator<3>()+=50000)<<endl;
/*
	//iterator i = integral<3>(&g,d)|vegas::integrator<3>()+=50000;
	//i+=50000;
	gsl::monte_carlo::result mcr = *i;
	boost::numeric::interval<double> inter = mcr.interval();	
  std::cout<<"inter ="<<inter<<std::endl;
	i+=50000;
	inter = (*i).interval();	
	std::cout<<"inter ="<<inter<<endl;
  std::cout<<"inter ="<<(*i).estimate()<<" "<<(*i).error()<<std::endl;
*/
	std::cout<<"exact 120.174"<<std::endl;
	return 0;
}
#endif
// ~/usr/bin-hera/astyle --brackets=attach --indent=tab --indent-col1-comments --pad-oper --delete-empty-lines --add-brackets --keep-one-line-statements --convert-tabs --align-pointer=type --unpad-paren  qbox.hpp
// Editor modelines  -  http://www.wireshark.org/tools/modelines.html
// Local variables:
// c-basic-offset: 4
// tab-width: 4
// indent-tabs-mode: t
// truncate-lines: 1
// End:
/* vim:set ft=cpp ts=4 sw=4 sts=4 cindent: */
