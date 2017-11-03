#if 0
ln -sf $0 $0.cpp && c++ -std=c++0x `echo -Wfatal-errors` $0.cpp -L$HOME/lib `pkg-config --libs gsl` -D_TEST_ODEIV_HPP -o ./$0.x && ./$0.x $@
rm -f $0.x $0.cpp
exit
#endif
#ifndef ODEIV_HPP
#define ODEIV_HPP

#include <cstdlib>
#include <gsl/gsl_errno.h> //for GSL_SUCCESS
#include <gsl/gsl_odeiv.h>
#include<boost/array.hpp>
#include<boost/function.hpp>
#include<boost/tuple/tuple.hpp> //for tie
#include<boost/numeric/interval.hpp>
#include<boost/operators.hpp> //for dereferenceable
#include "boost/concept/assert.hpp"
#include "boost/concept_check.hpp"
#include<iostream>
using std::clog;
using std::endl;

namespace gsl{
namespace odeiv{
	using std::pair;
//		using boost::array;
	template<size_t Dim>	
	struct array : boost::array<double, Dim>{
		typedef boost::array<double, Dim> type;
	};
	template<size_t Dim1, size_t Dim2=Dim1> 
	struct array2 : boost::array<boost::array<double, Dim2>, Dim1>{
		typedef boost::array<boost::array<double, Dim2>, Dim1> type;
	};
	template<size_t Dim>
	struct function_ : 
		boost::function<typename array<Dim>::type(typename array<Dim>::type const&, double)>{
		template<class T>
		function_(T function_object) : 
			boost::function<typename array<Dim>::type(typename array<Dim>::type const&, double)>(function_object){}
		//private:
		static int free_function_(double t, const double y[], double dydt[], void *self){
			int ret = GSL_SUCCESS;
//				array<Dim>::type* dydtarrp = (array<Dim>::type*)&dydt[0];
			try{
				typename array<Dim>::type* dydtarrp=(array<Dim>*)dydt;
				*dydtarrp = ((function_*)self)->operator()(*(array<Dim>*)y,t);
			}catch(...){
				ret++;
			}
			return ret;
		}
	};
	template<size_t Dim>
	struct jacobian_ : 
		boost::function<std::pair<typename array2<Dim,Dim>::type, typename array<Dim>::type>(typename array<Dim>::type const&, double t)>{
		template<class T>
		jacobian_(T jacobian_object) : 
			boost::function<std::pair<typename array2<Dim,Dim>::type, typename array<Dim>::type>(typename array<Dim>::type const&, double t)>(jacobian_object){}
		//private:
		static int free_jacobian_(double t, const double y[], double * dfdy, double dfdt[], void * self){
			int ret= GSL_SUCCESS;
			try{
				typename array2<Dim,Dim>::type*  dfdyarrp=(typename array2<Dim,Dim>::type*)dfdy;
				typename array<Dim>::type*       dfdtarrp=(typename array<Dim>::type*)dfdt;
				*dfdyarrp = ((jacobian_*)self)->operator()(*(array<Dim>*)y,t).first;
				*dfdtarrp = ((jacobian_*)self)->operator()(*(array<Dim>*)y,t).second;
				//boost::tie(*dfdyarrp, *dfdtarrp) = ((jacobian_*)self)->operator()(*(array<Dim>*)y,t);
			}catch(...){
				++ret;
			}
			return ret;
		}
	};
	template<unsigned > class evolve;

	template<size_t Dim>
	class system : 
		gsl_odeiv_system,
		boost::function<typename array<Dim>::type(typename array<Dim>::type const&, double)>	{
		template<unsigned > friend class evolve;
		void init(){
				gsl_odeiv_system::function  = system<Dim>::free_function_;
				gsl_odeiv_system::jacobian  = system<Dim>::free_jacobian_;
				gsl_odeiv_system::dimension = Dim;
				gsl_odeiv_system::params    = this;
		}
		public:
			boost::function<pair< typename array2<Dim,Dim>::type, typename array<Dim>::type >(typename array<Dim>::type const&,double) > jacobian;
			typedef boost::function<typename array<Dim>::type(typename array<Dim>::type const&, double)> base_;
			template<class F>
			system(
				F f //,typename boost::enable_if_c<(boost::UnaryFunction<F, typename array<Dim>::type, typename array<Dim>::type&>::value), void*>::type dummy=0
			) : 
				base_(f), 
				jacobian(0){
				//BOOST_CONCEPT_ASSERT((boost::UnaryFunction<F, typename array<Dim>::type, typename array<Dim>::type>));
				init();
			}
			template<class T, class J> //any vector function, any jacobian function
			system(T t, J j) : base_(t), jacobian(j){init();}
			static int free_function_(double t, const double y[], double dydt[], void *self){
				int ret = GSL_SUCCESS;
				try{
					typename array<Dim>::type* dydtarrp=(array<Dim>*)dydt;
					*dydtarrp = ((system*)self)->operator()(*(array<Dim>*)y,t);
				}catch(...){
					ret = GSL_FAILURE;
				}
				return ret;
			}
			static int free_jacobian_(double t, const double y[], double * dfdy, double dfdt[], void * self){
				int ret= GSL_SUCCESS;
				try{
					typename array2<Dim,Dim>::type*  dfdyarrp=(typename array2<Dim,Dim>::type*)dfdy;
					typename array<Dim>::type*       dfdtarrp=(typename array<Dim>::type*)dfdt;
					//*dfdyarrp = ((sys*)self)->jacobian(*(array<Dim>*)y,t).first;
					//*dfdtarrp = ((sys*)self)->jacobian(*(array<Dim>*)y,t).second;
					boost::tie(*dfdyarrp, *dfdtarrp) = ((system*)self)->jacobian(*(array<Dim>*)y,t);
				}catch(...){
					++ret;
				}
				return ret;
			}
	};
	template<size_t Dim>
	struct sys{
		typedef class system<Dim> type;
	};
	
///////////////////////////////////

		typedef gsl_odeiv_step_type const* step_type; //dynamic type tag variable
		template<unsigned Dim>
		class step{
			struct io : boost::array<array<Dim>, 3>{
				array<Dim> const& value() const{return (*this)[0];}
				array<Dim> const& error() const{return (*this)[1];}
				array<Dim> const& derivative() const{return (*this)[2];}
				array<Dim>& value()     {return (*this)[0];}
				array<Dim>& error()     {return (*this)[1];}
				array<Dim>& derivative(){return (*this)[2];}
				//bn::interval<double> interval() const{return bn::interval<double>(this->estimate()-error(),estimate()+error());}
				io(std::pair<array<Dim>, array<Dim> > const& p) : std::pair<array<Dim>, array<Dim> >(p){}
			};
	   		template <unsigned UDim> friend class evolve;
			gsl_odeiv_step* const s_;
			public:
			typedef step_type type;
			static type const rk2;
			static type const rk4;
			static type const rkf45;
			static type const rkck;
			static type const rk8pd;
			static type const rk2imp;
			static type const rk4imp;
			static type const bsimp; //uses jac
			static type const gear1; 
			static type const gear2; 
			step(step_type const& t, unsigned const dim=Dim) : s_(gsl_odeiv_step_alloc(t, dim)){}
			void reset(){gsl_odeiv_step_reset(s_);}
			std::string name() const{return gsl_odeiv_step_name (s_);}
			unsigned order() const{return gsl_odeiv_step_order(s_);}
			array<Dim> apply(io& y, array<Dim> const& dydt, double t, double h /*dt*/, system<Dim> const& s) const{
				array<Dim> dydt_out;
				gsl_odeiv_step_apply(s_, t, h, y.value().data(), y.error().data(), dydt.data(), &dydt_out[0], &s)?:throw std::runtime_error("gsl error");
				return dydt_out;
			}
			~step(){gsl_odeiv_step_free(s_);}
		};
		step_type const step_rk2    = gsl_odeiv_step_rk2   ; //dynamic type tags
		step_type const step_rk4    = gsl_odeiv_step_rk4   ; //dynamic type tags
		step_type const step_rkf45  = gsl_odeiv_step_rkf45 ; //dynamic type tags
		step_type const step_rkck   = gsl_odeiv_step_rkck  ; //dynamic type tags
		step_type const step_rk8pd  = gsl_odeiv_step_rk8pd ; //dynamic type tags
		step_type const step_rk2imp = gsl_odeiv_step_rk2imp; //dynamic type tags
		step_type const step_rk4imp = gsl_odeiv_step_rk4imp; //dynamic type tags
		step_type const step_bsimp  = gsl_odeiv_step_bsimp ; //dynamic type tags
		step_type const step_gear1  = gsl_odeiv_step_gear1 ; //dynamic type tags
		step_type const step_gear2  = gsl_odeiv_step_gear2 ; //dynamic type tags

		template<unsigned Dim> typename step<Dim>::type const step<Dim>::rk2     = step_rk2   ;
		template<unsigned Dim> typename step<Dim>::type const step<Dim>::rk4     = step_rk4   ;
		template<unsigned Dim> typename step<Dim>::type const step<Dim>::rkf45   = step_rkf45 ; 
		template<unsigned Dim> typename step<Dim>::type const step<Dim>::rkck    = step_rkck  ;
		template<unsigned Dim> typename step<Dim>::type const step<Dim>::rk8pd   = step_rk8pd ; //example case: same evals as rkf45 but less points
		template<unsigned Dim> typename step<Dim>::type const step<Dim>::rk2imp  = step_rk2imp;
		template<unsigned Dim> typename step<Dim>::type const step<Dim>::rk4imp  = step_rk4imp;
 	 	template<unsigned Dim> typename step<Dim>::type const step<Dim>::bsimp   = step_bsimp ;
		template<unsigned Dim> typename step<Dim>::type const step<Dim>::gear1   = step_gear1 ;
		template<unsigned Dim> typename step<Dim>::type const step<Dim>::gear2   = step_gear2 ;

		class control{
			struct epsilon{
				double abs_; //accuracy
				double rel_; //precision
				epsilon(double abs, double rel=0) : abs_(abs), rel_(rel){}
			};
   		template <unsigned UDim> friend class evolve;
			gsl_odeiv_control* const c_;
			public:
			control(epsilon const& eps) : c_(gsl_odeiv_control_y_new(eps.abs_, eps.rel_)){}
			~control(){gsl_odeiv_control_free(c_);}
			std::string name() const{return gsl_odeiv_control_name(c_);}
		};
		template<unsigned Dim>
		class evolve{
			gsl_odeiv_evolve* e_;
			public:
			evolve(unsigned const dim=Dim) : e_(gsl_odeiv_evolve_alloc(dim)){}
			void apply(control const& c, step<Dim> const& s, system<Dim> const& dydt, double& t, double const t1, double& h, boost::array<double, Dim>& y){
				gsl_odeiv_evolve_apply(e_, c.c_, s.s_, &dydt, &t, t1, &h, &y[0])?
					throw std::runtime_error("gsl error: was Jacobian provided? bsimp needs jacobian. Otherwise, function was lost in the middle"):0;
				return;
			}
			~evolve(){gsl_odeiv_evolve_free(e_);}
		};
		template<unsigned Dim>	class iterator;
		typedef boost::numeric::interval<double> interval;
		template<unsigned Dim>
		class problem{
			system<Dim> eqns_;
			typename array<Dim>::type initial_conditions_;
			interval domain_;
			template<unsigned>	friend class iterator;
			public:
			problem(
				system<Dim> const& s, 
				typename array<Dim>::type const& init, 
				interval const& dom
			) : eqns_(s), initial_conditions_(init), domain_(dom){}
			problem(
				system<Dim> const& s, 
				typename array<Dim>::type const& init, 
				double t_init=0.
			) : eqns_(s), initial_conditions_(init), domain_(interval(t_init, std::numeric_limits<double>::infinity() )){}
		};
		template<unsigned Dim>
		class iterator : 	public boost::dereferenceable<iterator<Dim>, std::pair<double, typename array<Dim>::type> const*>
			{
			evolve<Dim> e_;
			problem<Dim> p_;
			//step<Dim> s_;
			//double t_;
			//array<Dim>::type y_;
			double h_;
			std::pair<double, typename array<Dim>::type> state_;
			public:
			iterator(problem<Dim> const& p, double h=1E-6) : p_(p), h_(h),
				state_(
					lower(p_.domain_),
					p_.initial_conditions_
				){}
			iterator& operator++(){
				e_.apply(
					gsl::odeiv::control(1E-6), 
					gsl::odeiv::step<Dim>(gsl::odeiv::step<Dim>::rkf45), 
					p_.eqns_, 
					state_.first, 
					upper(p_.domain_), 
					h_, 
					state_.second
				);
				return *this;
			}
			std::pair<double, typename array<Dim>::type>& operator*(){return state_;	}
			std::pair<double, typename array<Dim>::type> const& operator*() const{return state_;	}
		};
	}//ns odeiv
}//ns gsl
#ifdef _TEST_ODEIV_HPP

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include<iostream>
using namespace std;

int func (double t, const double y[], double f[], void *params){
 double mu = *(double *)params;
 f[0] = y[1];
 f[1] = -y[0] - mu*y[1]*(y[0]*y[0] - 1);
 return GSL_SUCCESS;
}

struct fun{
	int count_;
	fun(double mu) : mu_(mu), count_(0){}
	gsl::odeiv::array<2>::type operator()(gsl::odeiv::array<2>::type const& y, double t=0){
	++count_;
	// f[{x_,y_}]:={y,-x-\[Mu] y (x^2-1)}/.\[Mu]->10	
	return (gsl::odeiv::array<2>::type){{
			 y[1],
			-y[0] - mu_*y[1]*(y[0]*y[0]-1)
		}};
	}
	double mu_;
};
gsl::odeiv::array<2>::type fu(gsl::odeiv::array<2>::type const& y, double t=0){
	return (gsl::odeiv::array<2>::type){{
			 y[1],
			-y[0] - 10.*y[1]*(y[0]*y[0]-1)
	}};
}

struct jac{
	jac(double mu) : mu_(mu){}
	std::pair<gsl::odeiv::array2<2,2>::type, gsl::odeiv::array<2>::type > operator()(gsl::odeiv::array<2>::type const& y, double t=0){
		std::pair<gsl::odeiv::array2<2,2>::type, gsl::odeiv::array<2>::type > ret;
		//In[28]:= D[f[{x, y}][[1]], x]
		//Out[28]= 0
		//In[30]:= D[f[{x, y}][[1]], y]
		//Out[30]= 1
		//In[31]:= D[f[{x, y}][[2]], x]
		//Out[31]= -1 - 2 x y \[Mu]
		//In[32]:= D[f[{x, y}][[2]], y]
		//Out[32]= -(-1 + x^2) \[Mu]
		ret.first[0][0]= 0.;
		ret.first[0][1]= 1.;
		ret.first[1][0]= -2.*mu_*y[0]*y[1]-1.;
		ret.first[1][1]= -mu_*(y[0]*y[0]-1.);
		ret.second[0]=0.;
		ret.second[1]=0.;
		return ret;
	}
	double mu_;
};

int main(){
	fun f(10);
	gsl::odeiv::system<2> const mysys(boost::ref(f));
	using namespace gsl::odeiv;
	for(gsl::odeiv::iterator<2> i(problem<2>(mysys, (array<2>::type){{1.,0}}, 0.), 1e-10); 
		i->first<100.; 
		++i
	){
		cout<< i->first <<' '<< i->second[0] <<' '<< i->second[1]<< endl;
	}
	return 0;
}
#endif //TEST

#endif
