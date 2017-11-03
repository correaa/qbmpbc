#ifdef compile_instructions
ln -sf $0 $0.cpp && c++ -Wfatal-errors $0.cpp -L$HOME/lib `pkg-config --libs gsl` -D_TEST_QUASI_RANDOM_HPP -o ./$0.x && ./$0.x $@ ; exit
#endif
#ifndef QUASI_RANDOM_HPP
#define QUASI_RANDOM_HPP

#include<gsl/gsl_qrng.h>
#include<boost/array.hpp>
namespace gsl{
	namespace quasi_random{
		template<unsigned Dim>
		class generator{
			gsl_qrng* pimpl_;
			public:
			typedef gsl_qrng_type const* type;
			static type const niederreiter_2; //it is valid up to 12 dimensions. 
			static type const sobol; //It is valid up to 40 dimensions. 
			static type const halton; //up to 1229
			static type const reverse_halton;
			generator(type const& t=sobol) : pimpl_(gsl_qrng_alloc(t, Dim) ){}
			generator(generator const& other) : pimpl_(gsl_qrng_clone(other.pimpl_)){} 
			boost::array<double, Dim> operator()() const{
				boost::array<double, Dim> ret;
				gsl_qrng_get(pimpl_, ret.c_array());
				return ret;
			}
			generator& operator=(generator const& other){
				if(this!=&other){
					//assert(this->name()==other.name());
					gsl_qrng_memcpy(this->pimpl_, other.pimp_);
				}
				return *this;
			}
			std::string name() const{return gsl_qrng_name(pimpl_);}
			size_t size() const{return gsl_qrng_size(pimpl_)/* = Dim*/;}
			//<ret> state(){ void * gsl_qrng_state (const gsl_qrng * q) }
			~generator(){gsl_qrng_free(pimpl_);}
		};
		template<unsigned Dim> typename generator<Dim>::type const generator<Dim>::niederreiter_2 = gsl_qrng_niederreiter_2;
		template<unsigned Dim> typename generator<Dim>::type const generator<Dim>::sobol          = gsl_qrng_sobol;
		template<unsigned Dim> typename generator<Dim>::type const generator<Dim>::halton         = gsl_qrng_halton;
		template<unsigned Dim> typename generator<Dim>::type const generator<Dim>::reverse_halton = gsl_qrng_reversehalton;
	}
}
#ifdef _TEST_QUASI_RANDOM_HPP
#include<iostream>
//#include<boost/fusion.hpp>
#include <boost/fusion/sequence.hpp>
//#include <boost/fusion/include/sequence.hpp>
#include<boost/fusion/adapted/array.hpp>
int main(){
	gsl::quasi_random::generator<1> g;             //(gsl::quasi_random::generator<1>::sobol);
	for(unsigned i=0; i<100; ++i){
		boost::array<double, 1> a=g();
		std::cout<<a[0]<<std::endl;
	}	
	return 0;
}
#endif
#endif

