#if compilation_instructions
ln -sf $0 .$0.cpp && c++ .$0.cpp `#-Wfatal-errors` -L$HOME/lib `pkg-config --libs gsl` -D_TEST_GSL_RANDOM_HPP -o ./.$0.x && ./.$0.x $@; exit
#endif
#ifndef GSL_RANDOM_HPP
#include<gsl/gsl_rng.h>
#include<string>
namespace gsl{
namespace random{
class generator{
	public:
	gsl_rng* pimpl_;
	typedef gsl_rng_type const* type;
	static type const mt19937;   //< Mersenne Twister
	static type const ranlxs0;   //< luxury random numbers single precission level 0
	static type const ranlxs1;   //< luxury random numbers single precission level 1
	static type const ranlxs2;   //< luxury random numbers single precission level 2
	static type const ranlxd1;   //< luxury random numbers double precission level 1
	static type const ranlxd2;   //< luxury random numbers double precission level 2
	static type const ranlux;    //< original luxury random numbers
	static type const ranlux389; //< original luxury random numbers with the highest level of randomness
	static type const cmrg;      //< L'Ecuyer
	static type const mrg;       //< fifth-order multiple recursive generator by L'Ecuyer, Blouin and Coutre
	static type const taus;      //< Tausworthe generator
	static type const taus2;     //< Tausworthe generator with improved seeding
	static type const gfsr4;     //< Lagged-fibonacci generator
	generator(type const& t=mt19937) : pimpl_(gsl_rng_alloc(t)){}
	generator(generator const& other) : pimpl_(gsl_rng_clone(other.pimpl_)){}
	generator& operator=(generator const& other){
		if(this==&other){
			//assert(this->name()==other.name());
			gsl_rng_memcpy(pimpl_, other.pimpl_);
		}
		return *this;
	}
	void set(unsigned long seed) const{gsl_rng_set(pimpl_, seed);}
	unsigned long get() const{return gsl_rng_get(pimpl_);}
	double uniform() const{return gsl_rng_uniform(pimpl_);}
	double uniform_pos() const{return gsl_rng_uniform_pos(pimpl_);}
	double positive_uniform() const{return uniform_pos();}
	unsigned long uniform_int(unsigned long int n){return gsl_rng_uniform_int(pimpl_, n);}
	~generator(){gsl_rng_free(pimpl_);}
	/// http://www.gnu.org/software/gsl/manual/html_node/Auxiliary-random-number-generator-functions.html
	std::string name() const{return gsl_rng_name(pimpl_);}
	unsigned long max() const{return gsl_rng_max(pimpl_);}
	unsigned long min() const{return gsl_rng_min(pimpl_);}
	static std::basic_string<type> types(){return gsl_rng_types_setup();}
	/// http://www.gnu.org/software/gsl/manual/html_node/Random-number-generator-algorithms.html
};
	generator::type const generator::mt19937   = gsl_rng_mt19937; //< Mersenne Twister
	generator::type const generator::ranlxs0   = gsl_rng_ranlxs0; //< luxury random numbers single precission level 0
	generator::type const generator::ranlxs1   = gsl_rng_ranlxs1; //< luxury random numbers single precission level 1
	generator::type const generator::ranlxs2   = gsl_rng_ranlxs2; //< luxury random numbers single precission level 2
	generator::type const generator::ranlxd1   = gsl_rng_ranlxd1; //< luxury random numbers double precission level 1
	generator::type const generator::ranlxd2   = gsl_rng_ranlxd2; //< luxury random numbers double precission level 2
	generator::type const generator::ranlux    = gsl_rng_ranlux; //< original luxury random numbers
	generator::type const generator::ranlux389 = gsl_rng_ranlux389; //<  original luxury random numbers with the highest level of randomness
	generator::type const generator::cmrg      = gsl_rng_cmrg; //< L'Ecuyer
	generator::type const generator::mrg       = gsl_rng_mrg; //< fifth-order multiple recursive generator by L'Ecuyer, Blouin and Coutre
	generator::type const generator::taus      = gsl_rng_taus; //< Tausworthe generator
	generator::type const generator::taus2     = gsl_rng_taus2; //< Tausworthe generator with improved seeding
	generator::type const generator::gfsr4     = gsl_rng_gfsr4; //<  Lagged-fibonacci generator
	// unix generators not included http://www.gnu.org/software/gsl/manual/html_node/Unix-random-number-generators.html
	// low quality generators not included http://www.gnu.org/software/gsl/manual/html_node/Other-random-number-generators.html
}
}

#include <gsl/gsl_randist.h>
namespace gsl{
namespace random{
	class gaussian{
		double sigma_;
		double mu_;
		generator g_;
		public:
		gaussian(double sigma=1, double mu=0) : sigma_(sigma), mu_(mu){}
		double operator()() const{return mu_ + gsl_ran_gaussian(g_.pimpl_, sigma_);}
		double pdf(double x) const{return gsl_ran_gaussian_pdf(x, sigma_);} //< probability density p
	};
}
}
#endif //GSL_RANDOM_HPP
#ifdef _TEST_GSL_RANDOM_HPP
#include<iostream>
using namespace std;
int main(){
	using namespace gsl::random;
	std::basic_string<generator::type> t = generator::types();
	cout<< "there are "<<t.size()<<" generators available"<<endl;
	for(std::basic_string<generator::type>::const_iterator it=t.begin(); it!=t.end(); ++it){
		cout << (*it)->name <<'\n';
	}
	generator g;
	cout << "one random number between [0, 1) " << g.uniform() << endl;
	cout << "one random number between [0, 1) " << g.uniform() << endl;
	cout << "one random number between [0, 1) " << g.uniform() << endl;

	gaussian G;
	cout << "one gaussian "<< G() << endl;
	cout << "one gaussian "<< G() << endl;
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

