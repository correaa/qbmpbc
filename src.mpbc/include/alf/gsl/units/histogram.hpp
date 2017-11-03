#if 0
ln -sf $0 .$0.cpp && c++ -std=c++0x -I$HOME/prj .$0.cpp -Wall -Wfatal-errors -L$HOME/lib `pkg-config --libs gsl` -D_TEST_HISTOGRAM_UNITS_HPP -o ./.$0.x && ./.$0.x $@; exit
#endif
#ifndef HISTOGRAM_UNITS_HPP
#define HISTOGRAM_UNITS_HPP
#include "../../gsl/histogram.hpp"
#include<boost/units/quantity.hpp>
#include"boost/units/interval.hpp" //compat boost::units/boost::interval
namespace gsl{
using namespace boost::units;
template<class Unit>
class histogram::units : histogram
{
	public:
	typedef boost::numeric::interval<quantity<Unit> > interval;
	units(
		interval const& iv, 
		unsigned partition
	) : 
		histogram(
			histogram::interval(lower(iv).value(), upper(iv).value()),
			partition
		)
	{
	}
	bool operator()(quantity<Unit> const& x){
		return histogram::operator()(reinterpret_cast<double const&>(x));
	}
	bool operator()(quantity<Unit> const& x, double weight){
		return histogram::operator()(reinterpret_cast<double const&>(x), weight);
	}
	quantity<Unit> mean () const{return quantity<Unit>::from_value(histogram::mean ());}
	quantity<Unit> sigma() const{return quantity<Unit>::from_value(histogram::sigma());}
	using histogram::size;
	interval get_range(size_t i){
		return interval(
			quantity<Unit>::from_value(histogram::get_range(i).lower()), 
			quantity<Unit>::from_value(histogram::get_range(i).upper())
		);
	}
	//using histogram::operator[];
	double operator[](size_t i) const{return get(i);}
	using histogram::sum;
};

}
#endif //HISTOGRAM_UNITS_HPP

#ifdef _TEST_HISTOGRAM_UNITS_HPP
#include<iostream>
#include<boost/units/systems/si.hpp>
#include "../../gsl/random.hpp"
#include<iostream>
#include"boost/units/interval.hpp"
using namespace std;
using namespace boost::units;
int main(){
	typedef gsl::histogram::units<si::length> histogram;
	histogram 
	h(
		histogram::interval(-10.*si::meter, +10.*si::meter),
		10 
	);
	gsl::random::gaussian g(0.5, 3.);
	for(unsigned i=0; i!=100000; ++i){
		h(g()*si::meter);
	}
	clog << "mean is " << h.mean() << " sigma is "<< h.sigma() <<endl;
	return 0;
}
#endif

