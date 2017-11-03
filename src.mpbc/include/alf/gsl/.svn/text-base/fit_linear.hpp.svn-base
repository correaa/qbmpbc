#if 0
ln -sf $0 .$0.cpp && export LD_LIBRARY_PATH=$HOME/usr/lib; 
c++ `pkg-config --libs gsl` -D_TEST_FIT_LINEAR_HPP .$0.cpp -o ./$0.cpp.x && ./$0.cpp.x $@; exit;
#endif
#ifndef FIT_LINEAR_HPP
#define FIT_LINEAR_HPP

#include "gsl/gsl_fit.h"
#include<map>
#include<vector>
#include<cassert>
#include<iostream> //for clog
using std::clog;
using std::endl;

namespace gsl{
namespace model{
	class linear_{ /*linear model with constant term*/
		double c0_;
		double c1_;
		public:
		linear_(double const& c0, double const& c1) : c0_(c0), c1_(c1){}
		double operator()(double const& x) const{return c0_+x*c1_;}
		double const& c0() const{return c0_;}
		double const& c1() const{return c1_;}
	};
}
namespace linear{
	typedef model::linear_ model;
}
namespace fit{
	typedef std::pair<double, double> pair;
	typedef std::map<double, double> map;
	linear::model linear(std::vector<pair> const& vp){
		double c0, c1;{
			double cov00, cov01, cov11, sumsq;
			gsl_fit_linear(
				&vp[0].first, /*const size_t xstride*/ 2, &vp[0].second, /*const size_t xstride*/ 2, vp.size(), 
					&c0   , &c1, 
					&cov00, /**/
					&cov01, &cov11, 
					&sumsq
			);
		}
		return linear::model(c0, c1);
	}
	linear::model linear(map const& m){ 
		assert(m.size()!=0);
		std::vector<pair> vp(m.size());
		std::copy(m.begin(), m.end(), vp.begin()); //gsl needs to have a linear array anyway, copying is not harmful
		return fit::linear(vp);
	}
}
namespace linear{
	model fit(std::vector<std::pair<double, double> > const& pairs){
		return fit::linear(pairs);
	}
}
}
#ifdef _TEST_FIT_LINEAR_HPP
int main() {
	std::vector<
		std::pair<double, double> 
	> data(30); //do not use push back
	for(int i = 0; i < data.size(); i++) {
		double t = i;
		data[i] = std::make_pair(
		              t,
		              1.1 + 5.2 * t  //+ (i%10)/100.
		          );
	};
	gsl::linear::model l=gsl::fit::linear(data);
	clog<<l.c0()<<" "<<l.c1()<<endl;
	return 0;
}
#endif //TEST
#endif
// ~/usr/bin-hera/style --brackets=attach --indent=tab --indent-col1-comments --delete-empty-lines --add-brackets --keep-one-line-statements --convert-tabs --align-pointer=type --unpad-paren fit_linear.hpp
// Editor modelines  -  http://www.wireshark.org/tools/modelines.html
// Local variables:
// c-basic-offset: 4
// tab-width: 4
// indent-tabs-mode: t
// truncate-lines: 1
// End:
/* vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent: */
