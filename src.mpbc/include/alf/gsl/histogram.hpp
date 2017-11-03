#ifdef compile_instructions 
ln -sf $0 .$0.cpp && c++ -D_HISTOGRAM_TEST_HPP -I$HOME/prj/mag/Qbox/src.wei/include/alf `pkg-config --libs gsl` .$0.cpp -o .$0.x -Wall && ./.$0.x $@; exit
#endif
#include<gsl/gsl_histogram.h>
#include<boost/numeric/interval.hpp>
#include<boost/noncopyable.hpp>
#include<vector>
#include<iostream>
using std::clog ; using std::endl;
namespace gsl{
	///< based (literally) http://www.gnu.org/software/gsl/manual/html_node/The-histogram-struct.html
	class histogram{
		public:
		gsl_histogram* pimpl_;
		typedef boost::numeric::interval<double> interval;
		//size_t n: This is the number of histogram bins
		//double * range: The ranges of the bins are stored in an array of n+1 elements pointed to by range.
		//double * bin: The counts for each bin are stored in an array of n elements pointed to by bin. The bins are floating-point numbers, so you can increment them by non-integer values if necessary. 
		histogram(size_t n) : pimpl_(gsl_histogram_alloc(n)){}
		histogram(interval const& iv, size_t n) : pimpl_(gsl_histogram_alloc(n)){
			set_ranges_uniform(iv);
		}
		histogram(histogram const& other) : pimpl_(gsl_histogram_clone(other.pimpl_)){}
		histogram& operator=(histogram const& other){
			if(&other!=this){ //todo: replace by safe swap
				//assert(size()==other.size());
				gsl_histogram_memcpy(pimpl_, other.pimpl_);
			}
			return *this;
		}
		size_t size() const{return pimpl_->n;}
		//std::vector<double> const& ranges() const{return ;}
		//std::vector<double> const& bin() const{return ;}
		void set_ranges(std::vector<double> const& r){
			assert(
				gsl_histogram_set_ranges(pimpl_, &r[0], r.size())
			);
		}
		void set_ranges_uniform(interval const& iv){
			gsl_histogram_set_ranges_uniform(pimpl_, lower(iv), upper(iv));
		}
		~histogram(){gsl_histogram_free(pimpl_);}
		//http://www.gnu.org/software/gsl/manual/html_node/Updating-and-accessing-histogram-elements.html
		bool increment(double x){return gsl_histogram_increment(pimpl_, x);}
		bool accumulate(double x, double weight){return gsl_histogram_accumulate(pimpl_, x, weight);}

		bool operator()(double const& x){return increment(x);}
		bool operator()(double const& x, double const& weight){return accumulate(x, weight);}

		private:
		struct proxy_histogram_subscript{
			histogram& self;
			double const& x;
			proxy_histogram_subscript(histogram& self, double const& x) : self(self), x(x){}
			proxy_histogram_subscript const& operator++() const{self.increment(x); return *this;}
			bool operator++(int) const{return self.increment(x);}
			bool operator+=(double const& weight) const{return self.accumulate(x, weight);}
		};
		public:
		proxy_histogram_subscript operator[](double const& x){
			//clog << "called "<<endl;
			return proxy_histogram_subscript(*this, x);
		}

		double get(size_t i) const{return gsl_histogram_get(pimpl_, i);}

		double operator[](size_t i) const{return get(i);}

		interval get_range(size_t i) const{
			double lower=1, upper=0;
			gsl_histogram_get_range (pimpl_, i, &lower, &upper);
			return interval(lower,upper);
		}
		double max() const{return gsl_histogram_max(pimpl_);}
		double min() const{return gsl_histogram_min(pimpl_);}
		size_t bins() const{return gsl_histogram_bins(pimpl_);}
		void   reset(){return gsl_histogram_reset(pimpl_);}
		histogram& operator=(int*){reset(); return *this;}
		size_t find(double x) const{
			size_t ret;
			assert(gsl_histogram_find(pimpl_, x, &ret));
			return ret;
		}
		interval range() const{return interval(min(), max());}

		// http://www.gnu.org/software/gsl/manual/html_node/Histogram-Statistics.html
		double max_val() const{return gsl_histogram_max_val(pimpl_);} //< max value contained (histogrammed?)
		size_t max_bin() const{return gsl_histogram_max_bin(pimpl_);} //< bin containing the max value
		double min_val() const{return gsl_histogram_min_val(pimpl_);} //< mininum value contained (histogrammed?)
		size_t min_bin() const{return gsl_histogram_min_bin(pimpl_);} //< bin containing the minimum value
		double mean() const{return gsl_histogram_mean(pimpl_);} //< mean of the histogrammed variable
		double sigma() const{return gsl_histogram_sigma(pimpl_);} //< standard deviation of the histogrammed variable
		double sum() const{return gsl_histogram_sum(pimpl_);} //< sum of all bin values (including negative ones)

		// http://www.gnu.org/software/gsl/manual/html_node/Histogram-Operations.html
		friend int equal_bins(histogram const& h1, histogram const& h2){return gsl_histogram_equal_bins_p(h1.pimpl_, h2.pimpl_);}
		friend int add(histogram& h1, histogram const& h2){return gsl_histogram_add(h1.pimpl_, h2.pimpl_);}
		friend int sub(histogram& h1, histogram const& h2){return gsl_histogram_sub(h1.pimpl_, h2.pimpl_);}
		friend int mul(histogram& h1, histogram const& h2){return gsl_histogram_mul(h1.pimpl_, h2.pimpl_);}
		friend int div(histogram& h1, histogram const& h2){return gsl_histogram_div(h1.pimpl_, h2.pimpl_);}
		friend int scale(histogram& h, double s){return gsl_histogram_scale(h.pimpl_, s);}
		friend int shift(histogram& h, double offset){return gsl_histogram_shift(h.pimpl_, offset);}
		histogram& operator+=(histogram const& other){
			assert(add(*this, other)); return *this;
		}
		histogram& operator-=(histogram const& other){
			assert(sub(*this, other)); return *this;
		}
		histogram& operator*=(histogram const& other){
			assert(mul(*this, other)); return *this;
		}
		histogram& operator/=(histogram const& other){
			assert(div(*this, other)); return *this;
		}
		histogram& operator*=(double s){
			assert(scale(*this, s)); return *this;
		}
		histogram& operator+=(double offset){
			assert(shift(*this, offset)); return *this;
		}
		/// Probability distribution function from histogram
		class pdf{
			gsl_histogram_pdf* const pimpl_;
			public:
			pdf(histogram const& h) : pimpl_(gsl_histogram_pdf_alloc(h.size())){
				gsl_histogram_pdf_init(pimpl_, h.pimpl_);
			}
			double sample(double random01) const{
				assert(random01>=0 and random01<1);
				return gsl_histogram_pdf_sample(pimpl_, random01);
			}
//			double operator()(double x){
//			}
			~pdf(){gsl_histogram_pdf_free(pimpl_);}
		};
		template<class Unit> class units;
	};
}
#ifdef _HISTOGRAM_TEST_HPP
#include "../gsl/random.hpp"
#include<iostream>
using namespace std;
int main(){
	gsl::histogram h(
		gsl::histogram::interval(-10,10), 
		100
	);
	gsl::random::gaussian g(0.5, 3.);
	for(unsigned i=0; i!=100000; ++i){
		//h(g());
		h[g()]++;
		h[g()]+=1.;
	}
	clog << "mean is " << h.mean() << " sigma is "<< h.sigma() <<endl;
	gsl::histogram::pdf p(h);
	for(unsigned i=0; i!=h.size(); ++i){
		//cout << median(h.get_range(i)) << ' ' << h[i] << ' ' << g.pdf(median(h.get_range(i))) << std::endl;
	}
	for(double x=0; x<1.; x+=0.01){
		cout << x << ' ' <<p.sample(x)<<'\n';
	}
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

