#ifdef COMPILATION_INSTRUCTIONS
	rm -f .$0.cpp.x && ln -sf $0 .$0.cpp && c++ -Wfatal-errors -Wall .$0.cpp -L$HOME/lib `pkg-config --libs gsl` -D_TEST_GSL_FUNCTION_HPP -o .$0.x && ./.$0.x $1 $2 $3 $4 $5 $6 $7 $8 $9
	rm -f $0.cpp.x $0.cpp
	exit
#endif
#ifndef GSL_FUNCTION_HPP_
#define GSL_FUNCTION_HPP_
#include<boost/function.hpp>
#include<gsl/gsl_math.h> //gsl_function
namespace gsl{
class function : public boost::function<double(double)>, public gsl_function{
	public:
	template<class F>
	function(F function_object) : 
		boost::function<double(double)>(function_object){
		gsl_function::function = &function::free_function_;
		gsl_function::params = this;
		assert(not (this->empty()));
	}
	function(function const& f) : 
		boost::function<double(double)>((boost::function<double(double)> const&)f){
		gsl_function::function = &function::free_function_;
		gsl_function::params = this;
		assert(not (this->empty()));
	}
	//private:
	static double free_function_(double x, void* self){
		assert(not ((gsl::function*)self)->empty());
		return ((function*)self)->operator()(x);
	}
};
}
#endif
#ifdef _TEST_GSL_FUNCTION_HPP
#include<iostream> 
using std::cout; using std::endl;
double freef(double d){return d*d-5.;}
int main(){
	gsl::function f(freef);
	cout << f(2.2) <<endl;
	return 0;
}

#endif

// astyle --brackets=attach --indent=tab --indent-col1-comments --delete-empty-lines --add-brackets --keep-one-line-statements --convert-tabs --align-pointer=type --unpad-paren fit_linear.hpp
// Editor modelines  -  http://www.wireshark.org/tools/modelines.html
// Local variables:
// c-basic-offset: 4
// tab-width: 4
// indent-tabs-mode: t
// truncate-lines: 1
// End:
/* vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent: */

