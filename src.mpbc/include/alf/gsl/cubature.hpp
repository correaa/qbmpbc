#ifdef compile_instructions 
	ln -sf $0 .$0.cpp; c++ -std=c++0x -D_CUBATURE_TEST_HPP -I$HOME/prj/mag/Qbox/src.wei/include/alf .$0.cpp -o .$0.x -Wall && ./.$0.x ; exit
#endif
#include "../gsl/epsilon.hpp"
#include "../gsl/result.hpp"
#include<boost/array.hpp>
#include<boost/function.hpp>
#include<boost/numeric/interval.hpp>
#include "../gsl/cubature-20091002/cubature.c"
#include<boost/lexical_cast.hpp>
#include<iostream> //debug

namespace cubature{

namespace bn = boost::numeric;

template<size_t Dim>
struct array : boost::array<double, Dim>{
	typedef boost::array<double, Dim> base;
	typedef base type;
	array(base const& arr) : base(arr){}
	array(double const& d);
	operator double&(); //add const version?
	operator double const&() const;
};
template<> array<1>::array(double const& d) : base((base){{d}}){}
template<> array<1>::operator double&(){return (*this)[0];}
template<> array<1>::operator double const&() const{return (*this)[0];}

template<size_t NDim, size_t FDim>
struct function : //aka integrand 
#define BASE_TYPE boost::function</*typename*/ array<FDim> /*::type*/ (/*typename*/ array<NDim> /*::type*/ const&)> 
BASE_TYPE
{
	typedef BASE_TYPE base_type;
	template<class T> function(T function_object) : base_type(function_object){}
	function(function const& other){assert(0);}
	//private:
	static void free_function_(unsigned ndim, const double *x, void *self /*fdata*/, unsigned fdim, double *fval){
		assert(ndim=NDim);
		assert(fdim=FDim);
		array<FDim>* fval_arr = (array<FDim>*)fval;
		//typename array<FDim>::type* fval_arr = (typename array<FDim>::type*)fval;
		*fval_arr = ((function<NDim, FDim>*)self)->operator()(*(array<NDim> const*)x);
		//*fval_arr = ((function*)self)->operator()(*((array<NDim>::type const*)x));
	}
};
#undef BASE_CLASS

//template<size_t ND

// todo: repeated code for monte_carlo.hpp
template<size_t NDim> 
struct limits : boost::array<bn::interval<double>, NDim>{
	typedef boost::array<bn::interval<double>, NDim> base;
	limits(base const& arr) : base(arr){}
};

/// Genz-Malik algorithm for multidimensional vector-valued integrals, interface to cubature.adapt_intergrate http://ab-initio.mit.edu/wiki/index.php/Cubature
template<size_t NDim, size_t FDim>
boost::array<gsl::result, FDim> adaptive_integrate(
	function<NDim, FDim> const& f, 
	limits<NDim> const& l,
	unsigned maximum_evaluations = 100000, 
	gsl::epsilon const& eps_requiered = gsl::epsilon( 0, 1e-5 ) //remember this is a multidimensional integral small tolerances forces the code to do many evaluations
){
	boost::array<double, FDim> val; 
	boost::array<double, FDim> err;
	boost::array<double, NDim> xmin;
	boost::array<double, NDim> xmax;
	for(size_t i=0; i!=NDim; ++i){
		xmin[i] = lower(l[i]); 
		xmax[i] = upper(l[i]);
		//std::clog<<" ["<<xmin[i]<<", "<<xmax[i]<<"]"<<std::endl;
	}
	assert(FDim==1); assert(NDim==3);
	int code = adapt_integrate(
		FDim /*unsigned fdim*/, 
		function<NDim, FDim>::free_function_ /*integrand f*/, 
		const_cast<void*>((void*)&f) /*void *fdata*/,
		NDim /*unsigned dim*/, 
		xmin.data() /*const double *xmin*/, 
		xmax.data() /*const double *xmax*/,
		maximum_evaluations /*unsigned maxEval*/, 
		eps_requiered.absolute() /*double reqAbsError*/, 
		eps_requiered.relative() /*double reqRelError*/,
		val.c_array() /*double *val*/, 
		err.c_array() /*double *err*/
	);
	if(code==-1) throw std::runtime_error("cubature.adapt_intergrate code -1: probably didn't achieve the requiered tolerance with the maximum amount of evaluations");
	if(code!=0) throw std::runtime_error("cubature.adapt_intergrate failed with code "+boost::lexical_cast<std::string>(code));
	boost::array<gsl::result, FDim> ret; for(unsigned i=0; i!=FDim; ++i) ret[i]=gsl::result(val[i], err[i]);
	return ret;
}

template<size_t NDim>
gsl::result adaptive_integrate(
	function<NDim,1> const& f,
	limits<NDim> const& l
){
//	boost::array<gsl::result, 1> ret(adaptive_integrate<NDim, 1>(f, l));
	return adaptive_integrate<NDim, 1>(f, l)[0] ; //ret[0];
}


}//ns cubature

#ifdef _CUBATURE_TEST_HPP
#include<iostream>
#include<cmath>
using std::cout; using std::endl;

void cgaussian(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    double sigma = *((double *) fdata); // we can pass σ via fdata argument
    double sum = 0;
    for (unsigned i = 0; i < ndim; ++i) sum += x[i] * x[i];
    // compute the output value: note that fdim should == 1 from below
    fval[0] = exp(-sigma * sum);
}

boost::array<double, 1> gaussian(boost::array<double, 3> const& x){
    double sigma = 0.5;
    double sum = 0;
    for (unsigned i = 0; i < 3; ++i) sum += x[i] * x[i];
    boost::array<double,1> ret;
    ret[0]=exp(-sigma * sum);
    return ret;
}

double dgaussian(boost::array<double, 3> const& x){
    double sigma = 0.5;
    double sum = 0;
    for (unsigned i = 0; i < 3; ++i) sum += x[i] * x[i];
    return exp(-sigma * sum);
}

int main(){
	{
		double xmin[3] = {-2,-2,-2}, xmax[3] = {2,2,2}, sigma = 0.5, val, err;
		adapt_integrate(1, cgaussian, &sigma, 3, xmin, xmax, 0, 0, 1e-4, &val, &err);
		printf("Computed integral = %0.10g +/- %g\n", val, err);
	}
	{
		double x[3] = {1.,2.,3.};
		double fval[1]={0.};
		double sigma = 0.5;
		cgaussian(3, x, &sigma, 1, fval);
		cout << "fval "<< fval[0] << std::endl;
	}
	cubature::function<3,1> f(&gaussian);
	{
		double x[3] = {1.,2.,3.};
		double fval[1]={0.};
		//double sigma = 0.5;
		cubature::function<3,1>::free_function_(3, x, &f, 1, fval);
		cout << "fval2 "<< fval[0] << std::endl;
	}
	{
		double xmin[3] = {-2,-2,-2}, xmax[3] = {2,2,2}, val, err;
		adapt_integrate(1, cubature::function<3,1>::free_function_, &f, 3, xmin, xmax, 0, 0, 1e-4, &val, &err);
		printf("2 Computed integral = %0.10g +/- %g\n", val, err);
	}



	cout << "lop " << gaussian( (boost::array<double,3>){{1.,2.,3.}} )[0] <<std::endl;
	cout
		<< "oip " <<std::endl
		<< cubature::adaptive_integrate<3,1>(f, 
			(boost::array<boost::numeric::interval<double>,3>)
			{{
				boost::numeric::interval<double>(-2,2), 
				boost::numeric::interval<double>(-2,2), 
				boost::numeric::interval<double>(-2,2)
			}}
			//0,
			//gsl::epsilon(0, 1e-4)
		)[0];
	cout<<
		cubature::adaptive_integrate<3>(f, 
			(boost::array<boost::numeric::interval<double>,3>)
			{{
				boost::numeric::interval<double>(-2,2), 
				boost::numeric::interval<double>(-2,2), 
				boost::numeric::interval<double>(-2,2)
			}}
		)
		;
	cubature::function<3,1> fd(&dgaussian);
	std::cout << 
		cubature::adaptive_integrate<3>(fd, 
			(boost::array<boost::numeric::interval<double>,3>)
			{{
				boost::numeric::interval<double>(-2,2), 
				boost::numeric::interval<double>(-2,2), 
				boost::numeric::interval<double>(-2,2)
			}}
		);
	return 0;
}

#endif

