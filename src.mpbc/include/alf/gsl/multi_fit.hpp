#if 0
ln -sf $0 $0.cpp && c++ $0.cpp `#-std=c++0x` `#echo -Wfatal-errors` -Wall `pkg-config --libs gsl` -I$HOME/prj -D_TEST_MULTI_FIT_HPP -o ./$0.x && ./$0.x $@
rm -f $0.cpp; exit
#endif
#ifndef MULTI_FIT_HPP
#define MULTI_FIT_HPP
#include "gsl/gsl_multifit_nlin.h"
#include<boost/array.hpp>
#include "boost/array_io.hpp"
#include<boost/function.hpp>
#include<boost/array.hpp>
#include<boost/noncopyable.hpp>
#include<boost/iterator/iterator_facade.hpp>
#include<vector>
#include<iostream> //for clog
#include<boost/math/special_functions/fpclassify.hpp>
using std::clog;
using std::endl;

namespace gsl{
namespace fit{

template<class Model>
class solver : 
	public boost::iterator_facade<
		solver<Model>, 
		Model const, 
		boost::forward_traversal_tag, 
		Model const /*noref*/ 
	>
{
	std::vector<std::pair<double,double> > const& data_;
	gsl_multifit_fdfsolver* pimpl_;
	gsl_multifit_function_fdf fdf_;
	Model cached_model_;
	public:
	solver(
		std::vector<std::pair<double,double> > const& data,
		boost::array<double, Model::static_size> const& x_init /*param*/ //or use the model itself? should use the model itself (see quantity case), but not implemented yet
	) : 
		data_(data), 
		pimpl_(gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, data_.size(), Model::static_size)),
		cached_model_(x_init)
		{
		assert(data_.size() > Model::static_size);
		gsl_vector_view x = gsl_vector_view_array(const_cast<double*>(&x_init[0]), Model::static_size);
		{
			fdf_.f      = &free_function_;
			fdf_.df     = &free_gradient_;
			fdf_.fdf    = &free_function_gradient_;
			fdf_.n      = data.size(); //data size, need to know size of data
			fdf_.p      = Model::static_size;
			fdf_.params = this;
		}
		gsl_multifit_fdfsolver_set(
			pimpl_, 
			&fdf_, //not copied! but always in scope
			&x.vector //copied (apparently)
		);
	}
	solver(solver const& other) :
		data_(other.data_),
		pimpl_(gsl_multifit_fdfsolver_alloc(gsl_multifit_fdfsolver_lmsder, data_.size(), Model::static_size)),
		fdf_(other.fdf_), //careful still need to change fdf_.param;
		cached_model_(other.position())
		{
		std::clog << "copied" << std::endl;
		gsl_vector_view x = gsl_vector_view_array(const_cast<double*>(&other.position()[0]), Model::static_size);
		fdf_.params = this;
		gsl_multifit_fdfsolver_set(
			pimpl_, 
			&fdf_, //not copied! but always in scope
			&x.vector
		);
	}
	~solver(){
		std::clog << "solver destroyed" << endl;
		gsl_multifit_fdfsolver_free(pimpl_);
	} //look at 
	boost::array<double, Model::static_size> const& position() const{
		return *(boost::array<double,Model::static_size>*)gsl_multifit_fdfsolver_position(pimpl_)->data;
	}
	int iterate(){return gsl_multifit_fdfsolver_iterate(pimpl_);}
	protected:
	Model dereference() const{ // iterator_facade implements operator*
		return Model(position());
	}
	void increment(){this->iterate();} //iterator_facade implements operator++
	friend class boost::iterator_core_access;
	public:
	std::string name(){return gsl_multifit_fdfsolver_name(pimpl_);}
	typedef gsl_multifit_fdfsolver_type* type;
	static type const lmder;
	static type const lmsder;
	static int free_function_(
		const gsl_vector* x /*adj param*/, 
		void* self_void     /*data and constant parameters*/, 
		gsl_vector* f       /*error at each point*/
	){
		//clog << "entering free_function" << endl;
		solver* self = (solver*)self_void;
		boost::array<double, Model::static_size> 
			param = *(boost::array<double, Model::static_size>*)(x->data);
		std::vector<double> f_as_vector(self->data_.size());
		Model m(param);
		for(unsigned i = 0; i != f_as_vector.size(); ++i){
			f_as_vector[i]=( m( self->data_[i].first) - self->data_[i].second );
			assert(not (boost::math::isnan)(f_as_vector));
		}
		for(unsigned i=0; i!=f_as_vector.size(); ++i){
			gsl_vector_set(f, i, f_as_vector[i]); // make it more efficient
		}
		return GSL_SUCCESS;
	}
	static int free_gradient_(
		const gsl_vector* x /*adj param*/, 
		void* self_void     /*data and constant parameters*/, 
		gsl_matrix* J       /*gradient of error at each point*/
	){
		//clog << "entering free_gradient" << endl;
		solver* self = (solver*)self_void;
		boost::array<double, Model::static_size> const& param = *(boost::array<double, Model::static_size>*)(x->data);
		Model m(param);
		std::vector<boost::array<double, Model::static_size> > J_as_vector(self->data_.size());
		for(unsigned i=0; i!=J_as_vector.size(); ++i){
			J_as_vector[i] = m.da(self->data_[i].first);
		}
		for(unsigned i=0; i!=J_as_vector.size(); ++i){
			for(unsigned j=0; j!= Model::static_size; ++j){
				gsl_matrix_set(J,i,j,J_as_vector[i][j]);
			}
		}
		return GSL_SUCCESS;
	}
	static int free_function_gradient_(
		const gsl_vector* x /*adj param*/, 
		void* self_void     /*data and constant parameters*/, 
		gsl_vector* f       /*error at each point*/,
		gsl_matrix* J       /*gradient of error at each point*/
	){
		//clog << "entering free_function_gradient" << endl;
		solver* self = (solver*)self_void;
		boost::array<double, Model::static_size> const& param = *(boost::array<double, Model::static_size>*)(x->data);
		Model m(param);
		std::vector<double> f_as_vector(self->data_.size());
		std::vector<boost::array<double, Model::static_size> > J_as_vector(self->data_.size());
		for(unsigned i = 0; i != f_as_vector.size(); ++i){
			f_as_vector[i]=( m( self->data_[i].first) - self->data_[i].second );
		}
		clog << endl;
		for(unsigned i=0; i!=f_as_vector.size(); ++i){
			gsl_vector_set(f, i, f_as_vector[i]); // make it more efficient
		}
		for(unsigned i=0; i!=J_as_vector.size(); ++i){
			J_as_vector[i] = m.da(self->data_[i].first);
		}
		for(unsigned i=0; i!=J_as_vector.size(); ++i){
			for(unsigned j=0; j!=Model::static_size; ++j){
				gsl_matrix_set(J,i,j,J_as_vector[i][j]);
			}
		}
		return GSL_SUCCESS;
	}
};
template<class Model> typename solver<Model>::type const solver<Model>::lmder  = gsl_multifit_fdfsolver_lmder;
template<class Model> typename solver<Model>::type const solver<Model>::lmsder = gsl_multifit_fdfsolver_lmsder;

template<unsigned NParam>
struct model_base{ //optionally you can inherit from this to make your model, to ensure proper syntax
	typedef boost::array<double, NParam> array;
	boost::array<double, NParam> const parameter;
	enum {static_size = NParam};
	model_base(boost::array<double, NParam> const& arr) : parameter(arr){}
	virtual double operator()(double const& x) const=0;
	virtual boost::array<double, NParam> da(double const& x) const=0;
};
}}
#endif

#ifdef _TEST_MULTI_FIT_HPP
#include<fstream>
#include<boost/assign.hpp>
using std::clog;
using std::endl;

struct exp_decay : gsl::fit::model_base<3>{
	double const& A_;
	double const& l_;
	double const& b_;
	exp_decay(boost::array<double, static_size> const& A_l_b) : gsl::fit::model_base<3>(A_l_b),
		A_(gsl::fit::model_base<3>::parameter[0]),
		l_(gsl::fit::model_base<3>::parameter[1]),
		b_(gsl::fit::model_base<3>::parameter[2]){
	}
	double operator()(double const& x) const{
		return A_*exp(-l_*x)+b_;
	}
	array da(double const& x) const{
		return (array){{
			exp(-l_*x),         // df/dA
			A_*exp(-l_*x)*(-x), // df/dlambda (NOT df/dx)
			1.                  // df/db
		}};
	}
};
/*
struct exp_decay{
	enum {static_size = 3};
	double A_;
	double l_;
	double b_;
	exp_decay(boost::array<double, static_size> const& A_l_b) :
		A_(A_l_b[0]),
		l_(A_l_b[1]),
		b_(A_l_b[2]){
	}
	double operator()(double const& x) const{
		return A_*exp(-l_*x)+b_;
	}
	boost::array<double,3> df(double const& x) const{
		return (boost::array<double,3>){
			exp(-l_*x),         // df/dA
			A_*exp(-l_*x)*(-x), // df/dlambda (NOT df/dx)
			1.                  // df/db
		};
	}
};*/
int main(){
	std::vector<std::pair<double, double> > data(30); //do not use push back
	{
		std::ofstream ofs("data.dat");
		for(unsigned i = 0; i < data.size(); i++){
			double t = i;
			data[i]=std::make_pair(
				t,
				5. * exp (-2. * t) + 1.0 
			);
		};
	}
	gsl::fit::solver<exp_decay> s(data, (boost::array<double, 3>){{5.2, 0.12, 1.0}});
	for(unsigned i=0; i<10; ++i){
		clog << "i = " << i << " {" << s.position()[0] << ", " << s.position()[1] << ", " << s.position()[2] << "}" << endl;
		++s;
		clog << (*s)(5.) << endl;
	}
	*boost::next(gsl::fit::solver<exp_decay>(data, (boost::array<double, 3>){{5.2, 0.12, 1.0}}),5.);
	return 0;
}

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
