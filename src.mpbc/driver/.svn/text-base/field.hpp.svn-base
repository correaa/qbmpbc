#ifndef FIELD_HPP_
#define FIELD_HPP_
#include "UnitCell.h"
#include "qbox/cell.hpp"
namespace qbox{
  //  typedef UnitCell cell;
  typedef D3vector vector;
}
using qbox::vector;

template<typename T, typename D=boost::multi_array<T, 3u> >			//needs multi_array.hpp
class field{
	public:
	class grid{ //fft grid
		public:
		grid(UnitCell const& uc, boost::array<unsigned, 3> const& s) : domain_(uc), shape_(s){}
		UnitCell domain_;
		boost::array<unsigned, 3> shape_;
	};

	field(grid const& g) : 
		domain_(g.domain_), 
		data_(boost::extents[g.shape_[0]][g.shape_[1]][g.shape_[2]]){
			for (T* it=data_.origin(); it!=data_.origin()+data_.num_elements(); ++it)
				*it=0;
		}
	T operator()(boost::array<double, 3> const& a) const{
		return operator()(qbox::vector(a[0],a[1],a[2]));
	}
	T	operator()(qbox::vector const& v) const{
		//trilinear interpolation
		boost::array<int,3> idxi;
		boost::array<double,3> idxd;
		for(size_t i=0; i!=3; ++i){
			boost::tie(idxi[i],idxd[i])=alf::modf((domain_.b(i)*v)*data_.shape()[i]);
		}
		int const& xf = idxi[0]; 
		int const& yf = idxi[1];
		int const& zf = idxi[2];
		int const  xc = (xf+1)%data_.shape[0];
		int const  yc = (yf+1)%data_.shape[1];
		int const  zc = (zf+1)%data_.shape[2];
		double const& xd = idxd[0];
		double const& yd = idxd[1];
		double const& zd = idxd[2];
		//from http://en.wikipedia.org/wiki/Trilinear_interpolation#Method
		T i1=data_[xf][yf][zf]*(1-zd) + data_[xf][yf][zc]*zd;
		T i2=data_[xf][yc][zf]*(1-zd) + data_[xf][yc][zc]*zd;
		T j1=data_[xc][yf][zf]*(1-zd) + data_[xc][yf][zc]*zd;
		T j2=data_[xc][yc][zf]*(1-zd) + data_[xc][yc][zc]*zd;
		T w1=i1*(1-yd)+i2*yd;
		T w2=j1*(1-yd)+j2*yd;
		return w1*(1-xd)+w2*xd;
	}
	UnitCell const& domain() const{return domain_;}	
	UnitCell domain_;
	boost::multi_array<T, 3>  data_;
};



template<class T>
hdf5::file const& operator<<(hdf5::file const& f, hdf5::nvp<field<T> > const& nf) {
	using hdf5::make_nvp;
	using hdf5::nvp;
	boost::const_multi_array_ref<double,2> cell_as_ma(nf.const_value().domain().amat(), boost::extents[3][3]);
	std::clog<<boost::shape(nf.const_value().data_);
	boost::const_multi_array_ref<double,3> den(nf.const_value().data_);
	//std::clog<<boost::shape(den);
	f["cell"]<<cell_as_ma;
	f["data"]<<den;
//	f
//		<<nvp<boost::const_multi_array_ref<double, 2> >("cell", cell_as_ma)
//		<<nvp<boost::const_multi_array_ref<double,3> >("data", den);
		
	return f;
}
#endif

