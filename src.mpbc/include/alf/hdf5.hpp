#if 0
cp $0 $0.cpp
c++ -Wall -Wl,-rpath=$HOME/usr/lib -I$HOME/usr/include -D_TEST_HDF5_HPP $0.cpp -L/home/correaa/usr/lib -lhdf5 -lboost_system && ./a.out && rm -f $0.cpp a.out
exit
#endif

#ifndef _HDF5_HPP_
#define _HDF5_HPP_

#include<iostream>
#include<boost/multi_array.hpp>
//#include "boost/array_io.hpp"
#include "boost/iostreams/device/array.hpp"
#include "boost/multi_array_utility.hpp"
#include "boost/array_io.hpp"
#include<boost/noncopyable.hpp>
#include<boost/serialization/nvp.hpp>
#include<boost/lexical_cast.hpp>
#include<boost/filesystem.hpp>
//namespace hdf5{
#include<hdf5.h>
//}

using std::clog;
using std::endl;
using boost::filesystem::path;

namespace boost{
  template<class MultiArray>
  array<hsize_t, MultiArray::dimensionality>
  shape_hsize_t(MultiArray const& ma){
    array<hsize_t, MultiArray::dimensionality> ret;
    ret=shape(ma);
    return ret;
  }
//template<class MultiArray>
//hsize_t* shape_hsize_t_carr(MultiArray const& ma){
//	return 
//}
}

namespace hdf5{
namespace{
	struct file;

  using boost::serialization::nvp;  //use as hdf5::nvp<type>(name, type_instance);
  using boost::serialization::make_nvp;  //use as hdf5::make_nvp(name, type_instance); //no need to spec type
	using boost::array;
  using boost::lexical_cast;

  template<class T>
  nvp<T> make_const_nvp(char const* s, T const& t){return nvp<T>(s, const_cast<T&>(t));}

  namespace access{
    unsigned const truncate=H5F_ACC_TRUNC;
    unsigned const read_only=H5F_ACC_RDONLY;
  }
  struct handle{
    handle(hid_t _id) : id(_id){
      if(id<0) throw std::runtime_error("hdf5::handle '"+boost::lexical_cast<std::string>(id)+"' could not be constructed");
    }
    hid_t id;
  };  
  struct type : handle{
    type(int const& _id) : handle(_id){}
    template<class T>
    type(T*) : handle(H5Tcreate (H5T_COMPOUND, sizeof(T))){}
    type const& insert(std::string name, size_t offset, const int& typid) const{
      H5Tinsert (id, name.c_str(), offset, typid);
      return *this;
    }
    template<typename T>
    type const& insert(std::string name, T* offset) const{
      return insert(name, (size_t)offset, type::native<T>::id()); 
    }
    template<typename T, class C>
        type const& insert(std::string name, T C::* ptr_toCmember) const{
      return insert(name, (size_t)ptr_toCmember, type::native<T>::id());
    }
    template<typename T, class C>
    type const& insert(std::string name,  T& (C::*memfunc)() ) const{ //new trick! using member function to access member offset //works? gcc4.4.2
      return insert(name, (size_t)(&(((C*)0)->*memfunc)()),type::native<T>::id());
    }
    template<typename T> 
    struct native{
	  	static int id();
      //static int const id;
    };
  };
	template<> struct type::native<int>   {static int id(){return H5T_NATIVE_INT;}};
	template<> struct type::native<double>{static int id(){return H5T_NATIVE_DOUBLE;}};
	template<> struct type::native<std::complex<double> >{
		static int id(){
			return 
				(hdf5::type((std::complex<double>*)0)
					.insert("real", (size_t)((double*)0+0), type::native<double>::id() )
					.insert("imag", (size_t)((double*)0+1), type::native<double>::id() )).id;
		}
	};
  using boost::detail::multi_array::extent_gen;
  struct property_list : handle{
    property_list(property_list const& other) : handle(H5Pcopy(other.id)){}
  //template<std::size_t Rank>
  //proterty_list& set_chunk(boost::array<hsize_t, Rank> const& shape){
  //	H5Pset_chunk(this->id, Rank, shape.data());
  //	return *this;
  //}
    property_list& set_deflate(int level){ //compression level as in gzip
      assert(level>=0 and level<=10);
      H5Pset_deflate(this->id, level);
      return *this;
    } 
    ~property_list(){H5Pclose(this->id);}
  	protected:
    property_list(int _id) : handle(_id){}
  };
  struct dataset_creation_property_list : property_list{
    dataset_creation_property_list() : property_list(H5Pcreate(H5P_DATASET_CREATE)){}
  };
	struct check_type{check_type(herr_t check_code){if(check_code<0) throw std::logic_error("HDF5 error code");}};
  struct dataspace : handle{
    dataspace(int _id) : handle(_id){}
		template<typename MultiArray>
		dataspace(MultiArray const& ma) try : handle(H5Screate_simple(MultiArray::dimensionality, (hsize_t*)(&(shape_hsize_t(ma)[0])), NULL)){
    }catch(...){throw std::logic_error("hdf5 can create dataspace");}
		~dataspace(){H5Sclose(id);}
		template<typename MultiArrayView>
		check_type select_set_hyperslab(MultiArrayView const& mav){
			array<hsize_t, MultiArrayView::dimensionality> offset; for(unsigned i=0; i!=MultiArrayView::dimensionality; ++i) offset[i]=0;//mav.offset();
			array<hsize_t, MultiArrayView::dimensionality> stride; for(unsigned i=0; i!=MultiArrayView::dimensionality; ++i) stride[i]=mav.strides()[i];
			array<hsize_t, MultiArrayView::dimensionality> count ; for(unsigned i=0; i!=MultiArrayView::dimensionality; ++i) count[i] =mav.shape()[i];
			array<hsize_t, MultiArrayView::dimensionality> block ; for(unsigned i=0; i!=MultiArrayView::dimensionality; ++i) block[i] =1;
			return H5Sselect_hyperslab(this->id, H5S_SELECT_SET, &(offset[0]), &(stride[0]), &(count[0]), &(block[0]));
		}
    template<std::size_t Rank>
    extent_gen<Rank> get_extent(){
      assert(H5Sget_simple_extent_ndims(id)==Rank);
      extent_gen<Rank> ret;
      boost::array<hsize_t, Rank> shape;
      H5Sget_simple_extent_dims(id, &shape[0], NULL);
      for(unsigned i=0; i!=Rank; ++i)
        ret.ranges_[i]=typename extent_gen<Rank>::range(0, shape[i]);
      return ret;
    }
		H5S_sel_type get_select_type() const{ //returns the type of selection (all, hyperslab, etc) not the type of selection elements
			return H5Sget_select_type(this->id);
		}
  };
  struct datatype : handle{
    datatype(int typid) : handle(typid){}
    datatype(datatype const& other) : handle(H5Tcopy(other.id)){}
    H5T_class_t get_class(){return H5Tget_class(id);}
    void        set_order(H5T_order_t orderid){assert(H5Tset_order(id, orderid));}
    H5T_order_t get_order(){return H5Tget_order(id);}
    size_t      get_size(){return H5Tget_size(id);}
  };
  struct dataset : handle{
    dataset(file const& f, std::string const& name, datatype const& t, dataspace const& s);
//  dataset(file const& f, std::string const& name, datatype const& t, dataspace const& s) 
//	: handle(H5Dcreate(f.id, name.c_str(), t.id, s.id, H5P_DEFAULT)){}
		template<class MultiArray>
		dataset(file const& f, std::string const& name, MultiArray const& ma);
    dataset(file const& f, std::string const& name, datatype const& t, dataspace const& s, dataset_creation_property_list const& p);
    dataset(dataset const& other);
    ~dataset(){H5Dclose(id);}
		void write(int typid, int mem_dataspace_id, int file_dataspace_id, void const* in){H5Dwrite(id, typid, mem_dataspace_id, file_dataspace_id, H5P_DEFAULT, in);}
		void write(int typid, int mem_dataspace_id, void const* in){H5Dwrite(id, typid, mem_dataspace_id, H5S_ALL, H5P_DEFAULT, in);}
		template<class MultiArray, class MultiArrayView>
		void write(MultiArray const& ma, MultiArrayView const& mav){
			dataspace mads(ma); mads.select_set_hyperslab(mav);
			H5Dwrite(this->id, type::native<typename MultiArray::element>::id(),
			mads.id, H5S_ALL, H5P_DEFAULT, mav.origin());
		}
		template<class MultiArray>
		void write(MultiArray const& ma){
			dataspace mads(ma);
			H5Dwrite(
				this->id, type::native<typename MultiArray::element>::id(),
				mads.id, H5S_ALL, H5P_DEFAULT, ma.origin()
			);
		}
    //void write(int typid, void const* in){
		//	if(H5Dwrite(id, typid, H5S_ALL, H5S_ALL, H5P_DEFAULT, in)) 
		//		throw std::runtime_error("cannot write dataset "+boost::lexical_cast<std::string>(id));
    //}
		/*
		template<typename MultiArray>
		void write(MultiArray const& ma){
			dataspace ds(ma);
			H5Dwrite(
				id, 
				type::native<typename MultiArray::element>::id(), 
				H5S_ALL, H5S_ALL, H5P_DEFAULT, ma.origin()
			);
		}*/
//		template<class T>
//    void write(T const* in){return write(hdf5::type::native<T>::id(), in);}
    void read(int const& typid, dataspace const& m, dataspace const& s, void* out) const{
      if(H5Dread(id, typid, m.id, s.id, H5P_DEFAULT, out)) 
				throw std::runtime_error("cannot read dataset "+boost::lexical_cast<std::string>(id));
    }
    template<unsigned Rank>
    void read(int const& typid, boost::array<hsize_t, Rank> const& shape, hdf5::dataspace const& s, void* out) const{
      return read(typid, dataspace(shape), s, out);
    }
//    static dataset open(file const& f, std::string const& name){return dataset(open_tag(), f, name);}
    datatype get_type() const{return H5Dget_type(id);}
    dataspace get_space() const{return H5Dget_space(id);}
  private:
    struct open_tag{};
//    dataset(open_tag, file const& f, std::string const& name)
//		: handle(H5Dopen(f.id, name.c_str())){
//			if(id<0) 
//				throw std::runtime_error("hdf5::dataset with name '"+name+"' not found in hdf5::file id "+boost::lexical_cast<std::string>(f.id));
//    }
  };
	//private use
	class raw_dataset{
		friend struct file;
		file const* f_;
		std::string name_;
		raw_dataset(file const& f, std::string const& name) : f_(&f), name_(name){}
		public:
		template<class MultiArray>
		dataset operator()(MultiArray const& ma) const{return dataset(*f_, name_, ma);}
		template<class MultiArray>
		raw_dataset& operator<<(MultiArray const& ma);
		template<class MultiArray, class MultiArrayView>
		raw_dataset& write(MultiArray const& ma, MultiArrayView const& mav);
	};
  struct file : handle{
    file(path const& p, unsigned access=access::truncate)
		: handle(H5Fcreate(p.string().c_str(), access, H5P_DEFAULT, H5P_DEFAULT)){}
    ~file(){H5Fclose(id);}
    static file open(boost::filesystem::path const& p, unsigned access=access::read_only){
      return file(open_tag(), p, access);
    }
		raw_dataset operator[](std::string const& name) const{return raw_dataset(*this, name);}
  private:
    struct open_tag{};
    file(open_tag, path const& p, unsigned access) : handle(H5Fopen(p.string().c_str(), access,H5P_DEFAULT)){
      if(id<0) throw std::runtime_error("hdf5::file '"+p.string()+"' could not be open");
    }
  };
  dataset::dataset(file const& f, std::string const& name, datatype const& t, dataspace const& s) 
	: handle(H5Dcreate(f.id, name.c_str(), t.id, s.id, H5P_DEFAULT)){}
	template<class MultiArray>
	dataset::dataset(file const& f, std::string const& name, MultiArray const& ma)
	: handle(H5Dcreate(f.id, name.c_str(), type::native<typename MultiArray::element>::id(), dataspace(ma).id, H5P_DEFAULT)){}
  dataset::dataset(file const& f, std::string const& name, datatype const& t, dataspace const& s, dataset_creation_property_list const& p) 
	: handle(H5Dcreate(f.id, name.c_str(), t.id, s.id, p.id)){}
	template<class MultiArray>
	raw_dataset& raw_dataset::operator<<(MultiArray const& ma){
		hdf5::dataset ds(*f_, name_, ma);
		ds.write(ma);
		return *this;
	}
	template<class MultiArray, class MultiArrayView>
	raw_dataset& raw_dataset::write(MultiArray const& ma, MultiArrayView const& mav){
		hdf5::dataset ds(*f_, name_, mav);
		ds.write(ma, mav);
		return *this;
	}
}
}

template<typename T, unsigned NumDim>
void operator>>(hdf5::dataset const& ds, boost::multi_array<T, NumDim>& m){
  m.resize(ds.get_space().get_extent<boost::multi_array<T, NumDim>::dimensionality>());
  ds.read(hdf5::type::native<typename boost::multi_array<T, NumDim>::element>::id(), boost::shape_hsize_t(m), ds.get_space(), m.origin()); 
}

#ifdef _TEST_HDF5_HPP
using namespace boost;
using namespace hdf5;
using namespace std;
int main(){
/*
        multi_array<double, 3> a(extents[5][5][5]); 
        multi_array<double, 3> b(extents[5][5][5]);
        multi_array<double, 3> c(extents[5][5][5]);
        for(int i=0; i!=a.shape()[0]; ++i){
                for(int j=0; j!=a.shape()[0]; ++j){
                        for(int k=0; k!=a.shape()[0]; ++k){
                                a[i][j][k] = i+j+k+0.1;
                                b[i][j][k] = i+j+k+0.2;
                                c[i][j][k] = i+j+k+0.3;
                        }
                }
        }
        file("a.h5")<<make_nvp("a",a);  
        file("b.h5")<<make_nvp("b",b);
        file("c.h5")<<make_nvp("c",c);

        multi_array<array<double, 3>,3 > d(extents[5][5][5]);
        for(int l=0; l!=d.num_elements(); ++l){
                d.origin()[l][0]=a.origin()[l];
                d.origin()[l][1]=b.origin()[l];
                d.origin()[l][2]=c.origin()[l];
        }
        multi_array_ref<double, 4> d_ref(d.origin()->c_array(), extents[5][5][5][3]);
        file("d_ref.h5")<<make_nvp("d_ref",d_ref);
*/
				clog<<sizeof(boost::multi_array_types::size_type)<<"  "<<sizeof(hsize_t)<<endl;
				{
					multi_array<double, 1> ma(extents[8]);
					for(int i=0; i!=8; ++i) ma[i]=i;
//					file("ma.h5")<<make_nvp("ma",ma);
					typedef boost::detail::multi_array::multi_array_view<double, 1 > arrayview;
					typedef multi_array_types::index_range range;
					arrayview ma_odd  = ma[indices[range(0,8,2)]];
					arrayview ma_even = ma[indices[range(1,8,2)]];
					dataspace mads(ma);
					mads.select_set_hyperslab(ma_odd);
					dataspace fileS(ma_odd); //(ma_dd); 
					file mavf("mav.h5");

				  multi_array<double,2> cell_as_ma(boost::extents[3][3]);
					for(int i=0; i!=3; ++i) for(int j=0; j!=3; ++j)  cell_as_ma[i][j]=i+j;
					const_multi_array_ref<double, 2> cell_as_ma_ref(cell_as_ma.origin(), boost::extents[3][3]);

					mavf["cello"]<<cell_as_ma_ref;

					mavf["ma"]<<ma;
					mavf["ma_odd"].write(ma, ma[indices[range(1,8,2)]]);
					mavf["ma_even"].write(ma, ma[indices[range(0,8,2)]]);
//				  hdf5::dataset dset(
//    				mavf, 
//   				"ma_odd",
//    				hdf5::type::native<double>::id(),
//						fileS
//					);
				//		hdf5::dataset(mavf, "ma_odd", ma[indices[range(0,8,2)]])
				//			.write(ma, ma[indices[range(0,8,2)]]);
//					dset.write(
//						hdf5::type::native<double>::id(), 
//						mads.id,
//						ma_odd.origin()
//					);
				}
        return 0;
}
#endif

#endif

