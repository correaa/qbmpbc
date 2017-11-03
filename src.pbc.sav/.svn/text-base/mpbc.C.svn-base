#include<boost/lambda/lambda.hpp>
#include<boost/lambda/bind.hpp>
#include<boost/array_io.hpp>
using namespace boost::lambda;
using std::complex;
using boost::multi_array;


template<class MultiArray1, class MultiArray2>
void sum_to(MultiArray1& ma, MultiArray2 const& other){
	complex<double> const* j=other.origin();
	for(complex<double>* i=ma.origin(); i!=ma.origin()+ma.num_elements(); ++j, ++i){
		(*i)+=(*j);
	}
}

template<class MultiArray>
void set_zero(MultiArray& ma){
	for(complex<double> *i=ma.origin(); i!=ma.origin()+ma.num_elements(); ++i)
		(*i)=0;
}

template<class MultiArray>
double apply_kx2(MultiArray const& in, multi_array<complex<double>,3>& out, double a){
	out=in;
	double kin_en=0;
	for(int kxi=0; kxi!=out.shape()[0]; ++kxi){
		int kxiw=kxi<out.shape()[0]/2? kxi: kxi-out.shape()[0];
		double kx=kxiw*(2.*M_PI)/a;
		for(complex<double>* i=out[kxi].origin(); i!=out[kxi].origin()+out[kxi].num_elements(); ++i){
			kin_en+=norm(*i)*kx*kx*0.5;
			(*i)*=(kx*kx*0.5);
		}
	}
	return kin_en;
}

template<class MultiArray>
double apply_ky2(MultiArray const& in, multi_array<complex<double>,3>& out, double a){
	out=in;
	double kin_en=0;
	for(int kxi=0; kxi!=out.shape()[0]; ++kxi){
		for(int kyi=0; kyi!=out.shape()[1]; ++kyi){
			int kyiw=kyi<out.shape()[1]/2? kyi: kyi-out.shape()[1];
			double ky=kyiw*(2.*M_PI)/a;
			for(complex<double>* i=out[kxi][kyi].origin(); i!=out[kxi][kyi].origin()+out[kxi][kyi].num_elements(); ++i){
				kin_en+=(norm(*i)*ky*ky*0.5);
				(*i)*=(ky*ky*0.5);
			}
		}
	}
	return kin_en;
}

template<class MultiArray>
double apply_kz2(MultiArray const& in, multi_array<complex<double>,3>& out, double a){
	out=in;
	double kin_en=0;
	for(int kxi=0; kxi!=out.shape()[0]; ++kxi){
		for(int kyi=0; kyi!=out.shape()[1]; ++kyi){
			for(int kzi=0; kzi!=out.shape()[2]; ++kzi){
				int kziw=kzi<out.shape()[0]/2? kzi : kzi-out.shape()[0];
				double kz=kziw*(2.*M_PI)/a;
				kin_en+=norm(out[kxi][kyi][kzi])*(kz*kz*0.5);
				out[kxi][kyi][kzi]*=(kz*kz*0.5);
			}
		}
	}
	return kin_en;
}


template<class MultiArray>
double apply_kymBz2(MultiArray const& in, multi_array<complex<double>,3>& out, double a, double B){
	out=in;
	fftw3::normalized_dft_many(1, out, out, fftw3::forward);
	double kin_en=0;
	for(int kxi=0; kxi!=out.shape()[0]; ++kxi){
		for(int kyi=0; kyi!=out.shape()[1]; ++kyi){
			int kyiw=kyi<out.shape()[1]/2? kyi: kyi-out.shape()[1];
			double ky=kyiw*(2.*M_PI)/a;
			for(int zi=0; zi!=out.shape()[2]; ++zi){
				double z=(zi*a/out.shape()[2]);
				kin_en+=norm(out[kxi][kyi][zi])*(pow(ky+B*z/137.,2)*0.5);
				out[kxi][kyi][zi]*=(pow(ky+B*z/137.,2)*0.5);
			}
		}
	}
	fftw3::normalized_dft_many(1, out, out, fftw3::backward);
	return kin_en;
}

template<class MultiArray>
double apply_mpbc_kymBz2pkz2(MultiArray const& in, multi_array<complex<double>,3>& out, double a, double B){
	clog<<"apply_mpbc_kymBz2pkz2"<<endl;
	out=in;
	double kin_en=0;
	
	multi_array<complex<double>,3> out3(boost::extents[out.shape()[0]][out.shape()[1]][out.shape()[2]]);
	out3.assign(out.origin(), out.origin()+out.num_elements());	
	fftw3::normalized_dft_many(1, out, out, fftw3::forward);
	multi_array<complex<double>,2> out2(boost::extents[out.shape()[0]][out.shape()[1]*out.shape()[2]]);
	out2.assign(out.origin(), out.origin()+out.num_elements());

	//std::cout<<"B="<<B<<" (line 101)"<<std::endl; 
	for(int kxi=0; kxi!=out.shape()[0]; ++kxi){
		for(int kyi=0; kyi!=out.shape()[1]; ++kyi){
			int kyiw=kyi<out.shape()[1]/2? kyi: kyi-out.shape()[1];
			double ky=kyiw*(2.*M_PI)/a;
			for(int zi=0; zi!=out.shape()[2]; ++zi){
				double z=(zi*a/out.shape()[2]);
				//double zhat = ky + B*z/137.;
				//kin_en+=norm(out3[kxi][out.shape()[1]*kyi+zi])*(pow(zhat,2)*0.5);
				//out3[kxi][out.shape()[1]*kyi+zi]*=(pow(zhat,2)*0.5);
				kin_en+=norm(out[kxi][kyi][zi])*pow(ky-z*B/137.,2)*0.5;
				out[kxi][kyi][zi]*=(pow(ky-z*B/137.,2)*0.5);
			}
		}
	}
	/*
	for(int kxi=0; kxi!=out.shape()[0]; ++kxi){
		for(int kyi=0; kyi!=out.shape()[1]; ++kyi){			
			for(int kzi=0; kzi!=out.shape()[2]; ++kzi){				
				int kziw=kzi<out.shape()[2]/2? kzi: kzi-out.shape()[1];
				double kz=kziw*(2.*M_PI)/a;
				kin_en+=norm(out3[kxi][kyi][kzi])*pow(kz,2)*0.5;
				out3[kxi][kyi][kzi]*=(pow(kz,2)*0.5);
			}
		}
	}*/

	fftw3::normalized_dft_many(1, out2, out2, fftw3::forward);
	std::cout<<"acc after  "<<accumulate(out2[20].begin(), out2[20].end(), double(0.), _1+bind(&std::norm<double>, _2))<<std::endl;
	//for(complex<double>* i=out2.origin(); i!=out2.origin()+out2.num_elements(); ++i){
	//	(*i)/=sqrt((double)out2.shape()[1]);
	//}
	for(int kxi=0; kxi!=out2.shape()[0]; ++kxi){		
		for(int kzhati=0; kzhati!=out2.shape()[1]; ++kzhati){
			int kzhatiw=kzhati; //>out2.shape()[1]/2? kzhati: kzhati-out2.shape()[1];
			double kzhat=kzhatiw*(2.*M_PI)/(a*out2.shape()[0]*out2.shape()[0]);
			kin_en+=norm(out2[kxi][kzhati])*(pow(kzhat,2)*0.5);
			out2[kxi][kzhati]*=(pow(kzhat,2)*0.5);
		}
	}
	fftw3::normalized_dft_many(1, out2, out2, fftw3::backward);	
	//set_zero(out2);
	
	//for(complex<double>* i=out2.origin(); i!=out2.origin()+out2.num_elements(); ++i){
	//	(*i)/=sqrt((double)out2.shape()[1]);
	//}
	//fftw3::normalized_dft_many(1, out2_notuw, out2_notuw, fftw3::backward);
	/*
	for(int kxi=0; kxi!=out3.shape()[0]; ++kxi){
		for(int kyi=0; kyi!=out.shape()[1]; ++kyi){
			int kyiw=kyi<out.shape()[1]/2? kyi: kyi-out.shape()[1];
			double ky=kyiw*(2.*M_PI)/a;
			for(int kzi=0; kzi!=out.shape()[2]; ++kzi){
				int kziw=kzi<out.shape()[1]/2? kzi: kzi-out.shape()[1];
				double kz=kziw*(2.*M_PI)/a;
				kin_en+=norm(out2_notuw[kxi][kyi][kzi])*(pow(kz,2)*0.5);
				out2_notuw[kxi][kyi][kzi]*=(pow(kz,2)*0.5);
			}
		}
	}
	fftw3::normalized_dft_many(1, out2_notuw, out2_notuw, fftw3::forward);
	*/
	//sum_to(out2, out3);
	//multi_array_ref<complex<double>,3> out2_ref(out2.origin(), boost::extents[out.shape()[0]][out.shape()[1]][out.shape()[2]]);
	//set_zero(out);
	//sum_to(out, out2_ref);
	//sum_to(out, out2_notuw);
	fftw3::normalized_dft_many(1, out, out, fftw3::backward);
	

	
	//sum_to(out, out3);
	sum_to(out, out2);
	return kin_en;
}
