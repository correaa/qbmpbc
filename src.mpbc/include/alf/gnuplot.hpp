#ifdef compile_instructions
cp $0 $0.cpp && c++ -Wall `#-Wfatal-errors` -I$HOME/usr/include -I$HOME/prj $0.cpp -Wl,-rpath=$HOME/usr/lib -L$HOME/usr/lib -lboost_filesystem -lboost_system -D_TEST_GNUPLOT_HPP -o ./$0.cpp.x && ./$0.cpp.x $1 $2 $3 $4 $5 $6 $7 $8 $9
rm -f $0.cpp.x $0.cpp
exit
#endif
#ifndef _GNUPLOT_HPP_
#define _GNUPLOT_HPP_

#ifndef _TEST_GNUPLOT_HPP
#define PHOENIX_LIMIT 18
#define FUSION_MAX_VECTOR_SIZE 18
#define BOOST_RESULT_OF_NUM_ARGS 18
#endif

#include<boost/filesystem.hpp>
#include<boost/filesystem/fstream.hpp>
#include<map>
#include<vector>
#include<string>
#include<sstream>
#include<boost/lexical_cast.hpp>
#include<boost/spirit/home/phoenix.hpp>
#include<boost/array.hpp>
#include<iostream>
#include<boost/tuple/tuple.hpp>
#include"boost/just.hpp"
#include "boost/units/phoenix.hpp"
#include<boost/math/special_functions/fpclassify.hpp>
namespace gnuplot{
	class file;
};
namespace gnuplot{

using namespace boost::phoenix;
using namespace boost::phoenix::arg_names;

using namespace boost::filesystem;
using std::pair;
using std::map;
using std::string;
using std::getline;
using boost::lexical_cast;
using std::runtime_error;
using std::clog;
using std::endl;
using boost::tuple;
typedef boost::array<double, 3> triple;

class file : path{
	public:
	file(path const& p):path(p){}
	protected:
	template<class Actor1, class Actor2>
	struct expression_use : private pair<actor<Actor1>,  actor<Actor2> >{
		//file const& const_self;
		mutable ifstream ifs_; //(const_self); 
		expression_use();
		expression_use(expression_use const&);
		expression_use(file const& f, pair<actor<Actor1>,  actor<Actor2> > ij) : pair<actor<Actor1>,  actor<Actor2> >(ij), ifs_(f){
			if(not ifs_) throw std::runtime_error("could not open file "+f.string());
		}
		template<class T1, class T2>
		std::pair<T1, T2> read_pair() const{
			std::string line;
			do{
				getline(ifs_, line); 
			}while(line[0]=='#');
			std::istringstream iline(line);
			boost::array<double, PHOENIX_ARG_LIMIT> d;
			for(unsigned idx=0; idx<d.size(); ++idx){
				d[idx] = std::numeric_limits<double>::signaling_NaN();
				if(not iline) continue;
				std::string entry;
				iline >> entry; //clog<<"read "<<idx<<": " <<entry<<" "; 
				try{
					d[idx] = boost::lexical_cast<double>(entry);
				}catch(...){
					d[idx] = std::numeric_limits<double>::signaling_NaN();
					//clog << "could not read "<<entry<<" as a number assuming "<<d[idx]<<endl;
				}
			}
			BOOST_STATIC_ASSERT(PHOENIX_LIMIT>=16); // eg. #define PHOENIX_LIMIT 18 #define FUSION_MAX_VECTOR_SIZE 18 #define BOOST_RESULT_OF_NUM_ARGS 18 at top of program
			T1 x ( (this->first )(d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9], d[10], d[11], d[12], d[13], d[14], d[15]) );
			T2 y ( (this->second)(d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9], d[10], d[11], d[12], d[13], d[14], d[15]) );
			return pair<T1,T2>(x,y);
		}
		//template<class T>
		//double value(T const& v){
		//	return v.value();
		//}
		double value(double const& d){return d;}
		template<class T1, class T2>
		operator std::vector<std::pair<T1, T2> >() const{
			std::vector<std::pair<T1, T2> > ret;
			while(ifs_){
				std::pair<T1, T2> p(read_pair<T1, T2>());
				if(not ifs_) break;
				//std::clog << "about to insert " << p.first << " " << p.second << std::endl;
				if((p.first!=p.first) or (p.second!=p.second)) continue;
				//std::clog << "inserted " << p.first << " " << p.second << std::endl;
				ret.push_back(p);
			}
			return ret;
		}
		template<class T1, class T2>
		operator std::map<T1, T2>() const{
			std::map<T1, T2> ret;
			while(ifs_){
				std::pair<T1, T2> p(read_pair<T1, T2>());
				if((boost::math::isnan)(p.first) or (boost::math::isnan)(p.second)) continue;
				ret.insert(p);
			}
			return ret;
		}
		operator map<double, double>() const{
			std::map<double,double> ret;
			while(ifs_){
				std::pair<double, double> p(read_pair<double, double>());
				if((boost::math::isnan)(p.first) or (boost::math::isnan)(p.second)) continue;
				ret.insert(p); //do not use ret[x]=y; it gets confused with NaNs
			}
			return ret;
		}
	};
	template<class Actor1, class Actor2, class Actor3>
	struct expression_use3 : private tuple<actor<Actor1>, actor<Actor2>, actor<Actor3> >, boost::just<file const&>::type{
		expression_use3(file const& f, tuple<actor<Actor1>, actor<Actor2>, actor<Actor3> > const& ijk) : tuple<actor<Actor1>, actor<Actor2>, actor<Actor3> >(ijk), boost::just<file const&>::type(f){}
		operator std::vector<triple>() const{
			std::vector<triple> ret;
			ifstream ifs(*this);
			for(std::string line; getline(ifs, line);){
				if(line[0]=='#' or 0==line.size()) continue;
				std::istringstream iline(line);
				boost::array<double, PHOENIX_ARG_LIMIT> d;
				for(unsigned idx=0; idx<d.size(); ++idx){
					d[idx] = std::numeric_limits<double>::signaling_NaN();
					if(not iline) continue;
					std::string entry;
					iline >> entry; //clog<<"read "<<idx<<": " <<entry<<" "; 
					try{
						d[idx] = boost::lexical_cast<double>(entry);
					}catch(...){
						d[idx] = std::numeric_limits<double>::signaling_NaN();
						//clog << "could not read "<<entry<<" as a number assuming "<<d[idx]<<endl;
					}
				}
				double x = (boost::get<0>(*this) )(d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9]);
				double y = (boost::get<1>(*this) )(d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9]);
				double z = (boost::get<2>(*this) )(d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9]);
				ret.push_back((triple){{x,y,z}});
			}
			return ret;
		}
	};
	struct use : protected pair<unsigned, unsigned>{
		file const& const_self;
		use(file const& f, pair<unsigned, unsigned> ij) : pair<unsigned, unsigned>(ij), const_self(f){}
		operator map<double, double>() const try{
		map<double, double> ret;
		ifstream ifs(const_self); if(not ifs) throw std::runtime_error("can open file '"+const_self.string()+"'");
		{
			unsigned line_counter=0;
			for(
				std::string line; 
				getline(ifs, line); 
				++line_counter)
			{
			pair<double, double> entry; entry.second=-1.; //for debug
			std::istringstream iline(line);
			if(this->first==0) entry.first=line_counter;
			if(this->second==0) entry.second=line_counter;
			{
				std::string in;
				unsigned column_counter=1;
				for(; ;++column_counter){
					if(this->first==column_counter){
						iline>>entry.first; //=lexical_cast<double>(in); //more permisive
					}else if(this->second==column_counter){
						iline>>entry.second; //=lexical_cast<double>(in);
					}else{
						iline>>in;
					}
					if(not iline){
						//if(column_counter<this->first or column_counter<this->second){
						//throw runtime_error("can not read fields from line");
						//}
						break;
					}
				}
				if(column_counter<=this->first or column_counter<=this->second){//		  clog<<"ignored line "<<line_counter<<" from "<<const_self.string()<<endl;
					continue;
				}
			}
			assert(entry.second!=-1.);
			ret[entry.first]=entry.second;
		}
		return ret;
	}
	}catch(std::runtime_error r){
		throw runtime_error("can not read file '"+const_self.string()+"', "+r.what());
      }
    };
  public:
	template<class A1, class A2>
	expression_use<A1, A2> 
	u(actor<A1> const& a1, actor<A2> const& a2){
		return expression_use<A1, A2>(*this, pair<actor<A1>, actor<A2> >(a1, a2));
	}
	template<class A1, class A2, class A3> expression_use3<A1, A2, A3> u(actor<A1> const& a1, actor<A2> const& a2, actor<A3> const& a3){return expression_use3<A1, A2, A3>(*this, tuple<actor<A1>, actor<A2>, actor<A3> >(a1, a2, a3));}
    use u(unsigned i, unsigned j) const{return use(*this, pair<unsigned, unsigned>(i,j) );}
    use using_(unsigned i, unsigned j) const{return u(i,j);}
  };
}
#endif
#ifdef _TEST_GNUPLOT_HPP

#include<boost/spirit/home/phoenix.hpp>
#include<boost/units/systems/si.hpp>
#include<boost/function.hpp>
using namespace boost::phoenix::arg_names;
using namespace boost::phoenix;
//using namespace boost::units;
using std::cout; using std::endl;
int main(){
	std::map<double, double> m = gnuplot::file("gnuplot_test.dat").u(_1,_3);
	for(std::map<double, double>::const_iterator it=m.begin(); it!=m.end(); ++it){
		cout << it->first <<" " <<it->second<<"\n";
	}
	std::vector<gnuplot::triple> v= gnuplot::file("gnuplot_test.dat").u(_1,_2,_3*2.);
	for(std::vector<gnuplot::triple>::const_iterator it=v.begin(); it!=v.end(); ++it){
		cout << (*it)[0] << " "<<(*it)[1]<<" "<<(*it)[2] <<endl;
	}
	boost::function<quantity<si::length>(double)> f = arg1 * (1.*si::meter);
	cout << f(1.) << endl;
	std::vector<std::pair<double, quantity<si::length> > > v2 = 
		gnuplot::file("gnuplot_test.dat").u(
			arg1, 
			arg3*(1.*si::meter)
		);
	std::map<double, quantity<si::length> > v2m = 
		gnuplot::file("gnuplot_test.dat").u(
			arg1, 
			arg3*(1.*si::meter)
		);
	cout << v2[2].first <<" " <<v2[2].second << endl;
	std::vector<std::pair<
		quantity<si::area>, 
		//double, 
		quantity<si::length> 
	> > v3 = 
		gnuplot::file("/home/correaa/prj/hydro/data/morales/pressure_T1000.dat").u(
			//arg1*arg1*(si::meter*si::meter), 
			arg1*(1.*si::meter*si::meter),
			arg3*(1.*si::meter)
		);
	std::clog << "read size " << v3.size() << std::endl;

	std::vector<std::pair<
		quantity<si::area>, 
		//double, 
		quantity<si::length> 
	> > v4 = 
		gnuplot::file("gnuplot.dat").u(
			//arg1*arg1*(si::meter*si::meter), 
			arg1*(1.*si::meter*si::meter),
			arg3*(1.*si::meter)
		);
	std::clog << "read size " << v4.size() << std::endl;
	for(unsigned i = 0; i != v4.size(); ++i){
		std::clog << i << " " << v4[i].first << " " << v4[i].second << std::endl;
	}
	return 0;
}
#endif

