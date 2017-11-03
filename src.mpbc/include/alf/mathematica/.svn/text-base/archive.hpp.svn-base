#ifndef MATHEMATICA_ARCHIVE_
#define MATHEMATICA_ARCHIVE_
#include <boost/multi_array.hpp>
#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <boost/serialization/nvp.hpp>
#include <boost/algorithm/string.hpp> // boost::replace_last
#include <boost/spirit/core.hpp>
#include <boost/spirit/iterator/multi_pass.hpp>
#include <boost/spirit/actor/push_back_actor.hpp>
#include <boost/spirit/phoenix.hpp>
#include <boost/spirit/attribute/closure.hpp>
#include <boost/spirit/utility/lists.hpp>
#include <boost/lexical_cast.hpp> 
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>
using boost::multi_array;
using boost::multi_array_ref;
using std::vector;
using std::clog;
using std::endl;

namespace mathematica{
	//parser code starst here (to be cleaned up)
	using namespace boost::spirit; 											// spirit parser is beautiful
	using namespace phoenix;

	struct real_prsr : parser<real_prsr>{
		typedef real_prsr self_t;										 // utility
		real_prsr(){}													 // default constructor inhibits copy constructor?
		template <typename ScannerT> struct result{						 // be explicit on the result type of the parser for other classes (otherwise char*)
        	typedef typename match_result<ScannerT, double>::type type;
		};
		template<typename ScannerT>										 // parser function
		typename parser_result<self_t, ScannerT>::type
		parse(ScannerT const& scan) const{
			typename ScannerT::iterator_t save = scan.first;			 //save begin of parsing
			typedef typename parser_result<boost::spirit::real_parser<>, ScannerT>::type real_match_t;
			double ret; int ret_exp=0;
			int count=(boost::spirit::real_p[assign_a(ret)] >> !(ch_p('*') >> '^' >> boost::spirit::real_p[assign_a(ret_exp)])).parse(scan).length();
			if(ret_exp!=0) ret*=pow(10.,ret_exp);
			return scan.create_match(count, ret, save, scan.first);
		}
    } real_p;
   	struct complex_prsr : parser<complex_prsr>{
		typedef complex_prsr self_t;										 // utility
		complex_prsr(){}													 // default constructor inhibits copy constructor?
		template <typename ScannerT> struct result{						 // be explicit on the result type of the parser for other classes (otherwise char*)
        	typedef typename match_result<ScannerT, std::complex<double> >::type type;
		};
		template<typename ScannerT>										 // parser function
		typename parser_result<self_t, ScannerT>::type
		parse(ScannerT const& scan) const{
			typename ScannerT::iterator_t save = scan.first;			 //save begin of parsing
			typedef typename parser_result<boost::spirit::real_parser<>, ScannerT>::type real_match_t;
			double real=0, imag=0;
			bool is_negative=false;
			double real_or_imag=0;
			real_prsr mreal_p;
			int count= (
				  (!sign_p[assign_a(is_negative)] >> (mreal_p[assign_a(real_or_imag)] >> '*' >> 'I'))
				| (mreal_p[assign_a(real)] >>  !( sign_p[assign_a(is_negative)] >> mreal_p[assign_a(imag)] >> '*' >> 'I'))
			  ).parse(scan).length();
			 
			std::complex<double> ret;
			if(real!=0){
				ret=real;
				if(is_negative) ret-=std::complex<double>(0,imag); else ret+=std::complex<double>(0,imag);
			}else{
				ret=real_or_imag;
			}
			return scan.create_match(count, ret, save, scan.first);
		}
    } complex_p;

	struct real_closure : boost::spirit::closure<real_closure, double, double>{
		member1 value;
		member2 exponent;
	};
	struct complex_closure : boost::spirit::closure<complex_closure, std::complex<double> >{
  		member1 value;
  	};
	//possible rename to minimal_head; minimal_head<double>::type <- Real
	template<typename T> 
	struct parser_map; //{};
	template<>
	struct parser_map<double>{typedef real_prsr type;};
	template<> 
	struct parser_map<std::complex<double> >{typedef complex_prsr type;};
	
	template<typename T, class Parser> struct list_prsr;

	template<typename T>
	struct parser_map<std::vector<T> >{typedef list_prsr<T, typename parser_map<T>::type> type;};

  	template<typename T, class Parser=typename parser_map<T>::type>
  	struct list_prsr : parser<list_prsr<T, Parser> >{
  		typedef list_prsr<T, Parser> self_t;
  		template <typename ScannerT> struct result{
  			typedef typename 
  				match_result<
  					ScannerT, 
  					std::vector<T> 
  				>::type type;
  		};
  		template<typename ScannerT>
  		typename parser_result<self_t, ScannerT>::type
  		parse(ScannerT const& scan) const{
			typename ScannerT::iterator_t save = scan.first;			 //save begin of parsing
			std::vector<T> ret;
			Parser item_p;
			int count=('{' >> item_p[push_back_a(ret)] >> (!(*(ch_p(',') >> item_p[push_back_a(ret)])) ) >> '}').parse(scan).length();
			return scan.create_match(count, ret, save, scan.first);  			
		}	
	};


	//import/export code starts here
	class iarchive{
		public:
		iarchive(std::istream& is) : is_(is){
			if(not is) throw std::runtime_error("stream is unreadeable");
		}
		iarchive& get(double& d){
			using namespace std;
			mathematica::real_prsr real_p;
			typedef char char_t;
			typedef multi_pass<istreambuf_iterator<char_t> > iterator_t;
			typedef skip_parser_iteration_policy<space_parser> iter_policy_t;
			typedef scanner_policies<iter_policy_t> scanner_policies_t;
			typedef scanner<iterator_t, scanner_policies_t> scanner_t;
			typedef rule<scanner_t> rule_t;
			iter_policy_t iter_policy(space_p);
			scanner_policies_t policies(iter_policy);
			iterator_t first(
				make_multi_pass(std::istreambuf_iterator<char_t>(is_)));
			scanner_t scan(
				first, make_multi_pass(std::istreambuf_iterator<char_t>()),
				policies);
			real_p[assign_a(d)].parse(scan);
			return *this;
		}
		iarchive& get(std::complex<double>& c){
			using namespace std;
			mathematica::complex_prsr complex_p;
			typedef char char_t;
			typedef multi_pass<istreambuf_iterator<char_t> > iterator_t;
			typedef skip_parser_iteration_policy<space_parser> iter_policy_t;
			typedef scanner_policies<iter_policy_t> scanner_policies_t;
			typedef scanner<iterator_t, scanner_policies_t> scanner_t;
			typedef rule<scanner_t> rule_t;
			iter_policy_t iter_policy(space_p);
			scanner_policies_t policies(iter_policy);
			iterator_t first(
				make_multi_pass(std::istreambuf_iterator<char_t>(is_)));
			scanner_t scan(
				first, make_multi_pass(std::istreambuf_iterator<char_t>()),
				policies);
			complex_p[assign_a(c)].parse(scan);
			return *this;
		}
		template<typename T>
		iarchive& get(std::vector<T>& v){
			using namespace std;
			mathematica::list_prsr<T> list_p_T;
			typedef char char_t;
			typedef multi_pass<istreambuf_iterator<char_t> > iterator_t;
			typedef skip_parser_iteration_policy<space_parser> iter_policy_t;
			typedef scanner_policies<iter_policy_t> scanner_policies_t;
			typedef scanner<iterator_t, scanner_policies_t> scanner_t;
			typedef rule<scanner_t> rule_t;
			iter_policy_t iter_policy(space_p);
			scanner_policies_t policies(iter_policy);
			iterator_t first(
				make_multi_pass(std::istreambuf_iterator<char_t>(is_)));
			scanner_t scan(
				first, make_multi_pass(std::istreambuf_iterator<char_t>()),
				policies);
			list_p_T[assign_a(v)].parse(scan);
			return *this;
		}
		std::istream& is_;
	};
	iarchive operator>>(iarchive& ia, double& d){
		return ia.get(d);
	}
	iarchive operator>>(iarchive& ia, std::complex<double>& c){
		return ia.get(c);
	}
	template<typename T>
	iarchive& operator>>(iarchive& ia, std::vector<T>& v){
		return ia.get(v);
	}
	template<typename T>
	iarchive& operator>>(iarchive& ia, typename boost::multi_array<T, 1>& a){
		std::vector<T> ret;
		ia>>ret;
		a.resize(boost::extents[ret.size()]); 							 //remember that array_refs can't be resized but only reshaped
		a.assign(ret.begin(), ret.end());
		return ia;
	}
	template<typename T>
	iarchive& operator>>(mathematica::iarchive& ia, typename boost::multi_array<T,2>& a){
		std::vector<std::vector<T> > ret;
		ia>>ret;
		a.resize(boost::extents[ret.size()][ret[0].size()]);
		for(unsigned i=0; i!=ret.size(); i++){
			if(ret[i].size()!=ret[0].size()) throw std::runtime_error("mathematica list is not square at [["+boost::lexical_cast<std::string>(i)+"]].");
			for(unsigned j=0; j!=ret[i].size(); j++){
				a[i][j]=ret[i][j];
			}
		}
		return ia;
	}
	template<typename T>
	iarchive& operator>>(mathematica::iarchive& ia, typename boost::multi_array<T,3>& a){
		std::vector<std::vector<std::vector<T> > > ret;
		ia>>ret;													 	 // why not ia>>ret? because >> is reserved for whole stream import
		a.resize(boost::extents[ret.size()][ret[0].size()][ret[0][0].size()]);
		for(unsigned i=0; i!=ret.size(); i++){
			if(ret[i].size()!=ret[0].size()) throw std::runtime_error("mathematica list is not rectangular at [["+boost::lexical_cast<std::string>(i)+"]].");
			for(unsigned j=0; j!=ret[i].size(); j++){
				if(ret[i][j].size()!=ret[i][0].size()) throw std::runtime_error(
					"mathematica list is not rectangular at [["
					+boost::lexical_cast<std::string>(i)+","+boost::lexical_cast<std::string>(j)+"]].");
				for(unsigned k=0; k!=ret[i][j].size(); k++){
					a[i][j][k]=ret[i][j][k];
				}
			}
		}
		return ia;
	}
	template<typename T>
	iarchive& operator>>(mathematica::iarchive& ia, typename boost::multi_array_ref<T,3>& a){
		std::vector<std::vector<std::vector<T> > > ret;
		ia>>ret;													 // why not ia>>ret? because >> is reserved for whole stream import
		boost::array<unsigned,3> shape={{ret.size(),ret[0].size(),ret[0][0].size()}};
		a.reshape(shape); //boost::extents[ret.size()][ret[0].size()][ret[0][0].size()]);
		for(unsigned i=0; i!=ret.size(); i++){
			if(ret[i].size()!=ret[0].size()) throw std::runtime_error("mathematica list is not rectangular at [["+boost::lexical_cast<std::string>(i)+"]].");
			for(unsigned j=0; j!=ret[i].size(); j++){
				if(ret[i][j].size()!=ret[i][0].size()) throw std::runtime_error(
					"mathematica list is not rectangular at [["
					+boost::lexical_cast<std::string>(i)+","+boost::lexical_cast<std::string>(j)+"]].");
				for(unsigned k=0; k!=ret[i][j].size(); k++){
					a[i][j][k]=ret[i][j][k];
					clog<<"assigned on "<<i<<' '<<j<<' '<<k<<endl;
				}
			}
		}
		return ia;
	}
	
	class oarchive{
		public:
		oarchive(std::ostream& os) : os_(os){
			if(not os) throw std::runtime_error("output stream is unwritable");
		}
		//private:
		friend oarchive& operator<<(oarchive&, double const& d);
		oarchive& put(std::string const& s){
			os_<<"\""<<s<<"\"";
			return *this;
		}
		oarchive& put(double const& d){
			using namespace std; 
			using boost::replace_last;
			if(d!=d){os_<<"NaN"; return *this;} //detects NaN, needs a mathematica package to make sense out of it;
			ostringstream oss; oss<<setprecision(15)<<d; string s=oss.str();
			replace_last(s, "e", "*^");         // doing this from scratch (with no strings is a mess)
			replace_last(s, "^+", "^");
			os_<<s; 
			return *this;
		}
		oarchive& put(std::complex<double> const& d){
			put(real(d)); if(not (imag(d)<=0)) os_<<"+";                        //otherwise will always automatically print sign
			if(imag(d)!=0){put(imag(d)); os_<<"*I";}
			return *this;}
		public:
		template<class T> oarchive& put(vector<T> const& v);
		template<typename T, unsigned D> oarchive& put(multi_array<T,D> const& ma){return shift_array(ma);}
		template<typename T, unsigned D> oarchive& put(multi_array_ref<T,D> const& mar){return shift_array(mar);}
		template<typename Element>
		struct shift_array_dispatch{
			template<class MultiArray>
			static oarchive& go(oarchive& mo, MultiArray const& ma){
				mo.os_<<"{";
				for(typename MultiArray::const_iterator i=ma.begin(); i!=ma.end(); i++){
					go(mo,*i);
					if(i+1!=ma.end()) mo.os_<<",";
				}
				mo.os_<<"}";
				return mo;
			}
			static oarchive& go(oarchive& mo,Element const& e){
				return mo.put(e);
			}
		};
		template<typename MultiArray> 
		oarchive& shift_array(MultiArray const& ma){return shift_array_dispatch<typename MultiArray::element>::go(*this,ma);}
		//private:
		std::ostream& os_;
	};
	template<typename T>
	oarchive& oarchive::put(std::vector<T> const& v){
		os_<<"{";
		for(typename std::vector<T>::const_iterator i=v.begin(); i!=v.end(); ++i){
			if(i!=v.begin()) os_<<", ";
			this->put(*i); 
		}
		os_<<"}";
		return *this;
	}

	//todo: make a template operator<< to capture this;
	oarchive& operator<<(oarchive& oa, std::string const& s){
		oa.put(s); oa.os_<<'\n'; return oa;
	}
	oarchive& operator<<(oarchive& oa, double const& d){
		oa.put(d); oa.os_<<'\n'; return oa;
	}
	template<typename T>
	oarchive& operator<<(oarchive& oa, std::vector<T> const& v){
		oa.put(v); oa.os_<<'\n'; return oa;
	}	
	template<typename T, unsigned D>
	oarchive& operator<<(oarchive& oa, boost::multi_array_ref<T,D> const& mar){
		oa.put(mar); oa.os_<<'\n'; return oa;
	}
	template<typename T, unsigned D>
	oarchive& operator<<(oarchive& oa, typename boost::detail::multi_array::multi_array_view<T,D> const& mav){
		oa.shift_array(mav); oa.os_<<'\n'; return oa;
	}
	template<typename T, unsigned D>
	oarchive& operator<<(oarchive& oa, typename boost::detail::multi_array::sub_array<T,D> const& sa){
		oa.shift_array(sa); oa.os_<<'\n'; return oa;
	}
	oarchive& operator<<(oarchive& oa, std::complex<double> const& c){
		oa.put(c); oa.os_<<'\n'; return oa;
	}
	template<class T>
	oarchive& operator<<(oarchive& oa, boost::serialization::nvp<T> const& nvp){
		oa.os_<<nvp.name()<<"=";
		oa.put(nvp.value());
		oa.os_<<";\n";
		return oa;
	}
	
	//basically syntax sugar
	class archive{
		public:
		archive(boost::filesystem::path const& p) : p_(p){}
		private:
		boost::filesystem::path const& p_;
		template<typename T> friend void operator<(archive const& ar, T const& t);
	};
	template<typename T>
	void operator<(archive const& ar, T const& t){
		boost::filesystem::ofstream os(ar.p_);
		mathematica::oarchive oa(os);
		oa<<t;
		return;
	}
}
#endif
