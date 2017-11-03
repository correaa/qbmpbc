#ifdef compile_instructions
ln -sf $0 $0.cpp && c++ -I$HOME/prj -I$HOME/usr/include -lboost_system -lboost_regex -laspell $0.cpp -o ./$0.x -D_TEST_LATEX_VECTOR_IO_HPP && ./$0.x $@
rm -f .$0.x .$0.cpp
exit
#endif

#include "../latex.hpp"

namespace latex{
	template<class T>
	ostream& operator<<(ostream& os, std::vector<T> const& v){
		os << "\\{";
		for(typename std::vector<T>::const_iterator it = v.begin(); it != v.end(); ++it){
		//	if(distance(v.begin(), it)>10){
		//		it=boost::prior(v.end());
		//		os << "..., " << *it;
		//		break;
		//	}
			os << *it;
			if(boost::next(it)!=v.end()) os << ", ";
		}
		os << "\\}";
		return os;
	}
	template<class T1, class T2>
	ostream& operator<<(ostream& os, std::pair<T1, T2> const& p){
		return os << "(" << p.first << ", " << p.second << ")";
	}
	template<class T1, class T2>
	ostream& operator<<(ostream& os, std::map<T1, T2> const& v){
		os << "\\{";
		for(typename std::map<T1, T2>::const_iterator it = v.begin(); it != v.end(); ++it){
			if(distance(v.begin(), it)>10){
				it=boost::prior(v.end());
				os << "..., " << *it;
				break;
			}
			os << "(" << it->first << "$\\rightarrow$" << it->second << ")";
			if(boost::next(it)!=v.end()) os << ", ";
		}
		os << "\\}";
		return os;
	}
}

#ifdef _TEST_LATEX_VECTOR_IO_HPP
#include<boost/assign.hpp>
int main(){
	std::vector<double> data = boost::assign::list_of<double>(5.)(6.)(7.)(8.);
	std::vector<std::pair<double, double> > data2 = boost::assign::list_of<std::pair<double, double> >(5.,1.)(6.,2.)(7.,2.)(8.,1.)(8.,1.)(8.,1.)(8.,1.);
	std::map<double,double> next = boost::assign::map_list_of(1,2)(2,3)(3,4)(4,5)(5,6);

	latex::ostream lout("vector_io.pdf");
	lout << "This is a vector " << data2 << latex::par;
	lout << "This is a map " << next;
	return 0;
}
#endif
// vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent:

