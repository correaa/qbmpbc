#ifdef COMPILE_INSTRUCTIONS
	ln -sf $0 $0.cpp && c++ -D_FORMAT_CAST_TEST $0.cpp -o $0.cpp.x && ./$0.cpp.x && rm -f $0.cpp.x
	exit
#endif
#ifndef FORMAT_CAST_HPP
#define FORMAT_CAST_HPP
#include<boost/format.hpp>
#include<boost/lexical_cast.hpp>
namespace boost{

template<class Target, class Source>
Target format_cast(std::string const& strFrmt, Source const& s){
  Target t;
  std::stringstream ss;
  ss << boost::format("%1$"+strFrmt) % s;
  ss >> t;
  if(not ss) throw boost::bad_lexical_cast();
  return t;
}
template<class Source>
std::string format_cast(std::string const& strFrmt, Source const& s){
  std::string t;
  std::stringstream ss;
  ss << boost::format("%1$"+strFrmt) % s;
	t = ss.str();
  if(not ss) throw boost::bad_lexical_cast();
  return t;
}
template<class Target>
struct myformat{
	template<class Source>
	static Target cast(std::string const& strFrmt, Source const& s){
		Target t;
		std::stringstream ss;
		ss << boost::format("%1$"+strFrmt) % s;
		ss >> t;
		if(not ss) throw boost::bad_lexical_cast();
		return t;
	}
};
template<>
struct myformat<std::string>{
	template<class Source>
	static std::string cast(std::string const& strFrmt, Source const& s){
		std::string t;
		std::stringstream ss;
		ss << boost::format("%1$"+strFrmt) % s;
		//ss >> t;
		t = ss.str();
		if(not ss) throw boost::bad_lexical_cast();
		return t;
	}
};

/*
template<class Source>
std::string format_cast<std::string, Source>(std::string const& strFrmt, Source const& s){
  std::string t;
  std::stringstream ss;
  ss << boost::format("%1$"+strFrmt) % s;
  t = ss.str();
  if(not ss) throw boost::bad_lexical_cast();
  return t;
}*/

}
//use as
//unsigned u = 4;
//std::string str = boost::format_cast<std::string>("02f", u); // str == 04;
#ifdef _FORMAT_CAST_TEST
#include<iostream>
#include<boost/units/systems/si.hpp>
#include<boost/units/io.hpp>

int main(){
	{
		unsigned u = 4;
		std::string str = boost::format_cast<std::string>("03u", u); // str == 04;
		std::cout<<str<<std::endl;
	}
	{
		double d = 234.567;
		std::string str = boost::format_cast<std::string>("+08.2f", d); // str == 04;
		std::cout<<str<<std::endl;
	}
	{
		using namespace boost::units;
		quantity<si::length> l = 14.56*si::meter;
		std::string str =
			boost::format_cast("011.2f", l);
			//boost::myformat<std::string>::cast("09.2f", l); // str == 04;
		std::cout<<str<<std::endl;
	}
	return 0;
}
#endif
#endif

