#ifdef compile_instructions 
ln -sf $0 .$0.cpp && c++ -std=c++0x -D_BOOST_IOMANIP_HPP_TEST -I$HOME/prj `pkg-config --libs gsl` -lboost_system .$0.cpp -o .$0.x -Wall && ./.$0.x $@; exit
#endif
/* The following code example is taken from the book
 * "The C++ Standard Library - A Tutorial and Reference"
 * by Nicolai M. Josuttis, Addison-Wesley, 1999
 *
 * (C) Copyright Nicolai M. Josuttis 1999.
 * Permission to copy, use, modify, sell and distribute this software
 * is granted provided this copyright notice appears in all copies.
 * This software is provided "as is" without express or implied
 * warranty, and with no claim as to its suitability for any purpose.
 */
#include <iostream>
//#include <vector>
//#include <set>
//#include <algorithm>
#include <limits>
//#include<iomanip>

namespace boost{
template <class charT, class traits>
inline
std::basic_istream<charT,traits>&
ignore_line (
	std::basic_istream<charT,traits>& strm
)
{
    // skip until end-of-line
    strm.ignore(std::numeric_limits<int>::max(),strm.widen('\n'));

    // return stream for concatenation
    return strm;
}
 // use as std::cin >> boost::ignore<double> >> boost::ignore<int>
template<class T 
	,class charT, class traits
>
inline
std::basic_istream<charT, traits>&
ignore(
	std::basic_ifstream<charT, traits>& strm
){
  T ignored;
  strm >> ignored;
  return strm;
}

template<class charT, class traits>
inline
std::basic_istream<charT, traits>&
ignore(std::basic_istream<charT, traits>& strm){
  std::string ignored;
  strm >> ignored;
  return strm;
}


}

#ifdef _BOOST_IOMANIP_HPP_TEST
int main()
{
    int i;
    std::cout << "read int and ignore rest of the line" << std::endl;
    std::cin >> i;

    // ignore the rest of the line
    std::cin >> boost::ignore_line; //ignoreLine;

    std::cout << "int: " << i << std::endl;

    std::cout << "read int and ignore two lines" << std::endl;
    std::cin >> i;
		std::cin >> boost::ignore;

    // ignore two lines
    std::cin >> boost::ignore_line >> boost::ignore_line; //ignoreLine >> ignoreLine;

    std::cout << "int: " << i << std::endl;
}
#endif
