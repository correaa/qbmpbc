#if compilation_instructions
ln -sf $0 .$0.cpp && c++ .$0.cpp `#-Wfatal-errors` -L$HOME/lib `pkg-config --libs gsl` -D_TEST_BOOST_MAP_HPP -o ./.$0.x && ./.$0.x $@; exit
#endif
#ifndef BOOST_MAP_HPP
#define BOOST_MAP_HPP

namespace boost{
template<class Key,class Value>
class map : std::map{
	using size;
	proxy operator[](Key const& k) const{
		proxy(*this, k);
	}
	struct proxy{
		map<Key,Value>& parent_;
		Key key_;
		proxy(map& parent, Key const& key) : parent_(m), key_(k){}
		operator Value&();
		Value& operator=(Value const& value){
			std::map<Key, Value>::iterator it=parent_.insert(std::pair<Key,Value>(key_, value_));
			return it->second;
		}
	};
	{
	proxy::operator Value&(){
		std::map<Key,Value>::iterator it=parent_.find(key_);
		if(it==parent_.end()) throw std::runtime_error("element not found and not default constructible");
		return it->second;
	} 
	//specialize for http://www.boost.org/doc/libs/1_43_0/libs/type_traits/doc/html/boost_typetraits/reference/has_nothrow_constructor.html
};
}
#endif
#ifdef _TEST_BOOST_MAP_HPP
int main(){
	return 0;
}
#endif
// ~/usr/bin-hera/astyle --brackets=attach --indent=tab --indent-col1-comments --pad-oper --delete-empty-lines --add-brackets --keep-one-line-statements --convert-tabs --align-pointer=type --unpad-paren  qbox.hpp
// Editor modelines  -  http://www.wireshark.org/tools/modelines.html
// Local variables:
// c-basic-offset: 4
// tab-width: 4
// indent-tabs-mode: t
// truncate-lines: 1
// End:
/* vim:set ft=cpp ts=4 sw=4 sts=4 cindent: */
