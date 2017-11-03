#ifdef compile_instructions
ln -sf $0 $0.cpp && c++ -lhunspell -lboost_regex $0.cpp -Wall -Wextra  -o $0.x  -D_HUNSPELL_TEST && ./$0.x $@ 
rm -rf $0.cpp $0.x
exit
#endif
#ifndef HUNSPELL_HPP
#define HUNSPELL_HPP
#include<list> //#include<vector>
#include<string>
#include<hunspell/hunspell.hxx>
namespace hunspell{
class speller{
	mutable Hunspell impl;
	public:
	speller(std::string const& dictionary = "en_US") : impl( ("/usr/share/hunspell/" + dictionary+".aff").c_str(),  ("/usr/share/hunspell/" + dictionary+".dic").c_str()){
	}
	int spell(std::string const& word) const{return impl.spell(word.c_str());}
	std::list<std::string> suggest(std::string const& word) const{
		char ** wlst;
		std::list<std::string> ret; //std::vector<std::string> ret(impl.suggest(&wlst,word.c_str()));
		unsigned size = impl.suggest(&wlst,word.c_str());
		for (unsigned i=0; i < size; i++){
			ret.push_back(wlst[i]); //ret[i] = wlst[i];
		}
		impl.free_list(&wlst, ret.size());
		return ret;
	}
	bool operator()(std::string const& word) const{
		return spell(word);
	}
};
}
#endif

#ifdef _HUNSPELL_TEST
#include<iostream>
#include<boost/tokenizer.hpp>
#include<boost/regex.hpp>
int main(){
	{
		hunspell::speller c;
		std::string word = "chouse";
		std::cout << c(word) << std::endl;   
		if(not c(word)){
			std::list<std::string> const l = c.suggest(word);
			for(std::list<std::string>::const_iterator it = l.begin(); it!=l.end(); ++it){
				std::cout << *it << std::endl;
			}
		}
	}
	{
		std::string sentence = "Shakespeare\\cite{caca} was born and \\emph{raised} in Stratford-upon-Avon.";
		typedef boost::tokenizer<boost::char_separator<char> > tokenizer; 
		 boost::char_separator<char> sep(" .\\"); 
		tokenizer tok(sentence, sep); 
		hunspell::speller c;
		for(tokenizer::iterator it = tok.begin(); it != tok.end(); ++it){
		    std::cout << *it <<" ";
			if(it->size()>2 and boost::regex_match(*it, boost::regex("^[a-z]*$") ) ){
				if(c(*it)){
					std::cout << "ok";
				}else{
					std::cout << "not ok suggestions: ";
					std::list<std::string> const s = c.suggest(*it);
					for(std::list<std::string>::const_iterator it2=s.begin(); it2 != s.end() ; ++it){
						std::cout << *it << " ";
					}
					std::cout << ".";
				}
			}else{
				std::cout << "is not a word";
			}
			std::cout<< std::endl;
		}
	}
	return 0;
}
#endif
/* vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent: */

