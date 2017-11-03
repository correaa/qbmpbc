#ifdef compile_instructions
ln -sf $0 $0.cpp && c++ -laspell -lhunspell $0.cpp -Wall -Wextra  -o $0.x  -D_ASPELL_TEST && ./$0.x $@ 
rm -rf $0.cpp $0.x
exit
#endif

#include<string>
#include<list>
#include<stdexcept>
namespace aspell{
#include<aspell.h>

class speller;
class config{
	AspellConfig * pimpl_;
	public:
	config() : pimpl_(new_aspell_config()){}
	config(config const& other) : pimpl_(aspell_config_clone(other.pimpl_)){}
	config& replace(std::string const& option = "lang",  std::string const& value = "en_US"){
	     aspell_config_replace(pimpl_, option.c_str(), value.c_str());
		return *this;
	}
	~config(){delete_aspell_config(pimpl_);}
	friend class speller;
};

class speller{
	AspellSpeller * pimpl_;
	public:
	speller(config const& c) : pimpl_(0){
		AspellCanHaveError * possible_err = new_aspell_speller(c.pimpl_);
		if(aspell_error_number(possible_err) != 0){
			throw std::runtime_error("cannot construct aspell::speller because " + std::string(aspell_error_message(possible_err)));
		}else{
			pimpl_ = to_aspell_speller(possible_err);
		}
	}
	int check(std::string const& word) const{return aspell_speller_check(pimpl_, word.c_str(), -1 /*word.size()*/);}
	// this two funcions are not very useful in a programmatic context
	//add  store_repl based on aspell_speller_store_repl(spell_checker, misspelled_word, size, correctly_spelled_word, size);
	//add  add_to_session|personal based on aspell_speller_add_to_session|personal(spell_checker, word, size);
	~speller(){delete_aspell_speller(pimpl_);}
	// universal interface
	speller(std::string const& dictionary = "en_US") : pimpl_(0){
		config c; c.replace("lang", dictionary);
		AspellCanHaveError * possible_err = new_aspell_speller(c.pimpl_);
		if(aspell_error_number(possible_err) != 0){
			throw std::runtime_error("cannot construct aspell::speller because " + std::string(aspell_error_message(possible_err)));
		}else{
			pimpl_ = to_aspell_speller(possible_err);
		}
	}
	bool operator()(std::string const& word) const{return check(word);}
	std::list<std::string> suggest(std::string const& word) const{
		std::list<std::string> ret;
		AspellWordList const * suggestions = aspell_speller_suggest(pimpl_, word.c_str(), -1 /*word.size()*/);
		AspellStringEnumeration * elements = aspell_word_list_elements(suggestions);
		const char * suggested_word;
		while ( (suggested_word = aspell_string_enumeration_next(elements)) != NULL ){
			ret.push_back(suggested_word);
		}
		delete_aspell_string_enumeration(elements);
		//suggestions pointer is not deleted
		return ret;
	}
};

}

#ifdef _ASPELL_TEST
#include<iostream>
#include<vector>
#include "hunspell.hpp"
int main(){
	{
		//aspell::config c;
		aspell::speller s("en_US"); //hunspell::speller s; //(c); 
		std::string word = "house";
		std::cout << s(word) << std::endl;
		if(not s(word))
		{
			std::list<std::string> l = s.suggest(word);
			for(std::list<std::string>::const_iterator it = l.begin(); it!=l.end(); ++it){
				std::cout << *it << std::endl;
			}
		}
	}

	return 0;
}
#endif

