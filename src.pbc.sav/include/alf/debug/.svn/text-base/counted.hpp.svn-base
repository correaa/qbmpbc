// counted object class 
// based on http://www.linuxtopia.org/online_books/programming_books/c++_practical_programming/c++_practical_programming_129.html

#include<iostream>
#include<typeinfo>
namespace debug{
	
//template<class Derived> class counted;	
	
template<class Derived>													// designed for class myclass : public counted<myclass>{};
class counted{
	public:
	static int count;
	protected:
	counted();
	counted(const counted<Derived>&);
	~counted();
	static int get_count(){return count;}
};
template<class Derived> int counted<Derived>::count = 0;				// may not be necessary, static ints are init to zero

template<class Derived>
counted<Derived>::counted(){
	++count;
	std::clog<<"counted<>::counted object of type "<<typeid(Derived).name()<<" constructed, total in use "<<count<<std::endl;
}

template<class Derived>
counted<Derived>::counted(const counted<Derived>&){
	std::clog<<"counted<>::counted object of type "<<typeid(Derived).name()<<" copy-constructed, total in use "<<count<<std::endl;
	++count;
}

template<class Derived>
counted<Derived>::~counted(){
	--count;
	std::clog<<"counted<>::~counted object of type "<<typeid(Derived).name()<<" destructed, total in use "<<count<<std::endl;
}

}
