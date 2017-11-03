// counted object class 
// based on http://www.linuxtopia.org/online_books/programming_books/c++_practical_programming/c++_practical_programming_129.html

template<class Derived>													// designed for class myclass : public counted<myclass>{};
class counted{
	static int count;
	protected:
	counted(){++count;}
	counted(const counted<Derived>&){++count;} 								// not sure why I would need this
	~counted(){--count;}
	static int get_count(){return count;}
};
template<class Derived> int counted<Derived>::count = 0;				// may not be necessary, static ints are init to zero
 
