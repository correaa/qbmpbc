#include <limits>
#include <iostream>

namespace fftw3{ //allocators based on http://www.josuttis.com/libbook/memory/myalloc.hpp.html
   template <class T>
   class allocator{														// example usage std::vector<double, fftw3::allocator<double> > vfft(100);
     public:
       typedef T        value_type;										// type definitions
       typedef T*       pointer;
       typedef const T* const_pointer;
       typedef T&       reference;
       typedef const T& const_reference;
       typedef std::size_t    size_type;
       typedef std::ptrdiff_t difference_type;
       template <class U>												// rebind allocator to type U
       struct rebind {
           typedef allocator<U> other;
       };
       pointer address (reference value) const {						// return address of values
           return &value;
       }
       const_pointer address (const_reference value) const {
           return &value;
       }
       allocator() throw() {}											// constructors and destructor, nothing to do, because the allocator has no state
       allocator(const allocator&) throw() {}
       template <class U>
         allocator (const allocator<U>&) throw(){}
       ~allocator() throw(){}
       size_type max_size () const throw() {							// return maximum number of elements that can be allocated
           return std::numeric_limits<std::size_t>::max() / sizeof(T);
       }
       pointer allocate (size_type num, const void* = 0){				// allocate but don't initialize num elements of type T
           // print message and allocate memory with global new
           std::cerr << "allocate " << num << " element(s)"
                     << " of size " << sizeof(T) << std::endl;
           pointer ret = (pointer) fftw_malloc(sizeof(T)*num);			// original code: pointer ret = (pointer)(::operator new(num*sizeof(T)));
           std::cerr << " allocated at: " << (void*)ret << std::endl;
           return ret;
       }
       void construct (pointer p, const T& value){						// initialize elements of allocated storage p with value value
           new((void*)p)T(value);										// initialize memory with placement new
       }
       void destroy (pointer p) {										// destroy elements of initialized storage p
           p->~T();														// destroy objects by calling their destructor
       }
       void deallocate (pointer p, size_type num) {						// deallocate storage p of deleted elements
           std::cerr << "deallocate " << num << " element(s)"			// print message and deallocate memory with global delete
                     << " of size " << sizeof(T)
                     << " at: " << (void*)p << std::endl;
           fftw_free((void*)p);											//original: ::operator delete((void*)p);
       }
   }; // end of class allocator
   template <class T1, class T2>										// return that all specializations of this allocator are interchangeable
   bool operator== (const allocator<T1>&,
                    const allocator<T2>&) throw() {
       return true;														// all allocators of this metatype are equal 
   }
   template <class T1, class T2>
   bool operator!= (const allocator<T1>&,
                    const allocator<T2>&) throw() {
       return false;
   }
} // end of namespace fftw3
