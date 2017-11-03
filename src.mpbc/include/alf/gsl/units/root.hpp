#ifdef compile_instructions
ln -sf $0 .$0.cpp && c++ -std=c++0x -I$HOME/prj .$0.cpp -Wall `#-Wfatal-errors` -L$HOME/lib `pkg-config --libs gsl` -D_TEST_ROOT_BOOST_UNITS_HPP -o ./.$0.x && ./.$0.x $@
rm -f .$0.x .$0.cpp
exit
#endif
// use as: #include "gsl/units/root.hpp"

#ifndef ROOT_BOOST_UNITS_HPP
#define ROOT_BOOST_UNITS_HPP
#include "../../gsl/root.hpp"
#include "boost/units/interval.hpp"
#include<boost/units/systems/si.hpp>
#include "../../gsl/units/adimensional.hpp"
#include<boost/numeric/interval/io.hpp>
#include<iostream>

namespace gsl{
/// Numerical root finders in one dimension with no derivatives http://www.gnu.org/software/gsl/manual/html_node/One-dimensional-Root_002dFinding.html
namespace root{
	namespace units{
		using namespace boost::units;
		template<class UnitDomain, class UnitFunction>
		class solver : 
			protected gsl::root::solver
		{
			//typedef boost::numeric::interval<quantity<UnitDomain> > interval;
			public:
			solver(
				boost::function<quantity<UnitFunction>(quantity<UnitDomain> const&)> const& f,
				boost::numeric::interval<quantity<UnitDomain> > const& iv
			) : 
			gsl::root::solver(
				boost::units::adimensionalize_function<UnitDomain, UnitFunction>(f), 
				boost::numeric::interval<double>(lower(iv).value(), upper(iv).value() ) 
			){
			}
			using gsl::root::solver::operator++;
			using gsl::root::solver::iterate;
			quantity<UnitDomain> root() const{return quantity<UnitDomain>::from_value(gsl::root::solver::root());}
			boost::numeric::interval<quantity<UnitDomain> > x() const{
				return boost::numeric::interval<quantity<UnitDomain> >(
					quantity<UnitDomain>::from_value(gsl::root::solver::x_lower()),
					quantity<UnitDomain>::from_value(gsl::root::solver::x_upper())
				);
			}
		};
	}
}
}
#endif//ROOT_BOOST_UNITS_HPP

#ifdef _TEST_ROOT_BOOST_UNITS_HPP
#include<boost/spirit/home/phoenix.hpp>
#include<boost/units/cmath.hpp>
// /function/function.hpp>

using namespace boost::units;
quantity<si::area> G(quantity<si::length> const& l){
	return l*l-2.*si::meter*si::meter;
	//return l-2.*si::meter;
}
struct Gs_impl{
    template <typename Arg>
    struct result{
        typedef quantity<si::area> type;
    };
	template <typename Arg>
	typename result<Arg>::type operator()(Arg l) const{
		return G(l);
    }
};


using namespace boost::units;
quantity<si::area> G2(quantity<si::length> const& l, quantity<si::area> const& a){
	return l*l-a;
}

struct G2s{
	int i;
	G2s(int i) : i(i){}
	quantity<si::area> operator()(quantity<si::length> const& l, quantity<si::area> const& a) const{
		return l*l-a;
	}
};
struct G2s_impl{
	G2s g2;
	G2s_impl(G2s g2) : g2(g2){}
	template<typename Arg1, typename Arg2>
	struct result{
		typedef quantity<si::area> type;
	};
	template <typename Arg1, typename Arg2>
	typename result<Arg1, Arg2>::type operator()(Arg1 l, Arg2 a) const{
		return g2(l, a);
	}
};
struct G2s_smart : boost::phoenix::function<G2s_smart const&>{
	int i;
	G2s_smart(int i) : boost::phoenix::function<G2s_smart const&>(*this), i(i){}
	quantity<si::area> operator()(quantity<si::length> const& l, quantity<si::area> const& a) const{
		return l*l-a;
	}
	using boost::phoenix::function<G2s_smart const&>::operator();
	private:
	template <typename Arg1, typename Arg2>
	struct result{
		typedef quantity<si::area> type;
	};
	template <typename Arg1, typename Arg2>
	typename result<Arg1, Arg2>::type operator()(Arg1 l, Arg2 a) const{
		return this->operator()((quantity<si::length> const&)l, (quantity<si::area> const&)a);
	}
};

struct c_function_G2{
	template<typename Arg1, typename Arg2>
	struct result{
		typedef quantity<si::area> type;
	};
	template <typename Arg1, typename Arg2>
	typename result<Arg1, Arg2>::type operator()(Arg1 l, Arg2 a) const{
		return G2(l, a);
    }
};
template<quantity<si::area>(f)(quantity<si::length> const&, quantity<si::area> const&)>
struct c_function{
    template<typename Arg1, typename Arg2> struct result{typedef quantity<si::area> type;};
	template <typename Arg1, typename Arg2>
	typename result<Arg1, Arg2>::type operator()(Arg1 l, Arg2 a) const{
		return f(l, a);
    }
};

int main(){
	using std::cout; using std::endl;
	{
		boost::function<quantity<si::area>(quantity<si::length> const&)> F(&G);
		cout << F(1.*si::meter) <<endl;
		boost::numeric::interval<quantity<si::length> > iv(0.*si::meter, 5.*si::meter);
		gsl::root::units::solver<si::length, si::area> sr(F, iv);
		for(; width(sr.x()) > 0.01*si::meter ; ++sr){
			++sr;
		}
		cout << sr.x() << endl;
	}
	{
		using namespace boost::phoenix; using namespace arg_names;
		function<Gs_impl> Gf;
		quantity<si::length> l = 1.*si::meter;
		cout << Gf(arg1)(l) << endl;
		boost::function<quantity<si::area>(quantity<si::length> const&)> F=Gf(arg1);
		boost::numeric::interval<quantity<si::length> > iv(0.*si::meter, 5.*si::meter);
		gsl::root::units::solver<si::length, si::area> sr(F, iv);
		for(; width(sr.x()) > 0.01*si::meter ; ++sr){
			++sr;
		}
		cout << sr.x() << endl;
	}
	{
		using namespace boost::phoenix; using namespace arg_names;
		function<c_function<G2> > G2f;
		quantity<si::length> l = 1.*si::meter;
		quantity<si::area> a = 2.*pow<2>(si::meter);
		cout << G2f(arg1, 2.*pow<2>(si::meter) )(l) << endl;
		boost::function<quantity<si::area>(quantity<si::length> const&)> F(G2f(arg1, a));
		boost::numeric::interval<quantity<si::length> > iv(0.*si::meter, 5.*si::meter);
		gsl::root::units::solver<si::length, si::area> sr(F, iv);
		for(; width(sr.x()) > 0.01*si::meter ; ++sr){
			++sr;
		}
		cout << sr.x() << endl;
	}
	{
		using namespace boost::phoenix; using namespace arg_names;
		function<G2s_impl> G2f(G2s(2));
		quantity<si::length> l = 1.*si::meter;
		quantity<si::area> a = 2.*pow<2>(si::meter);
		cout << G2f(arg1, 2.*pow<2>(si::meter) )(l) << endl;
		boost::function<quantity<si::area>(quantity<si::length> const&)> F(G2f(arg1, a));
		boost::numeric::interval<quantity<si::length> > iv(0.*si::meter, 5.*si::meter);
		gsl::root::units::solver<si::length, si::area> sr(F, iv);
		for(; width(sr.x()) > 0.01*si::meter ; ++sr){
			++sr;
		}
		cout << sr.x() << endl;
	}
	return 0;
}
#endif //_TEST_ROOT_BOOST_UNITS_HPP

