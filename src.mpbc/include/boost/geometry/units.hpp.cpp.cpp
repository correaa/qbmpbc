#ifdef compile_instructions 
ln -f $0 $0.cpp && c++ -std=c++0x -D_TEST_BOOST_GEOMETRY_UNITS_HPP -I$HOME/prj -I$HOME/usr/include -Wno-unused-variable `pkg-config --libs gsl` -lboost_filesystem $0.cpp -o $0.x -Wall && \
./$0.x $@; exit
#endif
#ifndef BOOST_GEOMETRY_UNITS_HPP
#define BOOST_GEOMETRY_UNITS_HPP
// units/tags.hpp, similar to core/tags.hpp
namespace boost{namespace geometry{
	template<class Unit>
	struct units_cartesian_tag {};
}}

// units/cs.hpp, similar to core/cs.hpp
#include<boost/geometry/core/cs.hpp>
namespace boost{namespace geometry{
namespace cs{
	template<class Unit>
	struct units_cartesian{};
}
namespace traits{
	template<class Unit>
	struct cs_tag<cs::units_cartesian<Unit> >{
		typedef units_cartesian_tag<Unit> type;
	};
}
}}
// similar to 
#include<boost/geometry/strategies/distance.hpp>
#include<boost/geometry/geometries/point.hpp>
#include<boost/units/systems/si.hpp>

#include<boost/geometry/geometry.hpp>
namespace boost { namespace geometry{
namespace model{
template<class T, std::size_t Dimension, class Unit>
struct point<T, Dimension, cs::units_cartesian<Unit> > : protected model::point<T, Dimension, cs::cartesian >{
};
template<class T, class Unit>
struct point<T, 2, cs::units_cartesian<Unit> > : protected model::point<T, 2, cs::cartesian >{
	point(
		boost::units::quantity<Unit, T> const& x, 
		boost::units::quantity<Unit, T> const& y
	) : model::point<T, 2, cs::cartesian >(x.value(), y.value()){}
};
template<class T, class Unit>
struct point<T, 3, cs::units_cartesian<Unit> > : protected model::point<T, 3, cs::cartesian >{
	point(
		boost::units::quantity<Unit, T> const& x, 
		boost::units::quantity<Unit, T> const& y, 
		boost::units::quantity<Unit, T> const& z
	) : model::point<T, 3, cs::cartesian >(x.value(), y.value(), z.value() ){}
};
}

template<class T, std::size_t Dimension, class Unit>
boost::units::quantity<Unit, T> distance(
	model::point<T, Dimension, cs::units_cartesian<Unit> > const& p1,
	model::point<T, Dimension, cs::units_cartesian<Unit> > const& p2
){
	return boost::units::quantity<Unit, T>::from_value(
		boost::geometry::distance(
			(model::point<T, Dimension, cs::cartesian > const&)p1,
			(model::point<T, Dimension, cs::cartesian > const&)p2
		)
	);
}

namespace strategy{namespace distance{namespace services{

template<class Unit> struct return_type<boost::geometry::units_cartesian_tag<Unit> >{
	typedef boost::units::quantity<Unit> type;
};

template<std::size_t Dimension, class Unit>
struct default_strategy<
	boost::geometry::point_tag, 
	boost::geometry::model::point<
		double, 
		Dimension, 
		boost::geometry::cs::units_cartesian<Unit> 
	>, 
	boost::geometry::model::point<
		double, 
		Dimension, 
		boost::geometry::cs::units_cartesian<Unit> 
	>,
	boost::geometry::units_cartesian_tag<Unit>, 
	boost::geometry::units_cartesian_tag<Unit>, 
	void
>{
	typedef boost::geometry::units_cartesian_tag<Unit> type;
};
}}}

}}

#endif
#ifdef _TEST_BOOST_GEOMETRY_UNITS_HPP
#include<boost/geometry/geometry.hpp>
#include<boost/units/systems/si.hpp>
using namespace boost::geometry;
using namespace boost::units;
int main(){
	model::point<double, 3, cs::cartesian > p1(1.,2.,3.);
	model::point<double, 3, cs::cartesian > p2(2.,3.,4.);
	std::cout << distance(p1,p2) << std::endl;

	model::point<double, 3, cs::units_cartesian<si::length> > p1u(1.*si::meter,2.*si::meter,3.*si::meter);
	model::point<double, 3, cs::units_cartesian<si::length> > p2u(2.*si::meter,3.*si::meter,4.*si::meter);


	std::cout << distance(p1u, p2u) << std::endl;

	return 0;
}
#endif
// vim:set ft=cpp ts=4 sw=4 sts=4 nowrap: cindent:
