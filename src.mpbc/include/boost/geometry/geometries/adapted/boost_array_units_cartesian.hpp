// Boost.Geometry (aka GGL, Generic Geometry Library)
//
// Copyright Alfredo Correa 2010
// Copyright Bruno Lalande 2008, 2009
// Copyright Barend Gehrels 2007-2009, Geodan, Amsterdam, the Netherlands.
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_GEOMETRIES_ADAPTED_BOOST_ARRAY_UNITS_CARTESIAN_HPP
#define BOOST_GEOMETRY_GEOMETRIES_ADAPTED_BOOST_ARRAY_UNITS_CARTESIAN_HPP

#ifdef BOOST_GEOMETRY_ADAPTED_BOOST_ARRAY_UNITS_COORDINATE_SYSTEM_DEFINED
#error Include only one headerfile to register coordinate coordinate_system for adapted boost array
#endif

#define BOOST_GEOMETRY_ADAPTED_BOOST_ARRAY_UNITS_COORDINATE_SYSTEM_DEFINED


#include <boost/geometry/geometries/adapted/boost_array_units.hpp>
#include <boost/geometry/units.hpp>

namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_TRAITS_SPECIALIZATIONS
namespace traits
{
    template <class Unit, std::size_t N>
    struct coordinate_system<boost::array<boost::units::quantity<Unit>, N> >
    { typedef cs::units_cartesian<Unit> type; };

}
#endif


}} // namespace boost::geometry


#endif // BOOST_GEOMETRY_GEOMETRIES_ADAPTED_BOOST_ARRAY_UNITS_CARTESIAN_HPP

