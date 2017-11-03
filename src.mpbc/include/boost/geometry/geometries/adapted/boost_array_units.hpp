// Boost.Geometry (aka GGL, Generic Geometry Library)
//
// Copyright Alfredo Correa 2010
// Copyright Bruno Lalande 2008, 2009
// Copyright Barend Gehrels 2007-2009, Geodan, Amsterdam, the Netherlands.
// Use, modification and distribution is subject to the Boost Software License,
// Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_GEOMETRY_GEOMETRIES_ADAPTED_BOOST_ARRAY_UNITS_HPP
#define BOOST_GEOMETRY_GEOMETRIES_ADAPTED_BOOST_ARRAY_UNITS_HPP


#ifdef BOOST_GEOMETRY_ADAPTED_BOOST_ARRAY_UNITS_TAG_DEFINED
#error Include either "boost_array_units_as_point" or \
    "boost_array_units_as_linestring" or "boost_array_units_as_ring" \
    or "boost_array_units_as_multi_point" to adapt a boost_array
#endif

#define BOOST_GEOMETRY_ADAPTED_BOOST_ARRAY_UNITS_TAG_DEFINED


#include <cstddef>

#include <boost/type_traits/is_arithmetic.hpp>

#include <boost/geometry/core/access.hpp>
#include <boost/geometry/core/cs.hpp>
#include <boost/geometry/core/coordinate_dimension.hpp>
#include <boost/geometry/core/coordinate_type.hpp>
#include <boost/geometry/core/tags.hpp>

#include <boost/array.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/unit.hpp>

#include <boost/geometry/units.hpp>
namespace boost { namespace geometry
{


#ifndef DOXYGEN_NO_TRAITS_SPECIALIZATIONS
namespace traits
{


#ifndef DOXYGEN_NO_DETAIL
namespace detail
{


// Create class and specialization to indicate the tag
// for normal cases and the case that the type of the c-array is arithmetic
template <bool>
struct boost_array_units_tag
{
    typedef geometry_not_recognized_tag type;
};


template <>
struct boost_array_units_tag<true>
{
    typedef point_tag type;
};


} // namespace detail
#endif // DOXYGEN_NO_DETAIL


// Assign the point-tag, preventing arrays of points getting a point-tag
template <typename CoordinateUnitType, std::size_t DimensionCount>
struct tag<boost::array< boost::units::quantity<CoordinateUnitType>, DimensionCount> >
    : detail::boost_array_units_tag<boost::is_arithmetic< boost::units::quantity<CoordinateUnitType> >::value> {};


template <typename CoordinateUnitType, std::size_t DimensionCount>
struct coordinate_type<boost::array<boost::units::quantity<CoordinateUnitType>, DimensionCount> >
{
    typedef boost::units::quantity<CoordinateUnitType> type;
};


template <typename CoordinateUnitType, std::size_t DimensionCount>
struct dimension<boost::array<boost::units::quantity<CoordinateUnitType>, DimensionCount> >: boost::mpl::int_<DimensionCount> {};


template <typename CoordinateUnitType, std::size_t DimensionCount, std::size_t Dimension>
struct access<boost::array<boost::units::quantity<CoordinateUnitType>, DimensionCount>, Dimension>
{
    static inline boost::units::quantity<CoordinateUnitType> get(boost::array<boost::units::quantity<CoordinateUnitType>, DimensionCount> const& a)
    {
        return a[Dimension];
    }

    static inline void set(boost::array<boost::units::quantity<CoordinateUnitType>, DimensionCount>& a,
        boost::units::quantity<CoordinateUnitType> const& value)
    {
        a[Dimension] = value;
    }
};

// The library user has
// 1) either to specify the coordinate system
// 2) or include <boost/geometry/geometries/adapted/boost_array_units_@.hpp> where @=cartesian,geographic,...

template<class Unit>
struct point_type<
	boost::array<boost::units::quantity<Unit>, 3ul> 
>{
	typedef model::point<double, 3, cs::units_cartesian<Unit> > type;
};

} // namespace traits
#endif // DOXYGEN_NO_TRAITS_SPECIALIZATIONS



//#ifndef DOXYGEN_NO_DISPATCH
namespace core_dispatch{
template<class  Unit>
struct geometry_id<
	boost::geometry::units_cartesian_tag<Unit> 
> : boost::mpl::int_<1>{
};
}
//#endif
}} // namespace boost::geometry

#endif // BOOST_GEOMETRY_GEOMETRIES_ADAPTED_BOOST_ARRAY_UNITS_HPP

