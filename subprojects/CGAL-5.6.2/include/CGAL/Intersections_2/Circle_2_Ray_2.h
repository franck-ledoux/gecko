// Copyright (c) 2019
// GeometryFactory.  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Intersections_2/include/CGAL/Intersections_2/Circle_2_Ray_2.h $
// $Id: Circle_2_Ray_2.h 8ba0b41f510 2022-11-22T12:35:10+01:00 Mael Rouxel-Labbé
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Maxime Gimeno

#ifndef CGAL_INTERSECTIONS_2_CIRCLE_2_RAY_2_H
#define CGAL_INTERSECTIONS_2_CIRCLE_2_RAY_2_H

#include <CGAL/Distance_2/Point_2_Ray_2.h>
#include <CGAL/Intersection_traits_2.h>

#include <CGAL/Circle_2.h>
#include <CGAL/Ray_2.h>

namespace CGAL {
namespace Intersections {
namespace internal {

template <class K>
typename K::Boolean
do_intersect(const typename K::Circle_2& c,
             const typename K::Ray_2& r,
             const K&)
{
  return squared_distance(c.center(), r) <= c.squared_radius();
}

template <class K>
typename K::Boolean
do_intersect(const typename K::Ray_2& r,
             const typename K::Circle_2& c,
             const K&)
{
  return squared_distance(c.center(), r) <= c.squared_radius();
}

} // namespace internal
} // namespace Intersections

CGAL_DO_INTERSECT_FUNCTION(Circle_2, Ray_2, 2)

} // namespace CGAL

#endif // CGAL_INTERSECTIONS_2_CIRCLE_2_RAY_2_H
