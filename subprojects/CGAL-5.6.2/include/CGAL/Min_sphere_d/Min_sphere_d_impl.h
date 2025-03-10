// Copyright (c) 1997  ETH Zurich (Switzerland).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Bounding_volumes/include/CGAL/Min_sphere_d/Min_sphere_d_impl.h $
// $Id: Min_sphere_d_impl.h 74e4d89cbc1 2022-09-27T10:42:05+01:00 Andreas Fabri
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Schoenherr <sven@inf.fu-berlin.de>
//                 Bernd Gaertner

#include <iterator>

namespace CGAL {

// Class implementation (continued)
// ================================
// I/O
// ---
template < class Traits >
std::ostream&
operator << ( std::ostream& os, const Min_sphere_d<Traits>& min_sphere)
{
    typedef typename Min_sphere_d<Traits>::Point  Point;

    switch ( IO::get_mode( os)) {

      case IO::PRETTY:
        os << std::endl;
        os << "Min_sphere_d( |P| = " << min_sphere.number_of_points()
           << ", |S| = " << min_sphere.number_of_support_points()
           << std::endl;
        os << "  P = {" << std::endl;
        os << "    ";
        std::copy( min_sphere.points_begin(), min_sphere.points_end(),
              std::ostream_iterator<Point>( os, ",\n    "));
        os << "}" << std::endl;
        os << "  S = {" << std::endl;
        os << "    ";
        std::copy( min_sphere.support_points_begin(),
              min_sphere.support_points_end(),
              std::ostream_iterator<Point>( os, ",\n    "));
        os << "}" << std::endl;
        os << "  center = " << min_sphere.center() << std::endl;
        os << "  squared radius = " << min_sphere.squared_radius()
           << std::endl;
        os << ")" << std::endl;
        break;

      case IO::ASCII:
        os << min_sphere.number_of_points() << std::endl;
        std::copy( min_sphere.points_begin(), min_sphere.points_end(),
              std::ostream_iterator<Point>( os, "\n"));
        break;

      case IO::BINARY:
        os << min_sphere.number_of_points() << " ";
        std::copy( min_sphere.points_begin(), min_sphere.points_end(),
              std::ostream_iterator<Point>( os, " "));
        break;

      default:
        CGAL_assertion_msg
            ( false, "IO::get_mode( os) invalid!");
        break; }

    return( os);
}

template < class Traits >
std::istream&
operator >> ( std::istream& is, Min_sphere_d<Traits>& min_sphere)
{
    switch ( IO::get_mode( is)) {

      case IO::PRETTY:
        std::cerr << std::endl;
        std::cerr << "Stream must be in ASCII or binary mode" << std::endl;
        break;

      case IO::ASCII:
      case IO::BINARY:
      {
        min_sphere.clear();
        int n; is >> n;
        typename Min_sphere_d<Traits>::Point p;
        for (int i=0; i<n; ++i) {
            is >> p;
            min_sphere.insert (p);
        }
      } break;

      default:
        CGAL_assertion_msg( false, "IO::mode invalid!");
        break;
 }

    return( is);
}



} //namespace CGAL

// ===== EOF ==================================================================
