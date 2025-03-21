// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Nef_S2/include/CGAL/Nef_S2/Sphere_direction.h $
// $Id: Sphere_direction.h 78e2d5e4d2b 2023-03-02T13:42:35+01:00 Laurent Rineau
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel  <seel@mpi-sb.mpg.de>

#ifndef CGAL_SPHERE_DIRECTION_H
#define CGAL_SPHERE_DIRECTION_H

#include <CGAL/license/Nef_S2.h>


#include <CGAL/basic.h>
#include <CGAL/Kernel/global_functions_3.h>

namespace CGAL {

/*{\Manpage{Sphere_direction}{R}{Directions on the unit sphere}{c}}*/

template <class R_> class Sphere_direction : public R_::Plane_3 {

/*{\Mdefinition An object |\Mvar| of type |\Mname| is a direction
on the surface of the unit sphere.  Such directions can be used to
describe walks that are part of great circles.}*/

public:

/*{\Mtypes 5}*/
typedef R_ R;
/*{\Mtypemember representation class.}*/
typedef typename R::RT RT;
/*{\Mtypemember ring type.}*/

typedef typename R_::Point_3 Point_3;
typedef typename R_::Plane_3 Plane_3;
typedef typename R_::Plane_3 Base;
typedef Sphere_direction<R_> Self;

/*{\Mcreation 5}*/
Sphere_direction() : Base() {}
/*{\Mcreate creates some direction.}*/

Sphere_direction(const Sphere_circle<R>& c)
/*{\Mcreate creates the direction corresponding to the circle |c|.}*/
  : Base(c) {}

Sphere_direction(const Sphere_point<R>& p, const Sphere_point<R>&q)
  : Base(CGAL::ORIGIN,p,q)
/*{\Mcreate creates a direction that describes the orientation of
the great circle through $p$ and $q$ (oriented such that the segment
$pq$ is the shorter one of the two possible ones. \precond $p$ and $q$
are not opposite on $S_2$.}*/
{ CGAL_assertion(p!=q.opposite());
  Point_3 p4 = CGAL::ORIGIN + ((Base*) this)->orthogonal_vector();
  if ( R().orientation_3_object()(CGAL::ORIGIN,p,q,p4) != CGAL::POSITIVE )
    *this = Sphere_direction(opposite());
}

Sphere_direction(const typename R::Plane_3& h)
/*{\Xcreate creates the direction corresponding to the plane |h|.
\precond |h| contains the origin.}*/
 : Base(h) { CGAL_assertion(h.d() == 0); }

/*{\Moperations 4 2}*/

Sphere_direction<R> opposite() const
/*{\Mop returns the opposite of |\Mvar|.}*/
{ return Base::opposite(); }

Plane_3 plane() const { return Base(*this); }
/*{\Xop returns the plane supporting |\Mvar|.}*/

}; // Sphere_direction<R>


/* We have:
   1) all directions fixed at p
   2) d1==d3 possible
   return true iff d1,d2,d3 are strictly ccw ordered around p
   Note: Sphere_directions are Plane_3
         we therefore compare the normal vectors of the planes
         that underly the directions d1,d2,d3 in the plane
         through 0 and orthogonal to the vector p-0
 */

template <typename R>
bool strictly_ordered_ccw_at(const Sphere_point<R>& p,
  const Sphere_direction<R>& d1,
  const Sphere_direction<R>& d2,
  const Sphere_direction<R>& d3)
{ CGAL_assertion(d1.has_on(p) && d2.has_on(p) && d3.has_on(p));
  typedef typename R::Vector_3 Vector_3;
  Vector_3 v0 = p - CGAL::ORIGIN;
  Vector_3 v1 = d1.orthogonal_vector();
  Vector_3 v2 = d2.orthogonal_vector();
  Vector_3 v3 = d3.orthogonal_vector();

  if ( d1 == d3 ) return false;
  typename R::Orientation_3 orientation = R().orientation_3_object();
  if ( orientation(v0,v1,v3) == CGAL::POSITIVE ) {
    return orientation(v0,v1,v2) == CGAL::POSITIVE &&
           orientation(v0,v3,v2) == CGAL::NEGATIVE;
  } else {
    return orientation(v0,v1,v2) == CGAL::POSITIVE ||
           orientation(v0,v3,v2) == CGAL::NEGATIVE;
  }
}


} //namespace CGAL
#endif //CGAL_SPHERE_DIRECTION_H
