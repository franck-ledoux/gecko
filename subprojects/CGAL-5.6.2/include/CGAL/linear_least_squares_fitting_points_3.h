// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Principal_component_analysis/include/CGAL/linear_least_squares_fitting_points_3.h $
// $Id: linear_least_squares_fitting_points_3.h 07793738355 2020-03-26T13:31:46+01:00 Sébastien Loriot
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s) : Pierre Alliez and Sylvain Pion and Ankit Gupta

#ifndef CGAL_LINEAR_LEAST_SQUARES_FITTING_POINTS_3_H
#define CGAL_LINEAR_LEAST_SQUARES_FITTING_POINTS_3_H

#include <CGAL/license/Principal_component_analysis.h>


#include <CGAL/basic.h>
#include <CGAL/centroid.h>
#include <CGAL/PCA_util.h>

#include <iterator>

namespace CGAL {

namespace internal {
// fits a plane to a 3D point set
// returns a fitting quality (1 - lambda_min/lambda_max):
//  1 is best (zero variance orthogonally to the fitting line)
//  0 is worst (isotropic case, returns a plane with default direction)
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Plane_3& plane, // best fit plane
                               typename K::Point_3& c,     // centroid
                               const typename K::Point_3*, // used for indirection
                               const K& k,                 // kernel
                               const CGAL::Dimension_tag<0>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Point_3  Point;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K(),tag);

  // assemble covariance matrix
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0., 0., 0., 0. }};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Point*) nullptr,tag, diagonalize_traits);

  // compute fitting plane
  return fitting_plane_3(covariance,c,plane,k,diagonalize_traits);
} // end fit plane to point set

// fits a line to a 3D point set
// returns a fitting quality (1 - lambda_min/lambda_max):
//  1 is best (zero variance orthogonally to the fitting line)
//  0 is worst (isotropic case, returns a line along x axis)
template < typename InputIterator,
           typename K,
           typename DiagonalizeTraits >
typename K::FT
linear_least_squares_fitting_3(InputIterator first,
                               InputIterator beyond,
                               typename K::Line_3& line,  // best fit line
                               typename K::Point_3& c,    // centroid
                               const typename K::Point_3*, // used for indirection
                               const K& k,                // kernel
                               const CGAL::Dimension_tag<0>& tag,
                               const DiagonalizeTraits& diagonalize_traits)
{
  typedef typename K::Point_3  Point;

  // precondition: at least one element in the container.
  CGAL_precondition(first != beyond);

  // compute centroid
  c = centroid(first,beyond,K(),tag);

  // assemble covariance matrix
  typename DiagonalizeTraits::Covariance_matrix covariance = {{ 0., 0., 0., 0., 0., 0. }};
  assemble_covariance_matrix_3(first,beyond,covariance,c,k,(Point*) nullptr,tag, diagonalize_traits);

  // compute fitting line
  return fitting_line_3(covariance,c,line,k,diagonalize_traits);
} // end fit line to point set

} // end namespace internal

} //namespace CGAL

#endif // CGAL_LINEAR_LEAST_SQUARES_FITTING_POINTS_3_H
