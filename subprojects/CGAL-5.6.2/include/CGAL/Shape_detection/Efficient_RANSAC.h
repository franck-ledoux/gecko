// Copyright (c) 2015 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Shape_detection/include/CGAL/Shape_detection/Efficient_RANSAC.h $
// $Id: Efficient_RANSAC.h 1c7ba098636 2021-06-30T12:46:41+02:00 Dmitry Anisimov
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sven Oesau, Yannick Verdie, Clément Jamin, Pierre Alliez
//

#ifndef CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_HEADERS_H
#define CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_HEADERS_H

/// \cond SKIP_IN_MANUAL
#include <CGAL/license/Shape_detection.h>
/// \endcond

/**
* \ingroup PkgShapeDetectionRef
* \file CGAL/Shape_detection/Efficient_RANSAC.h
* A convenience header that includes all classes related to the efficient RANSAC algorithm.
*/

#include <CGAL/Shape_detection/Efficient_RANSAC/Efficient_RANSAC.h>
#include <CGAL/Shape_detection/Efficient_RANSAC/Efficient_RANSAC_traits.h>

#include <CGAL/Shape_detection/Efficient_RANSAC/Cone.h>
#include <CGAL/Shape_detection/Efficient_RANSAC/Plane.h>
#include <CGAL/Shape_detection/Efficient_RANSAC/Torus.h>
#include <CGAL/Shape_detection/Efficient_RANSAC/Sphere.h>
#include <CGAL/Shape_detection/Efficient_RANSAC/Cylinder.h>
#include <CGAL/Shape_detection/Efficient_RANSAC/Shape_base.h>

#include <CGAL/Shape_detection/Efficient_RANSAC/Octree.h>
#include <CGAL/Shape_detection/Efficient_RANSAC/property_map.h>

#endif // CGAL_SHAPE_DETECTION_EFFICIENT_RANSAC_HEADERS_H
