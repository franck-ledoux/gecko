// Copyright (c) 2016  GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/Installation/include/CGAL/Installation/internal/enable_third_party_libraries.h $
// $Id: enable_third_party_libraries.h 6b71a94254c 2023-04-07T13:59:38+02:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_INTERNAL_ENABLE_THIRD_PARTY_LIBRARIES_H
#define CGAL_INTERNAL_ENABLE_THIRD_PARTY_LIBRARIES_H

// GMP and MPFR are highly recommended in CGAL.
#define CGAL_USE_GMP 1
#define CGAL_USE_MPFR 1

#if CGAL_DISABLE_GMP && ! defined(CGAL_NO_GMP)
#  define CGAL_NO_GMP 1
#endif

#if CGAL_NO_GMP || CGAL_NO_MPFR
#  undef CGAL_USE_MPFR
#  undef CGAL_USE_GMP
#endif

#if defined(__has_include) && ( ! defined _MSC_VER || _MSC_VER > 1900)
#  if CGAL_USE_GMP && ! __has_include(<gmp.h>)
#    undef CGAL_USE_GMP
#    undef CGAL_USE_MPFR
#  elif CGAL_USE_MPFR && ! __has_include(<mpfr.h>)
#    undef CGAL_USE_GMP
#    undef CGAL_USE_MPFR
#  endif // CGAL_USE_MPFR and no <mpfr.h>
#endif // __has_include

#if CGAL_USE_GMP && CGAL_USE_MPFR && ! CGAL_NO_CORE
#  define CGAL_USE_CORE 1
#endif

#endif // CGAL_INTERNAL_ENABLE_THIRD_PARTY_LIBRARIES_H
