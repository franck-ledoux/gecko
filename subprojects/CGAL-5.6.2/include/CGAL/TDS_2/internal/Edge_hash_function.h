// Copyright (c) 2005  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/TDS_2/include/CGAL/TDS_2/internal/Edge_hash_function.h $
// $Id: Edge_hash_function.h 98e471849bd 2021-08-26T11:33:39+02:00 Sébastien Loriot
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_INTERNAL_TDS_2_EDGE_HASH_FUNCTION_H
#define CGAL_INTERNAL_TDS_2_EDGE_HASH_FUNCTION_H

#include <CGAL/license/TDS_2.h>

#include <CGAL/basic.h>
#include <CGAL/Handle_hash_function.h>

namespace CGAL {


class Edge_hash_function
  : public Handle_hash_function
{
private:
  typedef Handle_hash_function     Base;

public:
  typedef Base::result_type        result_type;

  template<class Edge>
  result_type operator()(const Edge& e) const
  {
    return (Base::operator()(e.first)) << e.second;
  }
};


} //namespace CGAL


#endif // CGAL_INTERNAL_TDS_2_EDGE_HASH_FUNCTION_H
