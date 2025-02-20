// Copyright (c) 2008 GeometryFactory (France). All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL: https://github.com/CGAL/cgal/blob/v5.6.2/GraphicsView/include/CGAL/auto_link/Qt.h $
// $Id: Qt.h 07793738355 2020-03-26T13:31:46+01:00 Sébastien Loriot
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Laurent Rineau

#ifndef CGAL_AUTO_LINK_QT_H
#define CGAL_AUTO_LINK_QT_H

#include <CGAL/config.h>
#include <QtCore/qglobal.h>

#if (! defined (CGAL_NO_AUTOLINK_QT))
#if ( ! defined( CGAL_EXPORTS )  && (! defined ( CGAL_Qt5_EXPORTS )))

// If CGAL_EXPORTS is defined it means that we are building the CGAL
// library as a DLL. The CGAL.dll does not really depend on CGAL_Qt,
// whatever the header inclusion graph says.

#define CGAL_LIB_NAME CGAL_Qt5

#include <CGAL/auto_link/auto_link.h>

#endif // CGAL_EXPORTS
#endif // CGAL_NO_AUTOLINK_QT

#endif // CGAL_AUTO_LINK_QT_H
