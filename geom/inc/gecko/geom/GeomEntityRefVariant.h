#pragma once

#include <variant>
#include "gecko/geom/GeomEntities.h"

namespace gecko {

    /**
     * @brief Variant holding a const pointer to one of the 4 CAD entity types
     * (vertex, curve, surface or volume), used to reference any geometric entity uniformly.
     */
    using GeomEntityRefVariant =
        std::variant<const GeomVertex *, const GeomCurve *, const GeomSurface *, const GeomVolume *>;

} // namespace gecko