#pragma once

#include <variant>
#include <gecko/geom/FacetedEntities.h>

namespace gecko {

    /**
     * @brief Variant holding a const pointer to one of the 4 faceted CAD entity types
     * (vertex, curve, surface or volume), used to reference any faceted geometric entity uniformly.
     */
    using FacetedEntityRefVariant =
        std::variant<const FacetedVertex *, const FacetedCurve *, const FacetedSurface *, const FacetedVolume *>;

} // namespace gecko
