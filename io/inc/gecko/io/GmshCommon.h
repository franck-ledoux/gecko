#pragma once
#include <string_view>

#include <gecko/mesh/UnstructuredMesh.h>

namespace gecko::io {

    /**
     * @brief Name of the per-edge/per-face/per-cell `Variable<GroupId>` used by
     * GmshMeshReader/GmshMeshWriter to store each mesh element's Gmsh physical group (element tag 0).
     */
    inline constexpr std::string_view PHYSICAL_GROUP_VARIABLE = "gmsh_physical_group";

    /**
     * @brief Name of the per-edge/per-face/per-cell `Variable<Int>` used by
     * GmshMeshReader/GmshMeshWriter to store each mesh element's Gmsh elementary entity tag
     * (element tag 1: the id of the geometric model entity the element is classified on). A value
     * of 0 means the element had no elementary entity tag in the file.
     */
    inline constexpr std::string_view ENTITY_TAG_VARIABLE = "gmsh_entity_tag";

    /** @brief Gmsh MSH element type code for a 2-node line, used for mesh edges regardless of topology. */
    inline constexpr int EDGE_TYPE = 1;

    /** @brief Gmsh MSH element type code for a 1-node point, used for mesh nodes regardless of topology. */
    inline constexpr int POINT_TYPE = 15;

    /**
     * @brief Maps a mesh topology traits type to the Gmsh MSH element type codes used for its
     * Face and Cell elements. Only defined for gecko::SimplicialTraits and gecko::CubicTraits.
     * @tparam TopologyTraits SimplicialTraits or CubicTraits.
     */
    template<typename TopologyTraits>
    struct GmshElementTraits;

    /**
     * @brief Gmsh element codes for a simplicial mesh: 3-node triangle faces, 4-node tetrahedron cells.
     */
    template<>
    struct GmshElementTraits<SimplicialTraits> {
        static constexpr int face_type = 2; ///< Gmsh element type 2: 3-node triangle.
        static constexpr int cell_type = 4; ///< Gmsh element type 4: 4-node tetrahedron.
    };

    /**
     * @brief Gmsh element codes for a cubic mesh: 4-node quadrangle faces, 8-node hexahedron cells.
     */
    template<>
    struct GmshElementTraits<CubicTraits> {
        static constexpr int face_type = 3; ///< Gmsh element type 3: 4-node quadrangle.
        static constexpr int cell_type = 5; ///< Gmsh element type 5: 8-node hexahedron.
    };

} // namespace gecko::io
