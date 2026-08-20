#pragma once

#include <array>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gecko/geom/FacetedGeometry.h>

namespace gecko::python {

    /**
     * @class GeomModelFacade
     * @brief Python-facing façade over gecko::FacetedGeometry: every method takes/returns only
     * primitive types (int/double/string and collections thereof), never a gecko::Point3d, entity
     * pointer, or other internal C++ type — see docs/user-guide/python.md.
     *
     * Entities are identified to Python solely by their Gmsh-style entity_tag() (an int); groups
     * solely by their GroupId's underlying int value.
     */
    class GeomModelFacade {
    public:
        /**
         * @brief Loads the geometric model from a Gmsh MSH file.
         * @param path Path to the .msh file (triangles-only for a boundary rep, or tetrahedra to
         *        also get volumes) — see gecko::FacetedGeometry.
         * @throw std::runtime_error if the file cannot be read.
         */
        explicit GeomModelFacade(const std::string &path);

        /** @brief Gets the number of vertices in the model. @return The vertex count. */
        [[nodiscard]] std::size_t nb_vertices() const noexcept;
        /** @brief Gets the number of curves in the model. @return The curve count. */
        [[nodiscard]] std::size_t nb_curves() const noexcept;
        /** @brief Gets the number of surfaces in the model. @return The surface count. */
        [[nodiscard]] std::size_t nb_surfaces() const noexcept;
        /** @brief Gets the number of volumes in the model. @return The volume count. */
        [[nodiscard]] std::size_t nb_volumes() const noexcept;

        /** @brief Gets the entity_tag() of every vertex in the model. @return The vertex tags. */
        [[nodiscard]] std::vector<int> vertex_tags() const;
        /** @brief Gets the entity_tag() of every curve in the model. @return The curve tags. */
        [[nodiscard]] std::vector<int> curve_tags() const;
        /** @brief Gets the entity_tag() of every surface in the model. @return The surface tags. */
        [[nodiscard]] std::vector<int> surface_tags() const;
        /** @brief Gets the entity_tag() of every volume in the model. @return The volume tags. */
        [[nodiscard]] std::vector<int> volume_tags() const;

        /** @brief Gets the id of every physical group of the model. @return The group ids. */
        [[nodiscard]] std::vector<int> group_ids() const;
        /**
         * @brief Gets the name of a physical group.
         * @param id Group id, from group_ids().
         * @return The group's name.
         * @throw std::out_of_range if @p id is not a known group id.
         */
        [[nodiscard]] std::string group_name(int id) const;
        /**
         * @brief Gets the dimension a physical group applies to.
         * @param id Group id, from group_ids().
         * @return The group's dimension, in [0,3].
         * @throw std::out_of_range if @p id is not a known group id.
         */
        [[nodiscard]] int group_dim(int id) const;
        /**
         * @brief Gets the entities belonging to a physical group.
         * @param id Group id, from group_ids().
         * @return One (dimension, entity_tag) pair per entity of the group; empty if @p id is not
         *         a known group id.
         */
        [[nodiscard]] std::vector<std::pair<int, int>> group_entities(int id) const;

        /**
         * @brief Gets the positions of every node of the model's backing (faceted) mesh.
         * @return One (x,y,z) triple per node, indexable by the indices mesh_triangles()/
         *         mesh_tets() return.
         */
        [[nodiscard]] std::vector<std::array<double, 3>> mesh_vertices() const;
        /**
         * @brief Gets the triangles of the model's backing mesh.
         * @return One triple of mesh_vertices() indices per triangle; empty if the model has none.
         */
        [[nodiscard]] std::vector<std::array<int, 3>> mesh_triangles() const;
        /**
         * @brief Gets the tetrahedra of the model's backing mesh.
         * @return One quadruple of mesh_vertices() indices per tetrahedron; empty if the model was
         *         read from a boundary-representation (triangles-only) file.
         */
        [[nodiscard]] std::vector<std::array<int, 4>> mesh_tets() const;

        /**
         * @brief Gives access to the wrapped model. Not exposed to Python — used only by
         * BlockingFacade to build the real gecko::Blocking against this model.
         * @return The wrapped model.
         */
        [[nodiscard]] const FacetedGeometry &native() const { return *m_model; }

    private:
        std::unique_ptr<FacetedGeometry> m_model;
    };

} // namespace gecko::python
