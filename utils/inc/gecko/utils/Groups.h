#pragma once
#include <string>
#include <vector>
#include <gecko/utils/Types.h>

namespace gecko {

    // Identifiant fort pour éviter toute confusion d'index
    struct GroupTag {};
    using GroupId = StrongId<GroupTag, std::uint32_t>;

    /** @brief Dimension an entity group applies to. */
    enum class GroupDim : std::uint8_t { Dim0 = 0, Dim1 = 1, Dim2 = 2, Dim3 = 3, MultiDim = 4, Undefined = 255 };

    /** @brief Full definition of a named group of entities. */
    struct Group {
        GroupId id;                              ///< Unique id of the group.
        GroupDim dimension{GroupDim::Undefined}; ///< Dimension the group applies to.
        std::string name;                        ///< Name of the group (e.g. "Farfield", "Inlet_Velocity", "Fluid").
    };

    /**
     * @class GroupRegistry
     * @brief Stores and deduplicates the set of named groups defined over the model entities.
     */
    class GroupRegistry {
    public:
        /**
         * @brief Creates a group or retrieves the id of an already existing one.
         * @param name Name of the group.
         * @param dim Dimension the group applies to.
         * @return The id of the (possibly newly created) group matching @p name and @p dim.
         */
        GroupId add_group(std::string_view name, GroupDim dim) {
            // Vérifier s'il existe déjà
            for (const auto &grp : m_groups) {
                if (grp.name == name && grp.dimension == dim) {
                    return grp.id;
                }
            }

            // Sinon le créer
            auto new_id = GroupId{static_cast<std::uint32_t>(m_groups.size())};
            m_groups.push_back(Group{new_id, dim, std::string(name)});
            return new_id;
        }

        /**
         * @brief Accesses a group by its id.
         * @param id Id of the group to retrieve.
         * @return Const reference to the requested group.
         * @pre @p id must be a valid id previously returned by add_group().
         */
        [[nodiscard]] const Group &get_group(GroupId id) const { return m_groups[static_cast<std::uint32_t>(id)]; }

    private:
        std::vector<Group> m_groups;
    };
} // namespace gecko