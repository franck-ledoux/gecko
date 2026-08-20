#include "GeomModelFacade.h"

#include <stdexcept>

#include <gecko/utils/Groups.h>

namespace gecko::python {

    GeomModelFacade::GeomModelFacade(const std::string &path) : m_model(std::make_unique<FacetedGeometry>(path)) {}

    std::size_t GeomModelFacade::nb_vertices() const noexcept { return m_model->nb_vertices(); }
    std::size_t GeomModelFacade::nb_curves() const noexcept { return m_model->nb_curves(); }
    std::size_t GeomModelFacade::nb_surfaces() const noexcept { return m_model->nb_surfaces(); }
    std::size_t GeomModelFacade::nb_volumes() const noexcept { return m_model->nb_volumes(); }

    std::vector<int> GeomModelFacade::vertex_tags() const {
        std::vector<int> tags;
        tags.reserve(m_model->nb_vertices());
        for (const auto &v : m_model->vertices()) tags.push_back(v.entity_tag());
        return tags;
    }

    std::vector<int> GeomModelFacade::curve_tags() const {
        std::vector<int> tags;
        tags.reserve(m_model->nb_curves());
        for (const auto &c : m_model->curves()) tags.push_back(c.entity_tag());
        return tags;
    }

    std::vector<int> GeomModelFacade::surface_tags() const {
        std::vector<int> tags;
        tags.reserve(m_model->nb_surfaces());
        for (const auto &s : m_model->surfaces()) tags.push_back(s.entity_tag());
        return tags;
    }

    std::vector<int> GeomModelFacade::volume_tags() const {
        std::vector<int> tags;
        tags.reserve(m_model->nb_volumes());
        for (const auto &v : m_model->volumes()) tags.push_back(v.entity_tag());
        return tags;
    }

    std::vector<int> GeomModelFacade::group_ids() const {
        std::vector<int> ids;
        ids.reserve(m_model->groups().size());
        for (const auto &g : m_model->groups()) ids.push_back(static_cast<int>(g.id.value));
        return ids;
    }

    namespace {
        const Group &find_group(const FacetedGeometry &model, int id) {
            for (const auto &g : model.groups()) {
                if (static_cast<int>(g.id.value) == id) return g;
            }
            throw std::out_of_range("GeomModel: unknown group id " + std::to_string(id));
        }
    } // namespace

    std::string GeomModelFacade::group_name(int id) const { return find_group(*m_model, id).name; }

    int GeomModelFacade::group_dim(int id) const { return static_cast<int>(find_group(*m_model, id).dimension); }

    std::vector<std::pair<int, int>> GeomModelFacade::group_entities(int id) const {
        std::vector<std::pair<int, int>> result;
        if (id < 0) return result;
        for (const auto &entity : m_model->entities(GroupId{static_cast<std::uint32_t>(id)})) {
            std::visit([&result](const auto *e) { result.emplace_back(static_cast<int>(e->dimension()), e->entity_tag()); },
                       entity);
        }
        return result;
    }

} // namespace gecko::python
