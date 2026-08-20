#include "BlockingFacade.h"

#include <stdexcept>

#include <gecko/io/VtkMeshWriter.h>

namespace gecko::python {

    namespace {
        template<std::size_t K>
        std::array<Point3d, K> to_points(const std::vector<std::array<double, 3>> &corners, const char *who) {
            if (corners.size() != K) {
                throw std::invalid_argument(std::string(who) + ": expected " + std::to_string(K) + " corners, got " +
                                            std::to_string(corners.size()));
            }
            std::array<Point3d, K> points{};
            for (std::size_t i = 0; i < K; ++i) {
                points[i] = Point3d(corners[i][0], corners[i][1], corners[i][2]);
            }
            return points;
        }

        using ImplVariant = std::variant<BlockingFacade::Impl<1>, BlockingFacade::Impl<3>>;

        ImplVariant make_impl(const FacetedGeometry &geom, int degree) {
            if (degree == 1) return ImplVariant{std::in_place_type<BlockingFacade::Impl<1>>, geom};
            if (degree == 3) return ImplVariant{std::in_place_type<BlockingFacade::Impl<3>>, geom};
            throw std::invalid_argument("Blocking: degree must be 1 (straight) or 3 (cubic Bezier), got " +
                                        std::to_string(degree));
        }
    } // namespace

    BlockingFacade::BlockingFacade(const GeomModelFacade &model, int degree)
        : m_impl(make_impl(model.native(), degree)) {}

    int BlockingFacade::create_quad_block(const std::vector<std::array<double, 3>> &corners) {
        const auto points = to_points<4>(corners, "Blocking.create_quad_block");
        return std::visit(
            [&points](auto &impl) {
                const int id = impl.next_face_id++;
                impl.faces_by_id.emplace(id, impl.blocking.create_quad_block(points));
                return id;
            },
            m_impl);
    }

    int BlockingFacade::create_hex_block(const std::vector<std::array<double, 3>> &corners) {
        const auto points = to_points<8>(corners, "Blocking.create_hex_block");
        return std::visit(
            [&points](auto &impl) {
                impl.blocking.create_hex_block(points);
                return impl.next_block_id++;
            },
            m_impl);
    }

    void BlockingFacade::build_connectivity() {
        std::visit([](auto &impl) { impl.blocking.build_connectivity(); }, m_impl);
    }

    void BlockingFacade::classify(double tol_vertex, double tol_curve_surface) {
        std::visit(
            [tol_vertex, tol_curve_surface](auto &impl) { impl.blocking.classify(tol_vertex, tol_curve_surface); },
            m_impl);
    }

    std::size_t BlockingFacade::nb_cells(int dim) const {
        return std::visit(
            [dim](const auto &impl) -> std::size_t {
                switch (dim) {
                    case 0:
                        return impl.blocking.template nb_cells<0>();
                    case 1:
                        return impl.blocking.template nb_cells<1>();
                    case 2:
                        return impl.blocking.template nb_cells<2>();
                    case 3:
                        return impl.blocking.template nb_cells<3>();
                    default:
                        throw std::invalid_argument("Blocking.nb_cells: dim must be in [0,3], got " +
                                                    std::to_string(dim));
                }
            },
            m_impl);
    }

    bool BlockingFacade::is_valid_topology() const {
        return std::visit([](const auto &impl) { return impl.blocking.is_valid_topology(); }, m_impl);
    }

    bool BlockingFacade::is_purely_2d() const {
        return std::visit([](const auto &impl) { return impl.blocking.is_purely_2d(); }, m_impl);
    }

    namespace {
        template<typename TImpl>
        auto find_face(TImpl &impl, int face_id) {
            const auto it = impl.faces_by_id.find(face_id);
            if (it == impl.faces_by_id.end()) {
                throw std::out_of_range("Blocking: unknown face id " + std::to_string(face_id));
            }
            return it;
        }
    } // namespace

    bool BlockingFacade::can_delete_face(int face_id) const {
        return std::visit(
            [face_id](const auto &impl) { return impl.blocking.can_delete_face(find_face(impl, face_id)->second); },
            m_impl);
    }

    void BlockingFacade::delete_face(int face_id) {
        std::visit(
            [face_id](auto &impl) {
                const auto it = find_face(impl, face_id);
                impl.blocking.delete_face(it->second);
                impl.faces_by_id.erase(it);
            },
            m_impl);
    }

    void BlockingFacade::write_vtk(int subdivisions, const std::string &path) {
        if (subdivisions < 1) {
            throw std::invalid_argument("Blocking.write_vtk: subdivisions must be >= 1, got " +
                                        std::to_string(subdivisions));
        }
        std::visit(
            [subdivisions, &path](auto &impl) {
                const auto mesh = impl.blocking.to_mesh(static_cast<SizeT>(subdivisions));
                io::VtkMeshWriter<CubicTraits>::write(path, mesh);
            },
            m_impl);
    }

} // namespace gecko::python
