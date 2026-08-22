#include "BlockingFacade.h"

#include <stdexcept>
#include <type_traits>

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

        void check_subdivisions(int subdivisions, const char *who) {
            if (subdivisions < 1) {
                throw std::invalid_argument(std::string(who) + ": subdivisions must be >= 1, got " +
                                            std::to_string(subdivisions));
            }
        }

        using ImplVariant = std::
            variant<BlockingFacade::Impl<1>, BlockingFacade::Impl<2>, BlockingFacade::Impl<3>, BlockingFacade::Impl<4>>;

        ImplVariant make_impl(const FacetedGeometry &geom, int degree) {
            switch (degree) {
                case 1:
                    return ImplVariant{std::in_place_type<BlockingFacade::Impl<1>>, geom};
                case 2:
                    return ImplVariant{std::in_place_type<BlockingFacade::Impl<2>>, geom};
                case 3:
                    return ImplVariant{std::in_place_type<BlockingFacade::Impl<3>>, geom};
                case 4:
                    return ImplVariant{std::in_place_type<BlockingFacade::Impl<4>>, geom};
                default:
                    throw std::invalid_argument(
                        "Blocking: degree must be in [" + std::to_string(BlockingFacade::MIN_DEGREE) + "," +
                        std::to_string(BlockingFacade::MAX_DEGREE) + "] (1 = straight), got " + std::to_string(degree));
            }
        }
    } // namespace

    BlockingFacade::BlockingFacade(const GeomModelFacade &model, int degree)
        : m_impl(make_impl(model.native(), degree)) {}

    int BlockingFacade::degree() const {
        return std::visit([](const auto &impl) { return std::decay_t<decltype(impl)>::DEGREE; }, m_impl);
    }

    int BlockingFacade::create_quad_block(const std::vector<std::array<double, 3>> &corners) {
        const auto points = to_points<4>(corners, "Blocking.create_quad_block");
        return std::visit(
            [&points](auto &impl) {
                const int id = impl.next_face_id++;
                impl.faces_by_id.emplace(id, impl.blocking.create_quad_block(points));
                impl.index_new_nodes();
                return id;
            },
            m_impl);
    }

    int BlockingFacade::create_hex_block(const std::vector<std::array<double, 3>> &corners) {
        const auto points = to_points<8>(corners, "Blocking.create_hex_block");
        return std::visit(
            [&points](auto &impl) {
                const auto block = impl.blocking.create_hex_block(points);
                impl.index_new_nodes();
                // Reported as a position in the block traversal rather than a counter of its own:
                // that is the one blocks are addressed by everywhere else here, and a counter no
                // other method accepts is not an id, only a number.
                auto &map = impl.blocking.cmap();
                int index = 0;
                for (auto it = map.template attributes<3>().begin(), itend = map.template attributes<3>().end();
                     it != itend;
                     ++it, ++index) {
                    if (it == block) return index;
                }
                return -1;
            },
            m_impl);
    }

    void BlockingFacade::build_connectivity() {
        std::visit(
            [](auto &impl) {
                impl.blocking.build_connectivity();
                // Sewing merges each pair of coincident corners into one attribute, discarding the
                // other — ids still pointing at the discarded one have to go.
                impl.forget_stale_nodes();
            },
            m_impl);
    }

    void BlockingFacade::classify(double tol_vertex, double tol_curve, double tol_surface) {
        std::visit([tol_vertex, tol_curve, tol_surface](
                       auto &impl) { impl.blocking.classify(tol_vertex, tol_curve, tol_surface); },
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

    namespace {
        /** @brief The block sitting at position @p index of a blocking's own block traversal — the
         * order `mesh_hexes()` emits its cells in, and the one `block_volumes()` reports in.
         * @tparam TImpl The per-degree Impl alternative.
         * @param impl The blocking state.
         * @param index The position to resolve.
         * @return That block.
         * @throw std::out_of_range if the blocking has no such block. */
        template<typename TImpl>
        auto block_at(TImpl &impl, int index) {
            auto &map = impl.blocking.cmap();
            if (index >= 0) {
                int seen = 0;
                for (auto it = map.template attributes<3>().begin(), itend = map.template attributes<3>().end();
                     it != itend;
                     ++it, ++seen) {
                    if (seen == index) return it;
                }
            }
            throw std::out_of_range("Blocking: unknown block index " + std::to_string(index));
        }
    } // namespace

    void BlockingFacade::delete_block(int block_index) {
        std::visit(
            [block_index](auto &impl) {
                impl.blocking.delete_block(block_at(impl, block_index));
                // A deletion takes corner nodes with it, so ids pointing at them have to go too —
                // an id map holding freed attributes is worse than one that has forgotten them.
                impl.forget_stale_nodes();
            },
            m_impl);
    }

    void BlockingFacade::write_vtk(int subdivisions, const std::string &path) {
        check_subdivisions(subdivisions, "Blocking.write_vtk");
        std::visit(
            [subdivisions, &path](auto &impl) {
                using Blocking = std::decay_t<decltype(impl.blocking)>;
                const auto mesh = impl.blocking.to_mesh(static_cast<SizeT>(subdivisions));
                io::VtkMeshWriter<CubicTraits>::write(path,
                                                      mesh,
                                                      {std::string(Blocking::NODE_CLASSIFICATION_DIM_VARIABLE),
                                                       std::string(Blocking::NODE_CLASSIFICATION_TAG_VARIABLE)});
            },
            m_impl);
    }

    namespace {
        /** @brief The edge sitting at position @p index of a blocking's own edge traversal — the
         * order every edge-indexed accessor here shares.
         * @tparam TImpl The per-degree Impl alternative.
         * @param impl The blocking state.
         * @param index The position to resolve.
         * @return That edge.
         * @throw std::out_of_range if the blocking has no such edge. */
        template<typename TImpl>
        auto edge_at(TImpl &impl, int index) {
            auto &map = impl.blocking.cmap();
            if (index >= 0) {
                int seen = 0;
                for (auto it = map.template attributes<1>().begin(), itend = map.template attributes<1>().end();
                     it != itend;
                     ++it, ++seen) {
                    if (seen == index) return it;
                }
            }
            throw std::out_of_range("Blocking: unknown edge index " + std::to_string(index));
        }
    } // namespace

    std::vector<double> BlockingFacade::block_volumes(int subdivisions) {
        check_subdivisions(subdivisions, "Blocking.block_volumes");
        return std::visit([subdivisions](auto &impl) { return impl.blocking.block_volumes(subdivisions); }, m_impl);
    }

    std::vector<int> BlockingFacade::sheet_edges(int edge_index) {
        return std::visit(
            [edge_index](auto &impl) {
                const auto sheet = impl.blocking.find_sheet(edge_at(impl, edge_index));
                std::vector<int> indices;
                if (!sheet.has_value()) return indices;

                // Back from attribute handles to the positions the caller speaks in, in one pass.
                auto &map = impl.blocking.cmap();
                int position = 0;
                for (auto it = map.template attributes<1>().begin(), itend = map.template attributes<1>().end();
                     it != itend;
                     ++it, ++position) {
                    for (const auto &se : sheet->edges) {
                        if (se.edge == it) {
                            indices.push_back(position);
                            break;
                        }
                    }
                }
                return indices;
            },
            m_impl);
    }

    std::vector<std::array<double, 3>> BlockingFacade::sheet_cut_points(int edge_index, double param) {
        return std::visit(
            [edge_index, param](auto &impl) {
                std::vector<std::array<double, 3>> points;
                const auto sheet = impl.blocking.find_sheet(edge_at(impl, edge_index));
                if (!sheet.has_value()) return points;

                // Walked in the blocking's own edge order rather than the sheet's: find_sheet()
                // returns its edges keyed by attribute handle, whose ordering is the container's and
                // need not match the traversal sheet_edges() reports positions in. Re-walking is what
                // makes the 2 line up entry for entry, which is what a caller pairing them expects.
                auto &map = impl.blocking.cmap();
                points.reserve(sheet->edges.size());
                for (auto it = map.template attributes<1>().begin(), itend = map.template attributes<1>().end();
                     it != itend;
                     ++it) {
                    for (const auto &se : sheet->edges) {
                        if (se.edge != it) continue;
                        const auto p = impl.blocking.cut_point(se, param);
                        points.push_back({p.x(), p.y(), p.z()});
                        break;
                    }
                }
                return points;
            },
            m_impl);
    }

    bool BlockingFacade::cut_sheet(int edge_index, double param) {
        return std::visit(
            [edge_index, param](auto &impl) {
                if (!impl.blocking.cut_sheet(edge_at(impl, edge_index), param)) return false;
                // A cut creates corner nodes, which have to become addressable like any other.
                impl.index_new_nodes();
                return true;
            },
            m_impl);
    }
    namespace {
        template<typename TImpl>
        auto find_node(TImpl &impl, int node_id) {
            const auto it = impl.nodes_by_id.find(node_id);
            if (it == impl.nodes_by_id.end()) {
                throw std::out_of_range("Blocking: unknown node id " + std::to_string(node_id));
            }
            return it;
        }
    } // namespace

    std::vector<int> BlockingFacade::node_ids() const {
        return std::visit(
            [](const auto &impl) {
                std::vector<int> ids;
                ids.reserve(impl.nodes_by_id.size());
                for (const auto &[id, node] : impl.nodes_by_id)
                    ids.push_back(id);
                std::sort(ids.begin(), ids.end());
                return ids;
            },
            m_impl);
    }

    std::array<double, 3> BlockingFacade::node_position(int node_id) const {
        return std::visit(
            [node_id](const auto &impl) {
                const auto &p = find_node(impl, node_id)->second->info().point;
                return std::array<double, 3>{p.x(), p.y(), p.z()};
            },
            m_impl);
    }

    void BlockingFacade::move_node(int node_id, double x, double y, double z) {
        std::visit([node_id, x, y, z](
                       auto &impl) { impl.blocking.move_node(find_node(impl, node_id)->second, Point3d(x, y, z)); },
                   m_impl);
    }

    void BlockingFacade::snap_node(int node_id, double tol_vertex, double tol_curve, double tol_surface) {
        std::visit(
            [node_id, tol_vertex, tol_curve, tol_surface](auto &impl) {
                impl.blocking.snap_node(find_node(impl, node_id)->second, tol_vertex, tol_curve, tol_surface);
            },
            m_impl);
    }

    namespace {
        /** @brief The dimension a cell's `geom_targets` puts it on, or -1 when unclassified. Every
         * target of one cell sits at the same dimension (see `Blocking::classify()`), so the first
         * speaks for all. */
        int classification_dim(const std::vector<std::pair<GroupDim, Int>> &ATargets) {
            return ATargets.empty() ? -1 : static_cast<int>(ATargets.front().first);
        }
    } // namespace

    std::vector<int> BlockingFacade::node_classification_dims() const {
        return std::visit(
            [](const auto &impl) {
                std::vector<int> ids;
                ids.reserve(impl.nodes_by_id.size());
                for (const auto &[id, node] : impl.nodes_by_id)
                    ids.push_back(id);
                std::sort(ids.begin(), ids.end()); // node_ids()' order

                std::vector<int> dims;
                dims.reserve(ids.size());
                for (const int id : ids) {
                    dims.push_back(classification_dim(impl.nodes_by_id.at(id)->info().geom_targets));
                }
                return dims;
            },
            m_impl);
    }

    std::vector<int> BlockingFacade::edge_classification_dims() const {
        return std::visit(
            [](const auto &impl) {
                const auto &map = impl.blocking.cmap();
                std::vector<int> dims;
                for (auto it = map.template attributes<1>().begin(), itend = map.template attributes<1>().end();
                     it != itend;
                     ++it) {
                    dims.push_back(classification_dim(it->info().geom_targets));
                }
                return dims;
            },
            m_impl);
    }

    std::vector<int> BlockingFacade::face_classification_dims() const {
        return std::visit(
            [](const auto &impl) {
                const auto &map = impl.blocking.cmap();
                std::vector<int> dims;
                for (auto it = map.template attributes<2>().begin(), itend = map.template attributes<2>().end();
                     it != itend;
                     ++it) {
                    dims.push_back(classification_dim(it->info().geom_targets));
                }
                return dims;
            },
            m_impl);
    }

    std::vector<std::array<double, 3>> BlockingFacade::edge_control_points() const {
        return std::visit(
            [](const auto &impl) {
                const auto &map = impl.blocking.cmap();
                std::vector<std::array<double, 3>> points;
                for (auto it = map.template attributes<1>().begin(), itend = map.template attributes<1>().end();
                     it != itend;
                     ++it) {
                    for (const auto &cp : it->info().curve.control_points()) {
                        points.push_back({cp.x(), cp.y(), cp.z()});
                    }
                }
                return points;
            },
            m_impl);
    }

    std::vector<std::array<int, 2>> BlockingFacade::edge_control_polygons() const {
        return std::visit(
            [](const auto &impl) {
                constexpr int n = static_cast<int>(std::decay_t<decltype(impl)>::DEGREE) + 1;
                const auto &map = impl.blocking.cmap();
                std::vector<std::array<int, 2>> segments;
                int base = 0;
                for (auto it = map.template attributes<1>().begin(), itend = map.template attributes<1>().end();
                     it != itend;
                     ++it) {
                    for (int i = 0; i + 1 < n; ++i) {
                        segments.push_back({base + i, base + i + 1});
                    }
                    base += n;
                }
                return segments;
            },
            m_impl);
    }

    std::vector<std::array<double, 3>> BlockingFacade::face_control_points() const {
        return std::visit(
            [](const auto &impl) {
                const auto &map = impl.blocking.cmap();
                std::vector<std::array<double, 3>> points;
                for (auto it = map.template attributes<2>().begin(), itend = map.template attributes<2>().end();
                     it != itend;
                     ++it) {
                    const auto &grid = it->info().surface.control_points();
                    for (const auto &row : grid) {
                        for (const auto &cp : row) {
                            points.push_back({cp.x(), cp.y(), cp.z()});
                        }
                    }
                }
                return points;
            },
            m_impl);
    }

    std::vector<std::array<int, 2>> BlockingFacade::face_control_nets() const {
        return std::visit(
            [](const auto &impl) {
                constexpr int n = std::decay_t<decltype(impl)>::DEGREE + 1;
                const auto &map = impl.blocking.cmap();
                std::vector<std::array<int, 2>> segments;
                int base = 0;
                for (auto it = map.template attributes<2>().begin(), itend = map.template attributes<2>().end();
                     it != itend;
                     ++it) {
                    for (int i = 0; i < n; ++i) {
                        for (int j = 0; j < n; ++j) {
                            const int here = base + i * n + j;
                            if (i + 1 < n) segments.push_back({here, here + n}); // along u
                            if (j + 1 < n) segments.push_back({here, here + 1}); // along v
                        }
                    }
                    base += n * n;
                }
                return segments;
            },
            m_impl);
    }

    std::vector<std::array<double, 3>> BlockingFacade::block_control_points() const {
        return std::visit(
            [](const auto &impl) {
                const auto &map = impl.blocking.cmap();
                std::vector<std::array<double, 3>> points;
                for (auto it = map.template attributes<3>().begin(), itend = map.template attributes<3>().end();
                     it != itend;
                     ++it) {
                    const auto &grid = it->info().volume.control_points();
                    for (const auto &plane : grid) {
                        for (const auto &row : plane) {
                            for (const auto &cp : row) {
                                points.push_back({cp.x(), cp.y(), cp.z()});
                            }
                        }
                    }
                }
                return points;
            },
            m_impl);
    }

    std::vector<std::array<int, 2>> BlockingFacade::block_control_lattices() const {
        return std::visit(
            [](const auto &impl) {
                constexpr int n = std::decay_t<decltype(impl)>::DEGREE + 1;
                const auto &map = impl.blocking.cmap();
                std::vector<std::array<int, 2>> segments;
                int base = 0;
                for (auto it = map.template attributes<3>().begin(), itend = map.template attributes<3>().end();
                     it != itend;
                     ++it) {
                    for (int i = 0; i < n; ++i) {
                        for (int j = 0; j < n; ++j) {
                            for (int k = 0; k < n; ++k) {
                                const int here = base + (i * n + j) * n + k;
                                if (i + 1 < n) segments.push_back({here, here + n * n}); // along u
                                if (j + 1 < n) segments.push_back({here, here + n});     // along v
                                if (k + 1 < n) segments.push_back({here, here + 1});     // along w
                            }
                        }
                    }
                    base += n * n * n;
                }
                return segments;
            },
            m_impl);
    }

    std::vector<std::array<double, 3>> BlockingFacade::face_grid_vertices(int subdivisions) const {
        check_subdivisions(subdivisions, "Blocking.face_grid_vertices");
        return std::visit(
            [subdivisions](const auto &impl) {
                const auto &map = impl.blocking.cmap();
                std::vector<std::array<double, 3>> points;
                for (auto it = map.template attributes<2>().begin(), itend = map.template attributes<2>().end();
                     it != itend;
                     ++it) {
                    // Straight off the face's own stored surface, rather than out of to_mesh():
                    // to_mesh() only emits quads for standalone 2D blocks, so a hex block's own
                    // bounding faces would otherwise never be drawable.
                    for (int i = 0; i <= subdivisions; ++i) {
                        for (int j = 0; j <= subdivisions; ++j) {
                            const double u = static_cast<double>(i) / static_cast<double>(subdivisions);
                            const double v = static_cast<double>(j) / static_cast<double>(subdivisions);
                            const auto p = it->info().surface.value(u, v);
                            points.push_back({p.x(), p.y(), p.z()});
                        }
                    }
                }
                return points;
            },
            m_impl);
    }

    std::vector<std::array<int, 4>> BlockingFacade::face_grid_quads(int subdivisions) const {
        check_subdivisions(subdivisions, "Blocking.face_grid_quads");
        return std::visit(
            [subdivisions](const auto &impl) {
                const auto &map = impl.blocking.cmap();
                const int side = subdivisions + 1;
                std::vector<std::array<int, 4>> quads;
                int base = 0;
                for (auto it = map.template attributes<2>().begin(), itend = map.template attributes<2>().end();
                     it != itend;
                     ++it) {
                    for (int i = 0; i < subdivisions; ++i) {
                        for (int j = 0; j < subdivisions; ++j) {
                            const int a = base + i * side + j;
                            quads.push_back({a, a + side, a + side + 1, a + 1});
                        }
                    }
                    base += side * side;
                }
                return quads;
            },
            m_impl);
    }

    std::vector<int> BlockingFacade::face_grid_owners(int subdivisions) const {
        check_subdivisions(subdivisions, "Blocking.face_grid_owners");
        return std::visit(
            [subdivisions](const auto &impl) {
                const auto &map = impl.blocking.cmap();
                const int per_face = subdivisions * subdivisions;
                std::vector<int> owners;
                int face_index = 0;
                for (auto it = map.template attributes<2>().begin(), itend = map.template attributes<2>().end();
                     it != itend;
                     ++it) {
                    owners.insert(owners.end(), static_cast<std::size_t>(per_face), face_index);
                    ++face_index;
                }
                return owners;
            },
            m_impl);
    }

    std::vector<std::array<double, 3>> BlockingFacade::edge_vertices(int samples) const {
        check_subdivisions(samples, "Blocking.edge_vertices");
        return std::visit(
            [samples](const auto &impl) {
                const auto &map = impl.blocking.cmap();
                std::vector<std::array<double, 3>> points;
                for (auto it = map.template attributes<1>().begin(), itend = map.template attributes<1>().end();
                     it != itend;
                     ++it) {
                    for (int i = 0; i <= samples; ++i) {
                        const double t = static_cast<double>(i) / static_cast<double>(samples);
                        const auto p = it->info().curve.value(t);
                        points.push_back({p.x(), p.y(), p.z()});
                    }
                }
                return points;
            },
            m_impl);
    }

    std::vector<std::array<int, 2>> BlockingFacade::edge_segments(int samples) const {
        check_subdivisions(samples, "Blocking.edge_segments");
        return std::visit(
            [samples](const auto &impl) {
                const auto &map = impl.blocking.cmap();
                std::vector<std::array<int, 2>> segments;
                int base = 0;
                for (auto it = map.template attributes<1>().begin(), itend = map.template attributes<1>().end();
                     it != itend;
                     ++it) {
                    for (int i = 0; i < samples; ++i) {
                        segments.push_back({base + i, base + i + 1});
                    }
                    base += samples + 1;
                }
                return segments;
            },
            m_impl);
    }

    std::vector<std::array<double, 3>> BlockingFacade::mesh_vertices(int subdivisions) {
        check_subdivisions(subdivisions, "Blocking.mesh_vertices");
        return std::visit(
            [subdivisions](auto &impl) {
                const auto mesh = impl.blocking.to_mesh(static_cast<SizeT>(subdivisions));
                std::vector<std::array<double, 3>> vertices;
                vertices.reserve(mesh.nb_nodes());
                for (UInt i = 0; i < mesh.nb_nodes(); ++i) {
                    const auto &p = mesh.node(NodeId{i});
                    vertices.push_back({p.x(), p.y(), p.z()});
                }
                return vertices;
            },
            m_impl);
    }

    std::vector<std::array<int, 4>> BlockingFacade::mesh_quads(int subdivisions) {
        check_subdivisions(subdivisions, "Blocking.mesh_quads");
        return std::visit(
            [subdivisions](auto &impl) {
                const auto mesh = impl.blocking.to_mesh(static_cast<SizeT>(subdivisions));
                std::vector<std::array<int, 4>> quads;
                quads.reserve(mesh.nb_faces());
                for (UInt i = 0; i < mesh.nb_faces(); ++i) {
                    const auto n = mesh.face_nodes(FaceId{i});
                    quads.push_back({static_cast<int>(n[0].value),
                                     static_cast<int>(n[1].value),
                                     static_cast<int>(n[2].value),
                                     static_cast<int>(n[3].value)});
                }
                return quads;
            },
            m_impl);
    }

    std::vector<std::array<int, 8>> BlockingFacade::mesh_hexes(int subdivisions) {
        check_subdivisions(subdivisions, "Blocking.mesh_hexes");
        return std::visit(
            [subdivisions](auto &impl) {
                const auto mesh = impl.blocking.to_mesh(static_cast<SizeT>(subdivisions));
                std::vector<std::array<int, 8>> hexes;
                hexes.reserve(mesh.nb_cells());
                for (UInt i = 0; i < mesh.nb_cells(); ++i) {
                    const auto n = mesh.cell_nodes(CellId{i});
                    std::array<int, 8> hex{};
                    for (std::size_t k = 0; k < 8; ++k)
                        hex[k] = static_cast<int>(n[k].value);
                    hexes.push_back(hex);
                }
                return hexes;
            },
            m_impl);
    }

} // namespace gecko::python
