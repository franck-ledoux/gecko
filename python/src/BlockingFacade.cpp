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

        void check_degree(int degree, const char *who) {
            if (degree < BlockingFacade::MIN_DEGREE) {
                throw std::invalid_argument(std::string(who) + ": degree must be at least " +
                                            std::to_string(BlockingFacade::MIN_DEGREE) + " (1 = straight), got " +
                                            std::to_string(degree));
            }
        }
    } // namespace

    namespace {
        /** @brief Validates a degree and hands it back, so it can be checked inside a member
         * initializer list. */
        int checked_degree(int degree) {
            check_degree(degree, "Blocking");
            return degree;
        }
    } // namespace

    BlockingFacade::BlockingFacade(const GeomModelFacade &model, int degree)
        : m_impl(model.native(), checked_degree(degree)) {}

    int BlockingFacade::degree() const { return static_cast<int>(m_impl.blocking.degree()); }

    void BlockingFacade::set_degree(int degree, double tol_vertex, double tol_curve, double tol_surface) {
        check_degree(degree, "Blocking.set_degree");
        m_impl.blocking.set_degree(static_cast<std::size_t>(degree), tol_vertex, tol_curve, tol_surface);
    }

    int BlockingFacade::create_quad_block(const std::vector<std::array<double, 3>> &corners) {
        const auto points = to_points<4>(corners, "Blocking.create_quad_block");
        return static_cast<int>(m_impl.blocking.create_quad_block(points)->info().id);
    }

    int BlockingFacade::create_hex_block(const std::vector<std::array<double, 3>> &corners) {
        const auto points = to_points<8>(corners, "Blocking.create_hex_block");
        return static_cast<int>(m_impl.blocking.create_hex_block(points)->info().id);
    }

    void BlockingFacade::build_connectivity() { m_impl.blocking.build_connectivity(); }

    void BlockingFacade::classify(double tol_vertex, double tol_curve, double tol_surface) {
        m_impl.blocking.classify(tol_vertex, tol_curve, tol_surface);
    }

    std::size_t BlockingFacade::nb_cells(int dim) const {
        switch (dim) {
            case 0:
                return m_impl.blocking.nb_cells<0>();
            case 1:
                return m_impl.blocking.nb_cells<1>();
            case 2:
                return m_impl.blocking.nb_cells<2>();
            case 3:
                return m_impl.blocking.nb_cells<3>();
            default:
                throw std::invalid_argument("Blocking.nb_cells: dim must be in [0,3], got " + std::to_string(dim));
        }
    }

    bool BlockingFacade::is_valid_topology() const { return m_impl.blocking.is_valid_topology(); }

    bool BlockingFacade::is_purely_2d() const { return m_impl.blocking.is_purely_2d(); }

    namespace {
        /** @brief The cell of dimension @p TDim carrying @p id, or a throw naming what was not found.
         * @tparam TDim The cell dimension.
         * @tparam TImpl The blocking state type, const or not.
         * @param impl The blocking state.
         * @param id The id to resolve.
         * @param what What to call it in the message ("edge", "face", "block").
         * @return That cell.
         * @throw std::out_of_range if no cell of that dimension carries @p id — which is what a
         *        caller sees when the cell it was holding on to has since been deleted. */
        template<unsigned int TDim, typename TImpl>
        auto cell_or_throw(TImpl &impl, int id, const char *what) {
            const auto cell = impl.blocking.template cell_by_id<TDim>(id);
            if (cell == nullptr) {
                throw std::out_of_range(std::string("Blocking: unknown ") + what + " id " + std::to_string(id));
            }
            return cell;
        }
    } // namespace

    bool BlockingFacade::can_delete_face(int face_id) const {
        return m_impl.blocking.can_delete_face(cell_or_throw<2>(m_impl, face_id, "face"));
    }

    void BlockingFacade::delete_face(int face_id) {
        m_impl.blocking.delete_face(cell_or_throw<2>(m_impl, face_id, "face"));
    }

    void BlockingFacade::delete_block(int block_id) {
        m_impl.blocking.delete_block(cell_or_throw<3>(m_impl, block_id, "block"));
    }

    std::vector<double> BlockingFacade::edge_bends() const {
        const auto &map = m_impl.blocking.cmap();
        std::vector<double> bends;
        for (auto it = map.attributes<1>().begin(), itend = map.attributes<1>().end(); it != itend; ++it) {
            const auto &points = it->info().curve.control_points();
            const Point3d &start = points.front();
            const Point3d &end = points.back();
            const Vector3d chord(start, end);
            const double length2 = chord.dot(chord);
            double worst = 0.0;
            for (std::size_t k = 1; k + 1 < points.size(); ++k) {
                const Vector3d offset(start, points[k]);
                // Distance to the chord's *line*, clamped to the segment: an interior control
                // point of a straight edge projects inside it, and one that does not is
                // already the anomaly this is looking for.
                const double t = (length2 > 0.0) ? offset.dot(chord) / length2 : 0.0;
                worst = std::max(worst, (offset - chord * t).norm());
            }
            bends.push_back(worst);
        }
        return bends;
    }

    void BlockingFacade::write_vtk(int subdivisions, const std::string &path) {
        check_subdivisions(subdivisions, "Blocking.write_vtk");
        const auto mesh = m_impl.blocking.to_mesh(static_cast<SizeT>(subdivisions));
        io::VtkMeshWriter<CubicTraits>::write(path,
                                              mesh,
                                              {std::string(Impl::BlockingT::NODE_CLASSIFICATION_DIM_VARIABLE),
                                               std::string(Impl::BlockingT::NODE_CLASSIFICATION_TAG_VARIABLE)});
    }

    std::vector<double> BlockingFacade::block_volumes(int subdivisions) {
        check_subdivisions(subdivisions, "Blocking.block_volumes");
        return m_impl.blocking.block_volumes(subdivisions);
    }

    std::vector<int> BlockingFacade::sheet_edges(int edge_id) {
        const auto sheet = m_impl.blocking.find_sheet(cell_or_throw<1>(m_impl, edge_id, "edge"));
        std::vector<int> indices;
        if (!sheet.has_value()) return indices;

        // Back from attribute handles to the positions the caller speaks in, in one pass.
        auto &map = m_impl.blocking.cmap();
        int position = 0;
        for (auto it = map.attributes<1>().begin(), itend = map.attributes<1>().end(); it != itend; ++it, ++position) {
            for (const auto &se : sheet->edges) {
                if (se.edge == it) {
                    indices.push_back(position);
                    break;
                }
            }
        }
        return indices;
    }

    std::vector<std::array<double, 3>> BlockingFacade::sheet_cut_points(int edge_id, double param) {
        std::vector<std::array<double, 3>> points;
        const auto sheet = m_impl.blocking.find_sheet(cell_or_throw<1>(m_impl, edge_id, "edge"));
        if (!sheet.has_value()) return points;

        // Walked in the blocking's own edge order rather than the sheet's: find_sheet()
        // returns its edges keyed by attribute handle, whose ordering is the container's and
        // need not match the traversal sheet_edges() reports positions in. Re-walking is what
        // makes the 2 line up entry for entry, which is what a caller pairing them expects.
        auto &map = m_impl.blocking.cmap();
        points.reserve(sheet->edges.size());
        for (auto it = map.attributes<1>().begin(), itend = map.attributes<1>().end(); it != itend; ++it) {
            for (const auto &se : sheet->edges) {
                if (se.edge != it) continue;
                const auto p = m_impl.blocking.cut_point(se, param);
                points.push_back({p.x(), p.y(), p.z()});
                break;
            }
        }
        return points;
    }

    bool BlockingFacade::cut_sheet(int edge_id, double param) {
        if (!m_impl.blocking.cut_sheet(cell_or_throw<1>(m_impl, edge_id, "edge"), param)) return false;
        return true;
    }

    bool BlockingFacade::delete_sheet(int edge_id, double tol_vertex, double tol_curve, double tol_surface) {
        if (!m_impl.blocking.delete_sheet(
                cell_or_throw<1>(m_impl, edge_id, "edge"), tol_vertex, tol_curve, tol_surface)) {
            return false;
        }
        // Collapsing merges corners pairwise, so ids naming the side that went have to go — and the
        // survivors can come back as attributes CGAL rebuilt, which no id then names. Prune, then
        // re-index, for the reason delete_block() spells out.
        return true;
    }

    namespace {
        /** @brief The node carrying @p node_id, straight from the map.
         *
         * There used to be a book of ids kept alongside the blocking here, refreshed after every
         * operation: pruned of the corners a deletion had collected, then extended with the ones an
         * edit had created. The kernel names its own cells now (`CellInfo::id`), so the book was
         * one more thing that had to be right rather than a source of truth.
         */
        template<typename TImpl>
        auto find_node(TImpl &impl, int node_id) {
            const auto node = impl.blocking.template cell_by_id<0>(node_id);
            if (node == nullptr) {
                throw std::out_of_range("Blocking: unknown node id " + std::to_string(node_id));
            }
            return node;
        }
    } // namespace

    namespace {
        /** @brief Every cell id of one dimension, in the map's own traversal order.
         * @tparam TDim The cell dimension.
         * @param map The blocking's map.
         * @return One id per cell, in the order the display accessors of that dimension emit. */
        template<unsigned int TDim, typename TMap>
        std::vector<int> ids_in_traversal_order(const TMap &map) {
            std::vector<int> ids;
            for (auto it = map.template attributes<TDim>().begin(), itend = map.template attributes<TDim>().end();
                 it != itend;
                 ++it) {
                ids.push_back(static_cast<int>(it->info().id));
            }
            return ids;
        }
    } // namespace

    std::vector<int> BlockingFacade::edge_ids() const { return ids_in_traversal_order<1>(m_impl.blocking.cmap()); }

    std::vector<int> BlockingFacade::face_ids() const { return ids_in_traversal_order<2>(m_impl.blocking.cmap()); }

    std::vector<int> BlockingFacade::block_ids() const { return ids_in_traversal_order<3>(m_impl.blocking.cmap()); }

    std::vector<int> BlockingFacade::node_ids() const {
        const auto &map = m_impl.blocking.cmap();
        std::vector<int> ids;
        for (auto it = map.attributes<0>().begin(), itend = map.attributes<0>().end(); it != itend; ++it) {
            ids.push_back(static_cast<int>(it->info().id));
        }
        std::sort(ids.begin(), ids.end());
        return ids;
    }

    std::array<double, 3> BlockingFacade::node_position(int node_id) const {
        const auto &p = find_node(m_impl, node_id)->info().point;
        return std::array<double, 3>{p.x(), p.y(), p.z()};
    }

    void BlockingFacade::move_node(int node_id, double x, double y, double z) {
        m_impl.blocking.move_node(find_node(m_impl, node_id), Point3d(x, y, z));
    }

    void BlockingFacade::snap_node(int node_id, double tol_vertex, double tol_curve, double tol_surface) {
        m_impl.blocking.snap_node(find_node(m_impl, node_id), tol_vertex, tol_curve, tol_surface);
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
        const auto &map = m_impl.blocking.cmap();
        // Paired with the id and sorted on it, so the result lines up with node_ids() whatever order
        // the map hands its attributes back in.
        std::vector<std::pair<int, int>> by_id;
        for (auto it = map.attributes<0>().begin(), itend = map.attributes<0>().end(); it != itend; ++it) {
            by_id.emplace_back(static_cast<int>(it->info().id), classification_dim(it->info().geom_targets));
        }
        std::sort(by_id.begin(), by_id.end());

        std::vector<int> dims;
        dims.reserve(by_id.size());
        for (const auto &[id, dim] : by_id) {
            dims.push_back(dim);
        }
        return dims;
    }

    std::vector<int> BlockingFacade::edge_classification_dims() const {
        const auto &map = m_impl.blocking.cmap();
        std::vector<int> dims;
        for (auto it = map.attributes<1>().begin(), itend = map.attributes<1>().end(); it != itend; ++it) {
            dims.push_back(classification_dim(it->info().geom_targets));
        }
        return dims;
    }

    std::vector<int> BlockingFacade::face_classification_dims() const {
        const auto &map = m_impl.blocking.cmap();
        std::vector<int> dims;
        for (auto it = map.attributes<2>().begin(), itend = map.attributes<2>().end(); it != itend; ++it) {
            dims.push_back(classification_dim(it->info().geom_targets));
        }
        return dims;
    }

    std::vector<std::array<double, 3>> BlockingFacade::edge_control_points() const {
        const auto &map = m_impl.blocking.cmap();
        std::vector<std::array<double, 3>> points;
        for (auto it = map.attributes<1>().begin(), itend = map.attributes<1>().end(); it != itend; ++it) {
            for (const auto &cp : it->info().curve.control_points()) {
                points.push_back({cp.x(), cp.y(), cp.z()});
            }
        }
        return points;
    }

    std::vector<std::array<int, 2>> BlockingFacade::edge_control_polygons() const {
        const int n = static_cast<int>(m_impl.blocking.degree()) + 1;
        const auto &map = m_impl.blocking.cmap();
        std::vector<std::array<int, 2>> segments;
        int base = 0;
        for (auto it = map.attributes<1>().begin(), itend = map.attributes<1>().end(); it != itend; ++it) {
            for (int i = 0; i + 1 < n; ++i) {
                segments.push_back({base + i, base + i + 1});
            }
            base += n;
        }
        return segments;
    }

    std::vector<std::array<double, 3>> BlockingFacade::face_control_points() const {
        const auto &map = m_impl.blocking.cmap();
        std::vector<std::array<double, 3>> points;
        for (auto it = map.attributes<2>().begin(), itend = map.attributes<2>().end(); it != itend; ++it) {
            const auto &grid = it->info().surface.control_points();
            for (const auto &row : grid) {
                for (const auto &cp : row) {
                    points.push_back({cp.x(), cp.y(), cp.z()});
                }
            }
        }
        return points;
    }

    std::vector<std::array<int, 2>> BlockingFacade::face_control_nets() const {
        const int n = static_cast<int>(m_impl.blocking.degree()) + 1;
        const auto &map = m_impl.blocking.cmap();
        std::vector<std::array<int, 2>> segments;
        int base = 0;
        for (auto it = map.attributes<2>().begin(), itend = map.attributes<2>().end(); it != itend; ++it) {
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
    }

    std::vector<std::array<double, 3>> BlockingFacade::block_control_points() const {
        const auto &map = m_impl.blocking.cmap();
        std::vector<std::array<double, 3>> points;
        for (auto it = map.attributes<3>().begin(), itend = map.attributes<3>().end(); it != itend; ++it) {
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
    }

    std::vector<std::array<int, 2>> BlockingFacade::block_control_lattices() const {
        const int n = static_cast<int>(m_impl.blocking.degree()) + 1;
        const auto &map = m_impl.blocking.cmap();
        std::vector<std::array<int, 2>> segments;
        int base = 0;
        for (auto it = map.attributes<3>().begin(), itend = map.attributes<3>().end(); it != itend; ++it) {
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
    }

    std::vector<std::array<double, 3>> BlockingFacade::face_grid_vertices(int subdivisions) const {
        check_subdivisions(subdivisions, "Blocking.face_grid_vertices");
        const auto &map = m_impl.blocking.cmap();
        std::vector<std::array<double, 3>> points;
        for (auto it = map.attributes<2>().begin(), itend = map.attributes<2>().end(); it != itend; ++it) {
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
    }

    std::vector<std::array<int, 4>> BlockingFacade::face_grid_quads(int subdivisions) const {
        check_subdivisions(subdivisions, "Blocking.face_grid_quads");
        const auto &map = m_impl.blocking.cmap();
        const int side = subdivisions + 1;
        std::vector<std::array<int, 4>> quads;
        int base = 0;
        for (auto it = map.attributes<2>().begin(), itend = map.attributes<2>().end(); it != itend; ++it) {
            for (int i = 0; i < subdivisions; ++i) {
                for (int j = 0; j < subdivisions; ++j) {
                    const int a = base + i * side + j;
                    quads.push_back({a, a + side, a + side + 1, a + 1});
                }
            }
            base += side * side;
        }
        return quads;
    }

    std::vector<int> BlockingFacade::face_grid_owners(int subdivisions) const {
        check_subdivisions(subdivisions, "Blocking.face_grid_owners");
        const auto &map = m_impl.blocking.cmap();
        const int per_face = subdivisions * subdivisions;
        std::vector<int> owners;
        int face_index = 0;
        for (auto it = map.attributes<2>().begin(), itend = map.attributes<2>().end(); it != itend; ++it) {
            owners.insert(owners.end(), static_cast<std::size_t>(per_face), face_index);
            ++face_index;
        }
        return owners;
    }

    std::vector<std::array<double, 3>> BlockingFacade::edge_vertices(int samples) const {
        check_subdivisions(samples, "Blocking.edge_vertices");
        const auto &map = m_impl.blocking.cmap();
        std::vector<std::array<double, 3>> points;
        for (auto it = map.attributes<1>().begin(), itend = map.attributes<1>().end(); it != itend; ++it) {
            for (int i = 0; i <= samples; ++i) {
                const double t = static_cast<double>(i) / static_cast<double>(samples);
                const auto p = it->info().curve.value(t);
                points.push_back({p.x(), p.y(), p.z()});
            }
        }
        return points;
    }

    std::vector<std::array<int, 2>> BlockingFacade::edge_segments(int samples) const {
        check_subdivisions(samples, "Blocking.edge_segments");
        const auto &map = m_impl.blocking.cmap();
        std::vector<std::array<int, 2>> segments;
        int base = 0;
        for (auto it = map.attributes<1>().begin(), itend = map.attributes<1>().end(); it != itend; ++it) {
            for (int i = 0; i < samples; ++i) {
                segments.push_back({base + i, base + i + 1});
            }
            base += samples + 1;
        }
        return segments;
    }

    std::vector<std::array<double, 3>> BlockingFacade::mesh_vertices(int subdivisions) {
        check_subdivisions(subdivisions, "Blocking.mesh_vertices");
        const auto mesh = m_impl.blocking.to_mesh(static_cast<SizeT>(subdivisions));
        std::vector<std::array<double, 3>> vertices;
        vertices.reserve(mesh.nb_nodes());
        for (UInt i = 0; i < mesh.nb_nodes(); ++i) {
            const auto &p = mesh.node(NodeId{i});
            vertices.push_back({p.x(), p.y(), p.z()});
        }
        return vertices;
    }

    std::vector<std::array<int, 4>> BlockingFacade::mesh_quads(int subdivisions) {
        check_subdivisions(subdivisions, "Blocking.mesh_quads");
        const auto mesh = m_impl.blocking.to_mesh(static_cast<SizeT>(subdivisions));
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
    }

    std::vector<std::array<int, 8>> BlockingFacade::mesh_hexes(int subdivisions) {
        check_subdivisions(subdivisions, "Blocking.mesh_hexes");
        const auto mesh = m_impl.blocking.to_mesh(static_cast<SizeT>(subdivisions));
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
    }

    std::vector<int> BlockingFacade::mesh_hex_owners(int subdivisions) {
        check_subdivisions(subdivisions, "Blocking.mesh_hex_owners");
        // to_mesh() walks the blocks in their own traversal order and emits exactly one grid of
        // subdivisions^3 cells per block, so the owner is a matter of counting rather than of asking
        // the kernel. Stated here, once, so that whoever colours a mesh by block does not have to
        // know it — and so a change to that emission order fails a test rather than silently
        // recolouring the picture.
        const auto per_block = static_cast<std::size_t>(subdivisions) * static_cast<std::size_t>(subdivisions) *
                               static_cast<std::size_t>(subdivisions);
        const auto &map = m_impl.blocking.cmap();
        std::vector<int> owners;
        int block_index = 0;
        for (auto it = map.attributes<3>().begin(), itend = map.attributes<3>().end(); it != itend; ++it) {
            owners.insert(owners.end(), per_block, block_index);
            ++block_index;
        }
        return owners;
    }

    std::vector<int> BlockingFacade::mesh_quad_owners(int subdivisions) {
        check_subdivisions(subdivisions, "Blocking.mesh_quad_owners");
        // Same counting, one dimension down — except that only a face belonging to no block emits
        // mesh quads at all: a hex's own bounding faces are part of its cells, not quads of their
        // own. So the faces are walked in the same order and the ones a block owns are skipped,
        // exactly as to_mesh() skips them.
        const auto per_face = static_cast<std::size_t>(subdivisions) * static_cast<std::size_t>(subdivisions);
        const auto &map = m_impl.blocking.cmap();
        std::vector<int> owners;
        int quad_block_index = 0;
        for (auto it = map.attributes<2>().begin(), itend = map.attributes<2>().end(); it != itend; ++it) {
            if (map.attribute<3>(it->dart()) != nullptr) continue;
            owners.insert(owners.end(), per_face, quad_block_index);
            ++quad_block_index;
        }
        return owners;
    }

} // namespace gecko::python
