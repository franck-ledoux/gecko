#include "PyGeckoBindings.h"

#include <string>

#include <pybind11/stl.h>

#include "BlockingFacade.h"
#include "GeomModelFacade.h"

namespace py = pybind11;

namespace gecko::python {

    namespace {
        std::string hello() { return "Hello from gecko!"; }
    } // namespace

    void register_gecko_module(py::module_ m) {
        m.doc() = "Gecko Python bindings: façade managers over the gecko C++ API.";
        m.def("hello", &hello, "Returns a greeting string, confirming the gecko extension loads and runs.");

        py::class_<GeomModelFacade>(m, "GeomModel", "A geometric model, loaded from a Gmsh MSH file.")
            .def(py::init<const std::string &>(),
                 py::arg("path"),
                 "Loads the model from a Gmsh MSH file (triangles for a boundary rep, tetrahedra to also get volumes).")
            .def("nb_vertices", &GeomModelFacade::nb_vertices, "Number of vertices in the model.")
            .def("nb_curves", &GeomModelFacade::nb_curves, "Number of curves in the model.")
            .def("nb_surfaces", &GeomModelFacade::nb_surfaces, "Number of surfaces in the model.")
            .def("nb_volumes", &GeomModelFacade::nb_volumes, "Number of volumes in the model.")
            .def("vertex_tags", &GeomModelFacade::vertex_tags, "Entity tags of every vertex in the model.")
            .def("curve_tags", &GeomModelFacade::curve_tags, "Entity tags of every curve in the model.")
            .def("surface_tags", &GeomModelFacade::surface_tags, "Entity tags of every surface in the model.")
            .def("volume_tags", &GeomModelFacade::volume_tags, "Entity tags of every volume in the model.")
            .def("group_ids", &GeomModelFacade::group_ids, "Ids of every physical group defined on the model.")
            .def("group_name", &GeomModelFacade::group_name, py::arg("id"), "Name of the physical group with this id.")
            .def("group_dim",
                 &GeomModelFacade::group_dim,
                 py::arg("id"),
                 "Dimension (0..3) the physical group with this id applies to.")
            .def("group_entities",
                 &GeomModelFacade::group_entities,
                 py::arg("id"),
                 "(dimension, entity_tag) pairs of every entity belonging to the physical group with this id.")
            .def("mesh_vertices",
                 &GeomModelFacade::mesh_vertices,
                 "(x,y,z) position of every node of the model's backing faceted mesh.")
            .def("mesh_triangles",
                 &GeomModelFacade::mesh_triangles,
                 "Node index triples of every triangle of the model's backing mesh.")
            .def("mesh_tets",
                 &GeomModelFacade::mesh_tets,
                 "Node index quadruples of every tetrahedron of the model's backing mesh (empty for a boundary rep).")
            .def("volume_boundary_triangles",
                 &GeomModelFacade::volume_boundary_triangles,
                 "Node index triples of the tetrahedral mesh's own outer boundary — its skin, not its interior tets.")
            .def("curve_vertices",
                 &GeomModelFacade::curve_vertices,
                 "(x,y,z) points the model's curves pass through, indexed by curve_segments().")
            .def("curve_segments",
                 &GeomModelFacade::curve_segments,
                 "Index pairs into curve_vertices() for every segment making up the model's curves.")
            .def("vertex_positions",
                 &GeomModelFacade::vertex_positions,
                 "(x,y,z) position of each of the model's vertices, in vertex_tags() order.");

        py::class_<BlockingFacade>(m, "Blocking", "A structured (quad/hex) blocking of a GeomModel.")
            .def(py::init<const GeomModelFacade &, int>(),
                 py::arg("model"),
                 py::arg("degree") = 1,
                 py::keep_alive<1, 2>(),
                 "Builds an empty blocking against model (degree=1: straight edges, degree=3: cubic Bezier edges). "
                 "The degree can be changed later with set_degree(). model must be kept alive for as long as this "
                 "Blocking is used.")
            .def("create_quad_block",
                 &BlockingFacade::create_quad_block,
                 py::arg("corners"),
                 "Creates a standalone quad block from its 4 (x,y,z) corners, in perimeter order. Returns the new "
                 "face id.")
            .def("create_hex_block",
                 &BlockingFacade::create_hex_block,
                 py::arg("corners"),
                 "Creates a standalone hex block from its 8 (x,y,z) corners (HEX8 ordering). Returns where it sits in "
                 "the block traversal — a position, not a lasting id: cutting or deleting renumbers what follows.")
            .def("build_connectivity",
                 &BlockingFacade::build_connectivity,
                 "Auto-detects and sews coincident blocks created so far. Not incremental.")
            .def("degree", &BlockingFacade::degree, "Edge curve degree this blocking currently uses.")
            .def("set_degree",
                 &BlockingFacade::set_degree,
                 py::arg("degree"),
                 py::arg("tol_vertex"),
                 py::arg("tol_curve") = -1.0,
                 py::arg("tol_surface") = -1.0,
                 "Rebuilds every edge, face and block at a new degree and refits them onto the model. Topology and "
                 "classification are untouched; raising the order lets an edge that lies on a curved model curve "
                 "actually follow it instead of cutting across as a chord.")
            .def("classify",
                 &BlockingFacade::classify,
                 py::arg("tol_vertex"),
                 py::arg("tol_curve") = -1.0,
                 py::arg("tol_surface") = -1.0,
                 "Classifies every cell onto the geometric model and refits geometry accordingly. An omitted "
                 "curve tolerance reuses the vertex one, an omitted surface tolerance the curve one.")
            .def("snap_node",
                 &BlockingFacade::snap_node,
                 py::arg("node_id"),
                 py::arg("tol_vertex"),
                 py::arg("tol_curve") = -1.0,
                 py::arg("tol_surface") = -1.0,
                 "Snaps one corner onto the model, reclassifying and refitting only the cells touching it.")
            .def("nb_cells",
                 &BlockingFacade::nb_cells,
                 py::arg("dim"),
                 "Number of dim-cells (dim in [0,3]) in the blocking.")
            .def("is_valid_topology",
                 &BlockingFacade::is_valid_topology,
                 "Checks the topological validity of the blocking.")
            .def(
                "is_purely_2d", &BlockingFacade::is_purely_2d, "Checks whether the blocking has no 3-cell (hex block).")
            .def("can_delete_face",
                 &BlockingFacade::can_delete_face,
                 py::arg("face_id"),
                 "Checks whether the face with this id can be deleted.")
            .def("delete_face", &BlockingFacade::delete_face, py::arg("face_id"), "Deletes the face with this id.")
            .def("node_ids", &BlockingFacade::node_ids, "Ids of every corner node of the block structure.")
            .def("node_position",
                 &BlockingFacade::node_position,
                 py::arg("node_id"),
                 "(x,y,z) position of the corner node with this id.")
            .def("move_node",
                 &BlockingFacade::move_node,
                 py::arg("node_id"),
                 py::arg("x"),
                 py::arg("y"),
                 py::arg("z"),
                 "Moves the corner node with this id, refitting every incident edge/face/block geometry.")
            .def("node_classification_dims",
                 &BlockingFacade::node_classification_dims,
                 "Dimension of the entity each corner node is classified on (0 vertex, 1 curve, 2 surface, "
                 "3 volume), or -1 if unclassified, in node_ids() order.")
            .def("edge_classification_dims",
                 &BlockingFacade::edge_classification_dims,
                 "Dimension of the entity each block edge is classified on, or -1 if unclassified.")
            .def("face_classification_dims",
                 &BlockingFacade::face_classification_dims,
                 "Dimension of the entity each block face is classified on, or -1 if unclassified.")
            .def("edge_control_points",
                 &BlockingFacade::edge_control_points,
                 "(x,y,z) control points of every block edge, ``degree + 1`` per edge.")
            .def("edge_control_polygons",
                 &BlockingFacade::edge_control_polygons,
                 "Index pairs joining consecutive control points within each edge.")
            .def("face_control_points",
                 &BlockingFacade::face_control_points,
                 "(x,y,z) control points of every block face, ``(degree+1)**2`` per face.")
            .def("face_control_nets",
                 &BlockingFacade::face_control_nets,
                 "Index pairs joining each face control point to its neighbours along u and v.")
            .def("block_control_points",
                 &BlockingFacade::block_control_points,
                 "(x,y,z) control points of every block, ``(degree+1)**3`` per block.")
            .def("block_control_lattices",
                 &BlockingFacade::block_control_lattices,
                 "Index pairs joining each block control point to its neighbours along u, v and w.")
            .def("face_grid_vertices",
                 &BlockingFacade::face_grid_vertices,
                 py::arg("subdivisions") = 1,
                 "(x,y,z) points sampling every block face into a grid, for display.")
            .def("face_grid_quads",
                 &BlockingFacade::face_grid_quads,
                 py::arg("subdivisions") = 1,
                 "Quads joining the points face_grid_vertices() returns.")
            .def("face_grid_owners",
                 &BlockingFacade::face_grid_owners,
                 py::arg("subdivisions") = 1,
                 "Which block face each quad of face_grid_quads() came from.")
            .def("edge_vertices",
                 &BlockingFacade::edge_vertices,
                 py::arg("samples") = 1,
                 "(x,y,z) sample points along every edge of the block structure, ``samples + 1`` per edge.")
            .def("edge_segments",
                 &BlockingFacade::edge_segments,
                 py::arg("samples") = 1,
                 "Index pairs joining the points edge_vertices() returns, in the same order.")
            .def("mesh_vertices",
                 &BlockingFacade::mesh_vertices,
                 py::arg("subdivisions") = 1,
                 "(x,y,z) position of every node of the mesh this blocking generates.")
            .def("mesh_quads",
                 &BlockingFacade::mesh_quads,
                 py::arg("subdivisions") = 1,
                 "Node index quadruples of every quad of the mesh this blocking generates.")
            .def("mesh_hexes",
                 &BlockingFacade::mesh_hexes,
                 py::arg("subdivisions") = 1,
                 "Node index 8-tuples of every hex of the mesh this blocking generates.")
            .def("write_vtk",
                 &BlockingFacade::write_vtk,
                 py::arg("subdivisions"),
                 py::arg("path"),
                 "Generates a mesh reproducing the blocking (subdivided ``subdivisions`` times per axis) and writes "
                 "it as a VTK legacy ASCII file.")
            .def("block_volumes",
                 &BlockingFacade::block_volumes,
                 py::arg("subdivisions") = 1,
                 "Signed volume of every block, from its own stored geometry. Exact at 1 subdivision for a block "
                 "with planar faces; converges as it grows for a curved one. A negative value means that block's "
                 "frame is inverted.")
            .def("edge_bends",
                 &BlockingFacade::edge_bends,
                 "How far each edge departs from a straight line — the largest distance from one of its interior "
                 "control points to its own chord. A straight blocking reads as zero here whatever is done to it.")
            .def("delete_block",
                 &BlockingFacade::delete_block,
                 py::arg("block_index"),
                 "Deletes one block, along with every face, edge and corner that existed only because of it. "
                 "What it shared with a neighbouring block stays, as that neighbour's boundary.")
            .def("sheet_edges",
                 &BlockingFacade::sheet_edges,
                 py::arg("edge_index"),
                 "Every edge ``cut_sheet`` would split if aimed at ``edge_index`` — the whole sheet, as positions in "
                 "the same order ``edge_vertices``/``edge_segments`` use. Empty when the sheet cannot be cut.")
            .def("sheet_cut_points",
                 &BlockingFacade::sheet_cut_points,
                 py::arg("edge_index"),
                 py::arg("param"),
                 "Where a cut would land: one point per sheet edge, in the same order ``sheet_edges`` reports, each "
                 "on the side the whole sheet agrees on. Empty when the sheet cannot be cut.")
            .def("cut_sheet",
                 &BlockingFacade::cut_sheet,
                 py::arg("edge_index"),
                 py::arg("param"),
                 "Cuts the blocking along the whole sheet through ``edge_index``, at ``param`` along that edge "
                 "(strictly between 0 and 1), splitting every block the sheet crosses in 2 and keeping the geometry "
                 "it cut through exactly. Returns False, changing nothing, if the cut is not possible.")
            .def("mesh_hex_owners",
                 &BlockingFacade::mesh_hex_owners,
                 py::arg("subdivisions") = 1,
                 "Which block each hex of ``mesh_hexes`` came from, as a position in the block order "
                 "``block_volumes``/``delete_block`` speak — one entry per hex, so a per-block value can be spread "
                 "onto every cell subdividing it.")
            .def("mesh_quad_owners",
                 &BlockingFacade::mesh_quad_owners,
                 py::arg("subdivisions") = 1,
                 "Which block each quad of ``mesh_quads`` came from, counting standalone quad blocks in the order "
                 "they are emitted. Not an index into ``face_classification_dims``: a hex's own bounding faces "
                 "generate no mesh quads and are not counted.")
            .def("delete_sheet",
                 &BlockingFacade::delete_sheet,
                 py::arg("edge_index"),
                 py::arg("tol_vertex"),
                 py::arg("tol_curve") = -1.0,
                 py::arg("tol_surface") = -1.0,
                 "Deletes the whole sheet through ``edge_index``, collapsing every block it crosses and gluing back "
                 "what was either side of it — the inverse of ``cut_sheet``. Where 2 corners merge, the more "
                 "constrained classification wins (a model vertex over a curve, a curve over a surface) and the "
                 "merged corner is projected onto it. Returns False, changing nothing, if the sheet cannot be "
                 "collapsed — in particular when it is the whole blocking, leaving nothing either side to glue.");
    }

} // namespace gecko::python
