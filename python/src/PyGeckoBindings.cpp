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
                 "Node index quadruples of every tetrahedron of the model's backing mesh (empty for a boundary rep).");

        py::class_<BlockingFacade>(m, "Blocking", "A structured (quad/hex) blocking of a GeomModel.")
            .def(py::init<const GeomModelFacade &, int>(),
                 py::arg("model"),
                 py::arg("degree") = 1,
                 py::keep_alive<1, 2>(),
                 "Builds an empty blocking against model (degree=1: straight edges, degree=3: cubic Bezier edges). "
                 "model must be kept alive for as long as this Blocking is used.")
            .def("create_quad_block",
                 &BlockingFacade::create_quad_block,
                 py::arg("corners"),
                 "Creates a standalone quad block from its 4 (x,y,z) corners, in perimeter order. Returns the new "
                 "face id.")
            .def("create_hex_block",
                 &BlockingFacade::create_hex_block,
                 py::arg("corners"),
                 "Creates a standalone hex block from its 8 (x,y,z) corners (HEX8 ordering). Returns the new block id.")
            .def("build_connectivity",
                 &BlockingFacade::build_connectivity,
                 "Auto-detects and sews coincident blocks created so far. Not incremental.")
            .def("classify",
                 &BlockingFacade::classify,
                 py::arg("tol_vertex"),
                 py::arg("tol_curve_surface") = -1.0,
                 "Classifies every cell onto the geometric model and refits geometry accordingly.")
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
                 "it as a VTK legacy ASCII file.");
    }

} // namespace gecko::python
