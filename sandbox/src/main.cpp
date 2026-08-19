// Scratch playground for local development — not part of the public API, not installed, not
// covered by lint/doc-completeness checks. Edit freely: builds a couple of sewn hex blocks, meshes
// them, and writes a VTK file you can open in ParaView, as a working starting point to hack on.
#include <array>
#include <iostream>
#include <string>

#include <gecko/block/Blocking.h>
#include <gecko/geom/FacetedGeometry.h>
#include <gecko/io/GmshMeshWriter.h>
#include <gecko/io/VtkMeshWriter.h>

using namespace gecko;

namespace {
    // A minimal in-code geometric model (a single tagged triangle), built via the same
    // SimplicialMesh -> write -> FacetedGeometry read-back pattern used throughout block/tst.
    // Swap this out for `FacetedGeometry geom("path/to/model.msh");` to load a real model instead.
    FacetedGeometry make_demo_geom_model() {
        SimplicialMesh mesh;
        auto n0 = mesh.add_node(0, 0, 0);
        auto n1 = mesh.add_node(1, 0, 0);
        auto n2 = mesh.add_node(0, 1, 0);

        GroupRegistry groups;
        auto surf = groups.add_group("Surf", GroupDim::Dim2);
        auto &face_group = mesh.add_variable<GroupId, CellType::Face>(std::string(io::PHYSICAL_GROUP_VARIABLE));
        auto &face_entity = mesh.add_variable<Int, CellType::Face>(std::string(io::ENTITY_TAG_VARIABLE));
        auto f0 = mesh.add_face(n0, n1, n2);
        face_group[f0.value] = surf;
        face_entity[f0.value] = 1;

        const std::string path = "gecko_sandbox_geom.msh";
        io::SimplicialMeshWriter::write(path, mesh, groups);
        return FacetedGeometry(path);
    }
} // namespace

int main() {
    const FacetedGeometry geom = make_demo_geom_model();
    Blocking<FacetedGeometry> blocking(geom);

    const std::array<Point3d, 8> corners_a = {Point3d(0, 0, 0),
                                              Point3d(1, 0, 0),
                                              Point3d(1, 1, 0),
                                              Point3d(0, 1, 0),
                                              Point3d(0, 0, 1),
                                              Point3d(1, 0, 1),
                                              Point3d(1, 1, 1),
                                              Point3d(0, 1, 1)};
    const std::array<Point3d, 8> corners_b = {Point3d(1, 0, 0),
                                              Point3d(2, 0, 0),
                                              Point3d(2, 1, 0),
                                              Point3d(1, 1, 0),
                                              Point3d(1, 0, 1),
                                              Point3d(2, 0, 1),
                                              Point3d(2, 1, 1),
                                              Point3d(1, 1, 1)};
    blocking.create_hex_block(corners_a);
    blocking.create_hex_block(corners_b);
    blocking.build_connectivity();

    std::cout << "Blocking: " << blocking.nb_cells<0>() << " nodes, " << blocking.nb_cells<3>()
              << " blocks, valid=" << std::boolalpha << blocking.is_valid_topology() << "\n";

    auto mesh = blocking.to_mesh(4);
    std::cout << "Mesh: " << mesh.nb_nodes() << " nodes, " << mesh.nb_cells() << " hexes\n";

    const std::string out_path = "gecko_sandbox_output.vtk";
    io::CubicVtkWriter::write(out_path, mesh);
    std::cout << "Wrote " << out_path << " -- open it in ParaView to see the result.\n";

    return 0;
}
