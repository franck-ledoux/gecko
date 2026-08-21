"""Tests for the gecko.GeomModel façade — see docs/user-guide/python.md."""

import pytest

import gecko


def test_counts_and_tags(geom_model_path):
    model = gecko.GeomModel(geom_model_path)

    assert model.nb_vertices() == 1
    assert model.nb_curves() == 0
    assert model.nb_surfaces() == 1
    assert model.nb_volumes() == 0

    assert model.vertex_tags() == [1]
    assert model.curve_tags() == []
    assert model.surface_tags() == [1]
    assert model.volume_tags() == []


def test_groups_round_trip(geom_model_path):
    model = gecko.GeomModel(geom_model_path)

    ids = model.group_ids()
    assert len(ids) == 2

    names_by_dim = {model.group_dim(gid): model.group_name(gid) for gid in ids}
    assert names_by_dim == {0: "Corner", 2: "Surf"}

    corner_id = next(gid for gid in ids if model.group_name(gid) == "Corner")
    surf_id = next(gid for gid in ids if model.group_name(gid) == "Surf")
    assert model.group_entities(corner_id) == [(0, 1)]
    assert model.group_entities(surf_id) == [(2, 1)]


def test_unknown_group_id_raises(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    unknown = max(model.group_ids()) + 1
    with pytest.raises(IndexError):
        model.group_name(unknown)
    with pytest.raises(IndexError):
        model.group_dim(unknown)
    assert model.group_entities(unknown) == []


def test_curve_segments_and_vertex_positions(geom_model_path):
    model = gecko.GeomModel(geom_model_path)

    # One position per model vertex; the fixture declares a single one, at the origin.
    assert model.vertex_positions() == [[0.0, 0.0, 0.0]]
    assert len(model.vertex_positions()) == model.nb_vertices()

    # It declares no line elements, so it has no curves and therefore no segments.
    assert model.nb_curves() == 0
    assert model.curve_segments() == []


def test_curve_segments_index_the_backing_mesh(tmp_path):
    # A square whose 4 sides are tagged as curves, so there is something to segment.
    msh = """\
$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
2
1 1 "Sides"
2 2 "Surf"
$EndPhysicalNames
$Nodes
4
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
$EndNodes
$Elements
6
1 1 2 1 1 1 2
2 1 2 1 2 2 3
3 1 2 1 3 3 4
4 1 2 1 4 4 1
5 2 2 2 1 1 2 3
6 2 2 2 1 1 3 4
$EndElements
"""
    path = tmp_path / "square.msh"
    path.write_text(msh)
    model = gecko.GeomModel(str(path))

    assert model.nb_curves() == 4
    segments = model.curve_segments()
    assert len(segments) == 4  # one straight segment per side

    vertices = model.mesh_vertices()
    for a, b in segments:
        assert 0 <= a < len(vertices)
        assert 0 <= b < len(vertices)
        assert a != b


def test_volume_boundary_triangles(tmp_path):
    # 2 unit tets glued along one shared face, so the mesh has exactly 1 interior face (the shared
    # one) and 6 boundary faces (4 per tet, minus the 2 that make up the shared one).
    msh = """\
$MeshFormat
2.2 0 8
$EndMeshFormat
$Nodes
5
1 0 0 0
2 1 0 0
3 0 1 0
4 0 0 1
5 1 1 1
$EndNodes
$Elements
2
1 4 0 1 2 3 4
2 4 0 1 2 3 5
$EndElements
"""
    path = tmp_path / "two_tets.msh"
    path.write_text(msh)
    model = gecko.GeomModel(str(path))

    boundary = model.volume_boundary_triangles()
    assert len(boundary) == 6

    # No duplicate faces, and no index out of range.
    vertices = model.mesh_vertices()
    seen = set()
    for tri in boundary:
        assert all(0 <= i < len(vertices) for i in tri)
        key = tuple(sorted(tri))
        assert key not in seen
        seen.add(key)

    # The face shared by both tets (nodes 1,2,3, i.e. indices 0,1,2) must not appear.
    assert (0, 1, 2) not in seen


def test_volume_boundary_triangles_empty_without_tets(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    assert model.volume_boundary_triangles() == []


def test_missing_file_raises(tmp_path):
    with pytest.raises(RuntimeError):
        gecko.GeomModel(str(tmp_path / "does_not_exist.msh"))
