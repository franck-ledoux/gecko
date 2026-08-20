"""Tests for the gecko.Blocking façade — see docs/user-guide/python.md."""

import pytest

import gecko

_QUAD_A = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0)]
_QUAD_B = [(1.0, 0.0, 0.0), (2.0, 0.0, 0.0), (2.0, 1.0, 0.0), (1.0, 1.0, 0.0)]
_UNIT_HEX = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0), (0.0, 1.0, 1.0),
]


def test_quad_creation_connectivity_and_deletion(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)

    face_a = blocking.create_quad_block(_QUAD_A)
    face_b = blocking.create_quad_block(_QUAD_B)
    assert face_a != face_b

    blocking.build_connectivity()

    # 2 quads sharing 1 edge: 6 distinct corners, 7 distinct edges, 2 faces.
    assert blocking.nb_cells(0) == 6
    assert blocking.nb_cells(1) == 7
    assert blocking.nb_cells(2) == 2
    assert blocking.is_purely_2d()
    assert blocking.is_valid_topology()

    # Nothing in the (tiny, far-away) geom_model_path fixture is close enough to classify onto —
    # this only checks that classify() runs without raising.
    blocking.classify(1e-9)

    assert blocking.can_delete_face(face_a)
    blocking.delete_face(face_a)
    assert blocking.nb_cells(2) == 1

    with pytest.raises(IndexError):
        blocking.delete_face(face_a)
    with pytest.raises(IndexError):
        blocking.can_delete_face(face_a)


def test_create_quad_block_wrong_corner_count_raises(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    with pytest.raises(ValueError):
        blocking.create_quad_block(_QUAD_A[:3])


def test_write_vtk(tmp_path, geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_quad_block(_QUAD_A)

    out_path = tmp_path / "blocking.vtk"
    blocking.write_vtk(1, str(out_path))

    assert out_path.exists()
    assert out_path.stat().st_size > 0


def test_hex_block_with_cubic_degree(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model, degree=3)

    block_id = blocking.create_hex_block(_UNIT_HEX)
    assert isinstance(block_id, int)
    assert blocking.nb_cells(3) == 1
    assert not blocking.is_purely_2d()


def test_invalid_degree_raises(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    with pytest.raises(ValueError):
        gecko.Blocking(model, degree=2)


def test_nb_cells_invalid_dim_raises(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    with pytest.raises(ValueError):
        blocking.nb_cells(4)
