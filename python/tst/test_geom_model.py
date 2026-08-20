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


def test_missing_file_raises(tmp_path):
    with pytest.raises(RuntimeError):
        gecko.GeomModel(str(tmp_path / "does_not_exist.msh"))
