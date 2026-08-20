"""Shared fixtures for the gecko Python test suite — see docs/user-guide/python.md."""

import pytest

# A minimal Gmsh MSH 2.2 ASCII fixture: one vertex (physical group "Corner", dim 0) and one
# triangle (physical group "Surf", dim 2), mirroring the small hand-built fixtures used by the C++
# block/geom test suites (e.g. blocking_creation_tests.cpp's make_minimal_geom_model()).
_MIN_GEOM_MSH = """\
$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
2
0 1 "Corner"
2 2 "Surf"
$EndPhysicalNames
$Nodes
3
1 0 0 0
2 1 0 0
3 0 1 0
$EndNodes
$Elements
2
1 15 2 1 1 1
2 2 2 2 1 1 2 3
$EndElements
"""


@pytest.fixture
def geom_model_path(tmp_path):
    """Path to a minimal, on-disk GeomModel fixture: 1 vertex ("Corner", tag 1) + 1 surface
    ("Surf", tag 1)."""
    path = tmp_path / "min_geom.msh"
    path.write_text(_MIN_GEOM_MSH)
    return str(path)
