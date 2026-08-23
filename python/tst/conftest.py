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

# The unit square as a full B-Rep: 4 tagged vertices (elementary 1-4), 4 tagged boundary curves
# (10-13, one segment each) and 1 tagged surface (20, 2 triangles) — the Python mirror of the C++
# suite's make_square_geom_model(). Needed wherever a test has to tell a corner sitting on a model
# *vertex* from one merely on a curve.
_SQUARE_MSH = """\
$MeshFormat
2.2 0 8
$EndMeshFormat
$PhysicalNames
3
0 1 "Vertices"
1 2 "Curves"
2 3 "Surf"
$EndPhysicalNames
$Nodes
4
1 0 0 0
2 1 0 0
3 1 1 0
4 0 1 0
$EndNodes
$Elements
10
1 15 2 1 1 1
2 15 2 1 2 2
3 15 2 1 3 3
4 15 2 1 4 4
5 1 2 2 10 1 2
6 1 2 2 11 2 3
7 1 2 2 12 3 4
8 1 2 2 13 4 1
9 2 2 3 20 1 2 3
10 2 2 3 20 1 3 4
$EndElements
"""


@pytest.fixture
def square_model_path(tmp_path):
    """Path to a unit-square GeomModel fixture: 4 vertices, 4 boundary curves and 1 surface."""
    path = tmp_path / "square_geom.msh"
    path.write_text(_SQUARE_MSH)
    return str(path)
