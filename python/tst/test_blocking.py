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
    blocking.classify(1e-6)

    out_path = tmp_path / "blocking.vtk"
    blocking.write_vtk(1, str(out_path))

    assert out_path.exists()
    content = out_path.read_text()
    assert "POINT_DATA 4" in content
    assert "SCALARS classification_dim int 1" in content
    assert "SCALARS classification_tag int 1" in content

    # Same fixture/tolerance as test_node_classification_dims: 1 corner on the model's one vertex
    # (dim 0), 2 on its one surface (dim 2), 1 outside it and unclassified (dim -1) — order isn't
    # guaranteed to match write_vtk's own node order, so compare as a multiset.
    dim_section = content.split("SCALARS classification_dim int 1\nLOOKUP_TABLE default\n")[1]
    dims = [int(line) for line in dim_section.splitlines()[:4]]
    assert sorted(dims) == [-1, 0, 2, 2]


def test_hex_block_with_cubic_degree(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model, degree=3)

    block_id = blocking.create_hex_block(_UNIT_HEX)
    assert isinstance(block_id, int)
    assert blocking.nb_cells(3) == 1
    assert not blocking.is_purely_2d()


def test_node_classification_dims(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_quad_block(_QUAD_A)

    # Nothing is classified until classify() runs.
    assert blocking.node_classification_dims() == [-1, -1, -1, -1]
    assert len(blocking.node_classification_dims()) == len(blocking.node_ids())

    # The fixture declares exactly one model vertex (the origin) and one surface (the triangle
    # (0,0,0)-(1,0,0)-(0,1,0)). So of the unit square's corners: (0,0,0) hits the vertex, (1,0,0)
    # and (0,1,0) lie on the triangle itself and classify onto the surface, and (1,1,0) is outside
    # it and stays unclassified.
    blocking.classify(1e-6)
    assert blocking.node_classification_dims() == [0, 2, -1, 2]


def test_edge_polylines(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_quad_block(_QUAD_A)

    # A lone quad has 4 edges: samples+1 points and samples segments each.
    for samples in (1, 5):
        points = blocking.edge_vertices(samples)
        segments = blocking.edge_segments(samples)
        assert len(points) == 4 * (samples + 1)
        assert len(segments) == 4 * samples
        # Every segment must index a point that exists, and join consecutive samples.
        for a, b in segments:
            assert 0 <= a < len(points)
            assert b == a + 1


def test_edge_vertices_invalid_samples_raises(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_quad_block(_QUAD_A)
    with pytest.raises(ValueError):
        blocking.edge_vertices(0)


@pytest.mark.parametrize("degree", [1, 2, 3, 4])
def test_every_supported_degree(geom_model_path, degree):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model, degree=degree)
    assert blocking.degree() == degree

    blocking.create_quad_block(_QUAD_A)
    # A quad has 4 edges, each driven by degree+1 control points joined by degree segments.
    assert len(blocking.edge_control_points()) == 4 * (degree + 1)
    assert len(blocking.edge_control_polygons()) == 4 * degree


@pytest.mark.parametrize("degree", [1, 2, 3, 4])
def test_face_and_block_control_nets(geom_model_path, degree):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model, degree=degree)
    blocking.create_hex_block(_UNIT_HEX)
    n = degree + 1

    face_points = blocking.face_control_points()
    face_net = blocking.face_control_nets()
    block_points = blocking.block_control_points()
    block_net = blocking.block_control_lattices()

    # 6 faces with an (n x n) grid each; 1 block with an (n x n x n) grid.
    assert len(face_points) == 6 * n * n
    assert len(block_points) == n**3
    # Per face: n rows and n columns, each of n-1 segments. Per block: 3 axes of n*n lines.
    assert len(face_net) == 6 * 2 * n * (n - 1)
    assert len(block_net) == 3 * n * n * (n - 1)

    # No segment may join two different faces (or two different blocks).
    for a, b in face_net:
        assert a // (n * n) == b // (n * n)
    for a, b in block_net:
        assert max(a, b) < len(block_points)


def test_degree_one_block_control_points_are_its_corners(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model, degree=1)
    blocking.create_hex_block(_UNIT_HEX)

    assert sorted(tuple(p) for p in blocking.block_control_points()) == sorted(set(_UNIT_HEX))
    # The 12 lattice segments of a straight block are exactly the cube's 12 edges.
    assert len(blocking.block_control_lattices()) == 12


@pytest.mark.parametrize("degree", [0, -1])
def test_unusable_degree_raises(geom_model_path, degree):
    model = gecko.GeomModel(geom_model_path)
    with pytest.raises(ValueError):
        gecko.Blocking(model, degree=degree)


def test_degree_has_no_upper_bound(geom_model_path):
    # There used to be one, at 4, and it was an artefact of the implementation rather than a
    # modelling limit: the degree was a template parameter, so the façade had to hold one C++ type
    # per degree and could only offer the ones it had instantiated. Carrying the degree in the
    # object removes the ceiling along with the variant.
    model = gecko.GeomModel(geom_model_path)
    for degree in (5, 7):
        blocking = gecko.Blocking(model, degree=degree)
        assert blocking.degree() == degree
        blocking.create_quad_block(_QUAD_A)
        assert len(blocking.edge_control_points()) == 4 * (degree + 1)


def test_edge_and_face_classification_dims(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_quad_block(_QUAD_A)

    assert blocking.edge_classification_dims() == [-1] * 4
    assert blocking.face_classification_dims() == [-1]

    blocking.classify(1e-6)
    # The fixture's one triangle spans (0,0,0)-(1,0,0)-(0,1,0), so 2 of the quad's edges lie on the
    # surface and the face itself follows them onto it.
    assert len(blocking.edge_classification_dims()) == 4
    assert blocking.face_classification_dims() == [2]


def test_face_grid_covers_3d_blocks_that_mesh_quads_does_not(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_hex_block(_UNIT_HEX)

    # to_mesh() only emits quads for standalone 2D blocks, so a hex's 6 bounding faces would be
    # undrawable without the dedicated face-grid path.
    assert blocking.mesh_quads(1) == []
    assert len(blocking.face_grid_quads(1)) == 6
    assert len(blocking.face_classification_dims()) == 6

    for subdivisions in (1, 3):
        vertices = blocking.face_grid_vertices(subdivisions)
        quads = blocking.face_grid_quads(subdivisions)
        owners = blocking.face_grid_owners(subdivisions)
        assert len(vertices) == 6 * (subdivisions + 1) ** 2
        assert len(quads) == 6 * subdivisions**2
        assert len(owners) == len(quads)
        assert max(i for quad in quads for i in quad) < len(vertices)
        assert sorted(set(owners)) == list(range(6))


def test_snap_node(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_quad_block(_QUAD_A)
    blocking.classify(1e-6)

    # Node 0 sits on the fixture's only model vertex, the origin. Nudge it off, then snap it back.
    blocking.move_node(0, 0.02, 0.02, 0.0)
    assert blocking.node_position(0) == [0.02, 0.02, 0.0]

    blocking.snap_node(0, 0.1, 0.1, 0.1)
    assert blocking.node_position(0) == [0.0, 0.0, 0.0]
    assert blocking.node_classification_dims()[0] == 0


def test_snap_node_unknown_id_raises(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_quad_block(_QUAD_A)
    with pytest.raises(IndexError):
        blocking.snap_node(99, 0.1)


def test_nb_cells_invalid_dim_raises(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    with pytest.raises(ValueError):
        blocking.nb_cells(4)


def _unit_hex_at(x0, x1):
    return [
        (x0, 0.0, 0.0), (x1, 0.0, 0.0), (x1, 1.0, 0.0), (x0, 1.0, 0.0),
        (x0, 0.0, 1.0), (x1, 0.0, 1.0), (x1, 1.0, 1.0), (x0, 1.0, 1.0),
    ]


def test_sheet_edges_reports_the_whole_sheet_before_cutting(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_hex_block(_unit_hex_at(0.0, 1.0))
    blocking.create_hex_block(_unit_hex_at(1.0, 2.0))
    blocking.build_connectivity()

    # Every edge belongs to some sheet, and a sheet is always reported including the edge asked for.
    for index in range(blocking.nb_cells(1)):
        sheet = blocking.sheet_edges(index)
        assert index in sheet
        # A hex edge has 3 parallel siblings; crossing into the neighbour block adds more.
        assert len(sheet) in (4, 6)

    # The 2 blocks share 2 y-running edges, so some sheet must span both blocks (6 edges) and some
    # must stay inside one (4) — that difference is the propagation itself.
    sizes = {len(blocking.sheet_edges(i)) for i in range(blocking.nb_cells(1))}
    assert sizes == {4, 6}


def test_cut_sheet_splits_the_blocks_and_keeps_the_blocking_valid(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_hex_block(_unit_hex_at(0.0, 1.0))
    assert blocking.nb_cells(3) == 1

    sheet = blocking.sheet_edges(0)
    assert len(sheet) == 4
    assert blocking.cut_sheet(0, 0.5)

    assert blocking.nb_cells(3) == 2
    assert blocking.nb_cells(0) == 12
    assert blocking.is_valid_topology()
    # The cut's new corners must be addressable like any other, or they could not be dragged.
    assert len(blocking.node_ids()) == 12


def test_cut_sheet_refuses_out_of_range_parameters(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_hex_block(_unit_hex_at(0.0, 1.0))

    for bad in (0.0, 1.0, -0.5, 2.0):
        assert not blocking.cut_sheet(0, bad)
    assert blocking.nb_cells(3) == 1


def test_cut_sheet_rejects_an_unknown_edge_index(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_hex_block(_unit_hex_at(0.0, 1.0))

    with pytest.raises(IndexError):
        blocking.cut_sheet(999, 0.5)
    with pytest.raises(IndexError):
        blocking.sheet_edges(-1)


# The 1x1x2 block of the cut scenario, in create_hex_block's expected HEX8 order.
_TALL_BOX = [
    (0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (1.0, 1.0, 0.0), (0.0, 1.0, 0.0),
    (0.0, 0.0, 2.0), (1.0, 0.0, 2.0), (1.0, 1.0, 2.0), (0.0, 1.0, 2.0),
]


def _edge_between(blocking, p, q):
    """The edge joining 2 positions, as a position in the order the cut API speaks."""
    ends = blocking.edge_vertices(1)
    for i in range(blocking.nb_cells(1)):
        pair = {tuple(round(x, 9) for x in ends[i * 2]), tuple(round(x, 9) for x in ends[i * 2 + 1])}
        if pair == {p, q}:
            return i
    raise AssertionError("no edge between %r and %r" % (p, q))


def test_cutting_a_tall_block_in_half_gives_two_blocks_of_volume_one(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_hex_block(_TALL_BOX)
    assert blocking.block_volumes() == pytest.approx([2.0])

    edge = _edge_between(blocking, (0.0, 0.0, 0.0), (0.0, 0.0, 2.0))
    assert blocking.cut_sheet(edge, 0.5)

    assert blocking.nb_cells(3) == 2
    assert blocking.is_valid_topology()

    # Measured from each block's own stored geometry, coarsely and finely: a half carrying the wrong
    # sub-volume, or the right one in a rotated frame, still matches its corners and misses in
    # between. Volume is the check that notices.
    for subdivisions in (1, 2, 4):
        assert blocking.block_volumes(subdivisions) == pytest.approx([1.0, 1.0])

    # The 2 halves meet along one face carrying exactly the 4 midplane corners.
    midplane = sorted(tuple(round(c, 9) for c in blocking.node_position(i))
                      for i in blocking.node_ids()
                      if abs(blocking.node_position(i)[2] - 1.0) < 1e-9)
    assert midplane == [(0.0, 0.0, 1.0), (0.0, 1.0, 1.0), (1.0, 0.0, 1.0), (1.0, 1.0, 1.0)]


def test_cutting_off_centre_splits_the_volume_in_the_same_ratio(geom_model_path):
    # At exactly 0.5 a swapped side, a transposed frame and a mirrored half all still measure 1 and
    # 1; off-centre they cannot, which makes this the sharper check of the two.
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_hex_block(_TALL_BOX)

    edge = _edge_between(blocking, (0.0, 0.0, 0.0), (0.0, 0.0, 2.0))
    assert blocking.cut_sheet(edge, 0.25)

    assert sorted(blocking.block_volumes(4)) == pytest.approx([0.5, 1.5])
    # A quarter of the way from the edge's own start, so the cut plane lands at z=0.5.
    zs = sorted({round(blocking.node_position(i)[2], 9) for i in blocking.node_ids()})
    assert zs == [0.0, 0.5, 2.0]


@pytest.mark.parametrize("degree", [1, 3])
def test_every_edge_cuts_a_block_into_two_positive_volumes(geom_model_path, degree):
    # Which frame each half lands in is CGAL's choice, and frame bugs show on some edges while
    # sparing others — so sweep all 12 rather than trusting one.
    for index in range(12):
        model = gecko.GeomModel(geom_model_path)
        blocking = gecko.Blocking(model, degree=degree)
        blocking.create_hex_block(_TALL_BOX)
        assert blocking.cut_sheet(index, 0.25)

        volumes = blocking.block_volumes(4)
        assert len(volumes) == 2
        # Positive, so neither half came out inverted, and adding back up to the whole.
        assert min(volumes) > 0.0
        assert sum(volumes) == pytest.approx(2.0)
        assert sorted(volumes) == pytest.approx([0.5, 1.5])


@pytest.mark.parametrize("degree", [1, 3])
def test_a_sheet_stays_in_one_plane_after_an_earlier_cut(geom_model_path, degree):
    # Needs 3 things at once, which is why cutting once, or cutting at 0.5, never showed it: a
    # second cut, an off-centre parameter, and an edge the previous cut created inside a face. Such
    # an edge takes its curve from a surface while its nodes came from evaluating another curve —
    # equal as real numbers, not as doubles — so deciding which end the curve starts at by exact
    # position equality picks the wrong one, cutting that edge at 1-t and tearing the sheet.
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model, degree=degree)
    blocking.create_hex_block(_TALL_BOX)

    edge = _edge_between(blocking, (0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert blocking.cut_sheet(edge, 0.4)

    ends = blocking.edge_vertices(1)
    for i in range(blocking.nb_cells(1)):
        sheet = blocking.sheet_edges(i)
        assert sheet, "edge %d lost its sheet" % i
        a, b = ends[i * 2], ends[i * 2 + 1]
        # The cut runs perpendicular to the edges it splits, so on this axis-aligned block every cut
        # point shares the coordinate of the axis the aimed edge runs along.
        axis = max(range(3), key=lambda k: abs(a[k] - b[k]))
        along = [p[axis] for p in blocking.sheet_cut_points(i, 0.4)]
        assert max(along) - min(along) < 1e-12, (
            "sheet through edge %d is not planar: %r" % (i, sorted(along)))


@pytest.mark.parametrize("degree", [1, 3])
def test_repeated_off_centre_cuts_conserve_the_block(geom_model_path, degree):
    # The same failure from where a user sees it: an out-of-plane sheet still leaves a valid map, so
    # only the geometry notices — the mesh drifts off the block it came from.
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model, degree=degree)
    blocking.create_hex_block(_TALL_BOX)

    for a, b in (((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
                 ((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
                 ((0.0, 0.0, 0.0), (0.0, 0.0, 2.0))):
        assert blocking.cut_sheet(_edge_between(blocking, a, b), 0.4)
        assert blocking.is_valid_topology()

        volumes = blocking.block_volumes(2)
        assert min(volumes) > 0.0
        assert sum(volumes) == pytest.approx(2.0)
        # Nothing the cut generates may leave the box the blocking started as.
        for p in blocking.mesh_vertices(2):
            assert -1e-9 <= p[0] <= 1 + 1e-9
            assert -1e-9 <= p[1] <= 1 + 1e-9
            assert -1e-9 <= p[2] <= 2 + 1e-9

    assert blocking.nb_cells(3) == 8


def test_delete_block_takes_only_what_the_block_owned(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_hex_block(_unit_hex_at(0.0, 1.0))
    blocking.create_hex_block(_unit_hex_at(1.0, 2.0))
    blocking.build_connectivity()
    assert blocking.nb_cells(0) == 12
    assert blocking.nb_cells(3) == 2

    blocking.delete_block(0)

    # The shared face, its edges and its corners survive as the other block's boundary; the 4
    # corners the deleted block alone owned go with it.
    assert blocking.nb_cells(3) == 1
    assert blocking.nb_cells(0) == 8
    assert blocking.nb_cells(2) == 6
    assert blocking.is_valid_topology()
    assert len(blocking.mesh_hexes(1)) == 1

    # Node ids pointing at corners that no longer exist must be gone, or reading them would touch
    # freed attributes.
    assert len(blocking.node_ids()) == blocking.nb_cells(0)
    for node_id in blocking.node_ids():
        blocking.node_position(node_id)


def test_delete_block_after_a_cut(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model, degree=3)
    blocking.create_hex_block(_unit_hex_at(0.0, 1.0))
    assert blocking.cut_sheet(0, 0.4)
    assert blocking.nb_cells(3) == 2

    volumes_before = sorted(blocking.block_volumes(2))
    blocking.delete_block(0)

    assert blocking.nb_cells(3) == 1
    assert blocking.is_valid_topology()
    # Whichever half survived kept its own volume, untouched by the deletion.
    assert blocking.block_volumes(2)[0] == pytest.approx(volumes_before[0]) or \
           blocking.block_volumes(2)[0] == pytest.approx(volumes_before[1])
    # And no stray quads: every remaining face still belongs to the surviving block.
    assert blocking.mesh_quads(1) == []


def test_delete_every_block_then_start_again(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_hex_block(_unit_hex_at(0.0, 1.0))
    blocking.create_hex_block(_unit_hex_at(1.0, 2.0))
    blocking.build_connectivity()

    while blocking.nb_cells(3) > 0:
        blocking.delete_block(0)
        assert blocking.is_valid_topology()

    # An empty blocking is a state to be in, not a broken one.
    assert blocking.nb_cells(0) == 0
    assert blocking.node_ids() == []
    assert blocking.mesh_hexes(1) == []
    blocking.create_hex_block(_unit_hex_at(0.0, 1.0))
    assert blocking.nb_cells(3) == 1


def test_delete_block_rejects_an_unknown_index(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model)
    blocking.create_hex_block(_unit_hex_at(0.0, 1.0))

    with pytest.raises(IndexError):
        blocking.delete_block(7)
    with pytest.raises(IndexError):
        blocking.delete_block(-1)
    assert blocking.nb_cells(3) == 1


def test_deleting_blocks_keeps_every_corner_addressable(geom_model_path):
    # Erasing a block's darts disturbs the vertex orbits of the corners it shared, and CGAL rebuilds
    # an attribute rather than re-pointing it whenever it has to — so a corner that is still there
    # can come back as a different attribute. Pruning the id map without re-indexing then quietly
    # loses it, and node_ids() stops listing corners the blocking still has. Only shows up when the
    # blocks are removed in an order that leaves survivors touching at a corner or an edge alone.
    import random

    random.seed(7)
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model, degree=3)
    blocking.create_hex_block(_TALL_BOX)
    for _ in range(3):
        assert blocking.cut_sheet(random.randrange(blocking.nb_cells(1)), 0.5)

    while blocking.nb_cells(3) > 0:
        blocking.delete_block(random.randrange(blocking.nb_cells(3)))
        assert blocking.is_valid_topology()
        # Every corner the blocking has must still be addressable, and every id must still resolve.
        assert len(blocking.node_ids()) == blocking.nb_cells(0)
        for node_id in blocking.node_ids():
            blocking.node_position(node_id)
        assert len(blocking.node_classification_dims()) == blocking.nb_cells(0)


def test_deleting_blocks_never_bends_a_straight_grid(geom_model_path):
    # A cut of a straight block gives straight halves, and a deletion moves nothing at all — so a
    # grid built this way has no reason to bend, whatever order its blocks come out in. Measured on
    # the control points, which is what the viewer draws.
    import random

    random.seed(3)
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model, degree=3)
    blocking.create_hex_block(_TALL_BOX)
    for _ in range(3):
        assert blocking.cut_sheet(random.randrange(blocking.nb_cells(1)), 0.5)

    def worst_bend():
        points = blocking.edge_control_points()
        worst = 0.0
        for e in range(len(points) // 4):
            a, b, c, d = points[e * 4:e * 4 + 4]
            chord = [d[k] - a[k] for k in range(3)]
            length2 = sum(v * v for v in chord)
            for p in (b, c):
                ap = [p[k] - a[k] for k in range(3)]
                t = sum(ap[k] * chord[k] for k in range(3)) / length2 if length2 > 0 else 0.0
                worst = max(worst, sum((ap[k] - t * chord[k]) ** 2 for k in range(3)) ** 0.5)
        return worst

    assert worst_bend() < 1e-12
    while blocking.nb_cells(3) > 0:
        blocking.delete_block(random.randrange(blocking.nb_cells(3)))
        assert worst_bend() < 1e-12


def test_edge_bends_is_zero_on_a_straight_blocking(geom_model_path):
    # The diagnostic's own contract: a straight blocking reads as zero whatever is done to it, so a
    # non-zero reading is always worth chasing rather than sometimes expected.
    import random

    random.seed(11)
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model, degree=3)
    blocking.create_hex_block(_TALL_BOX)
    assert max(blocking.edge_bends()) < 1e-12

    for _ in range(3):
        assert blocking.cut_sheet(random.randrange(blocking.nb_cells(1)), round(random.uniform(0.1, 0.9), 3))
        assert max(blocking.edge_bends()) < 1e-12
    while blocking.nb_cells(3) > 0:
        blocking.delete_block(random.randrange(blocking.nb_cells(3)))
        assert not blocking.edge_bends() or max(blocking.edge_bends()) < 1e-12

    # And it reports one entry per edge, so a caller can say which edge is at fault.
    blocking.create_hex_block(_TALL_BOX)
    assert len(blocking.edge_bends()) == blocking.nb_cells(1)


def test_set_degree_refits_the_blocking_onto_the_model(geom_model_path):
    # The point of changing order is not to add control points but to use them: at degree 1 an edge
    # lying on a model curve can only be its chord, and raising the order lets it follow the curve.
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model, degree=1)
    blocking.create_quad_block(_QUAD_A)
    blocking.classify(1e-6)

    assert blocking.degree() == 1
    assert len(blocking.edge_control_points()) == 4 * 2
    # A degree-1 edge has no interior control point, so it cannot bend whatever it lies on.
    assert max(blocking.edge_bends()) == pytest.approx(0.0)

    topology = (blocking.nb_cells(0), blocking.nb_cells(1), blocking.nb_cells(2))
    classification = blocking.edge_classification_dims()

    blocking.set_degree(3, 1e-6)

    assert blocking.degree() == 3
    assert len(blocking.edge_control_points()) == 4 * 4
    # Topology and classification are carried across untouched: only the representation changed.
    assert (blocking.nb_cells(0), blocking.nb_cells(1), blocking.nb_cells(2)) == topology
    assert blocking.edge_classification_dims() == classification

    # And back down again, on the same object — which the old degree-in-the-type design could not do
    # at all, since a degree-3 blocking was a different C++ type from a degree-1 one.
    blocking.set_degree(1, 1e-6)
    assert blocking.degree() == 1
    assert len(blocking.edge_control_points()) == 4 * 2
    assert (blocking.nb_cells(0), blocking.nb_cells(1), blocking.nb_cells(2)) == topology


def test_set_degree_rejects_a_degree_below_one(geom_model_path):
    model = gecko.GeomModel(geom_model_path)
    blocking = gecko.Blocking(model, degree=3)
    blocking.create_quad_block(_QUAD_A)
    with pytest.raises(ValueError):
        blocking.set_degree(0, 1e-6)
    assert blocking.degree() == 3
