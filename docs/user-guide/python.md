# Python Bindings

Gecko ships a `gecko` Python extension module built with [pybind11](https://pybind11.readthedocs.io/).
Rather than binding every internal C++ class, it exposes a small number of **façade manager
classes** — `GeomModel` and `Blocking` today — in the spirit of
[gmsh's Python API](https://gmsh.info/doc/texinfo/gmsh.html#Gmsh-API): every façade method takes
and returns only primitive types (`str`, `int`, `float`, and lists/tuples of those). Internal C++
types (`Point3d`, CGAL dart/attribute handles, `UnstructuredMesh`, ...) never cross into Python;
faces and blocks are identified to Python solely by plain `int` ids assigned by the façade.

!!! note
    This page, and the examples on it, are exercised by [`python/tst/`](https://github.com/franck-ledoux/gecko/tree/main/python/tst),
    which CI runs on every push. If this page and the actual behavior ever drift apart, those tests
    (and therefore CI) fail — so what's documented here is always what the code actually does.

## Building

The Python module is opt-in (it needs a Python development environment — headers included, not
just an interpreter):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DGECKO_BUILD_PYTHON=ON
cmake --build build --target gecko
```

This produces a `gecko` extension module (e.g. `gecko.cpython-312-x86_64-linux-gnu.so`) under
`build/python/`.

## Using it

Point `PYTHONPATH` at the directory containing the built module, then import it normally:

```bash
PYTHONPATH=build/python python3 -c "import gecko; print(gecko.hello())"
```

```
Hello from gecko!
```

### GeomModel

`GeomModel` loads a geometric model from a Gmsh MSH file (triangles for a boundary representation,
tetrahedra to also reconstruct volumes) and lets you inspect its entities and physical groups:

```python
import gecko

model = gecko.GeomModel("my_model.msh")
print(model.nb_vertices(), model.nb_curves(), model.nb_surfaces(), model.nb_volumes())
print(model.surface_tags())          # e.g. [1, 2, 3]

for group_id in model.group_ids():
    print(model.group_name(group_id), model.group_dim(group_id))
    print(model.group_entities(group_id))  # [(dim, entity_tag), ...]
```

### Blocking

`Blocking` builds a structured quad/hex blocking against a `GeomModel`. Block corners are passed
as lists of `(x, y, z)` tuples; created faces/blocks are identified by the `int` id `Blocking`
hands back. `degree=1` (the default) gives straight edges; `degree=3` gives cubic-Bezier edges that
`classify()` can bend onto curved geometry. Any degree of 1 or more works, and `set_degree()` changes
it on an existing blocking.

```python
import gecko

model = gecko.GeomModel("my_model.msh")
blocking = gecko.Blocking(model)  # model must stay alive for as long as blocking is used

face_a = blocking.create_quad_block([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)])
face_b = blocking.create_quad_block([(1, 0, 0), (2, 0, 0), (2, 1, 0), (1, 1, 0)])
blocking.build_connectivity()  # sews the shared edge between face_a and face_b

# Snap onto the model. Corners are projected onto it, and every edge, face and block is then
# fitted to what it landed on — including the *inside* of each face, which is pulled onto its model
# surface rather than left as a blend of the 4 boundary curves.
#
# An edge classified on a *surface* is pinned to one definite curve of it rather than to whichever
# one projecting each sample separately happens to trace: the surface's section by the plane through
# the edge's 2 ends containing the surface's own normal. A surface holds infinitely many curves
# between one pair of its points, and without that, 2 edges meeting at a corner need not agree on
# which of them they take.
blocking.classify(tol_vertex=1e-6, tol_curve=1e-3, tol_surface=1e-2)

# The degree is carried by the geometry, not by its C++ type, so it can be raised or lowered on a
# blocking that already exists. Topology and classification are kept; only the representation
# changes — and raising the order does not just add control points, it uses them: an edge that could
# only be a chord at degree 1 comes back following the model curve it lies on.
blocking.set_degree(3, tol_vertex=1e-6, tol_curve=1e-3, tol_surface=1e-2)
print(blocking.degree())

# What each cell ended up on: -1 free, 0 vertex, 1 curve, 2 surface, 3 volume. Edges and faces are
# inferred from their own boundary, not from proximity — see the biy guide.
print(blocking.node_classification_dims(), blocking.edge_classification_dims())

# Snap one corner back onto the model after moving it, reclassifying only what touches it.
blocking.move_node(0, 0.01, 0.01, 0.0)
blocking.snap_node(0, tol_vertex=0.1)

print(blocking.nb_cells(2), blocking.is_valid_topology())

if blocking.can_delete_face(face_a):
    blocking.delete_face(face_a)

# How far each edge departs from a straight line: the largest distance from one of its interior
# control points to its own chord. A straight blocking reads as zero here whatever is done to it, so
# a non-zero entry says the geometry is at fault and a zero one says the drawing is.
print(max(blocking.edge_bends()))

# Cells are named two different ways here, and the difference matters. What is *drawn* comes back
# as flat arrays a renderer indexes by position — edge_vertices(), edge_bends(), block_volumes(),
# mesh_hexes(). What is *acted on* is named by id, because a position is only true until the next
# operation renumbers the traversal, and cutting or deleting is exactly what does that. edge_ids(),
# face_ids() and block_ids() bridge the two: one id per cell, in the order those arrays index.
edge_ids = blocking.edge_ids()
print(len(edge_ids) == len(blocking.edge_bends()))   # True — same order, same length

# Deleting a block takes with it every face, edge and corner that existed only because of it;
# whatever it shared with a neighbour stays, as that neighbour's boundary.
blocking.delete_block(blocking.block_ids()[0])

# Cutting: pick an edge by id, and a parameter along it. The cut runs through every edge parallel to
# that one, across every block they touch — sheet_edges() reports which, and sheet_cut_points()
# where, before anything is modified. sheet_edges() answers in positions, because what it is for is
# highlighting the sheet on screen.
target = blocking.edge_ids()[0]
sheet = blocking.sheet_edges(target)            # [] if the sheet cannot be cut evenly
preview = blocking.sheet_cut_points(target, 0.5)  # one point per edge of `sheet`, in the same order
blocking.cut_sheet(target, 0.5)                 # False, changing nothing, if the cut is impossible

# The inverse: take a whole layer out and glue back what was either side of it. Where 2 corners
# merge, the more constrained classification wins — a model vertex over a curve, a curve over a
# surface — and the merged corner is projected onto it, so a blocking fitted to a model does not
# drift off its features. A sheet holding every block collapses too, leaving the blocking empty, so
# an unclassified grid can be taken apart one sheet at a time. Returns False, changing nothing, when
# one of the sheet's edges joins 2 corners on 2 *different* model vertices — merging them would leave
# one of those vertices with no corner on it — or when the sheet cannot be collapsed at all.
blocking.delete_sheet(target, tol_vertex=1e-6, tol_curve=1e-3, tol_surface=1e-2)

# Pillowing: insert a layer of blocks along a *nappe* — a sheet of block faces cutting the blocking
# in two. It may close back on itself, isolating a set of blocks from everything around them, or run
# clean through the structure and out on its boundary. block_faces() and face_blocks() are what a
# nappe is named with: the faces of a block, and the 1 or 2 blocks either side of a face.
block = blocking.block_ids()[0]
nappe = blocking.block_faces(block)              # its 6 faces: the smallest closed nappe there is
sides = blocking.face_blocks(nappe[0])           # 2 blocks, or 1 where the face is on the boundary

# The second argument names the side that shrinks; the other one does not move at all. It has to be
# named — a nappe through the middle of a structure has 2 block sides and nothing tells them apart —
# and where the nappe lies on the boundary of the blocking, the side that stays *is* the model's own
# boundary. The thickness is a fraction of the mean edge length at each corner that moves, so the
# layer follows the size of the blocks it is inserted between. A corner the nappe cuts through
# becomes 2: the outside one keeps its classification, the inside one keeps only what it is still on
# after moving, so a blocking nobody classified stays unclassified. Returns False, changing nothing,
# when what was given is not a nappe: a face named twice, a standalone quad block, a nappe that does
# not separate the named side from the rest, or one that is not manifold along its own edges.
blocking.pillow(nappe, block, thickness=0.25, tol_vertex=1e-6, tol_curve=1e-3, tol_surface=1e-2)

# Collapsing a chord: a *chord* is the column of blocks strung together through opposite faces — the
# dual curve of the structure, where a sheet is its dual surface. Taking it out means folding it:
# each cross-section folds onto one of its 2 diagonals, the 2 corners off that diagonal meeting in
# the middle, so the 2 blocks that were only edge-neighbours across the fold end up sharing a face.
# Each edge of the hinge itself loses one of the blocks around it: a valence-4 edge comes out with
# valence 3, which is what the operation is for. Folding is the only way out that leaves a
# blocking behind — merging each block's 2 opposite side faces instead would contract edges shared
# with blocks *outside* the column and leave those degenerate.
#
# face_corners() runs the face's 4 corners round its perimeter, so [0]/[2] and [1]/[3] are its 2
# diagonals: the hinge is a corner, and the fold runs along the diagonal through it.
face = blocking.face_ids()[0]
hinge = blocking.face_corners(face)[0]

# Where 2 corners meet, the more constrained classification wins, exactly as for delete_sheet().
# Returns False, changing nothing, when the chord runs back through a block it has already been
# through — a chord closing into a ring or crossing itself has no single fold — when the 2 corners
# meeting are on 2 *different* model vertices, when they are already joined by an edge (folding
# would leave it a loop), or when a block outside the column has both of them as its own corners
# (folding would leave it degenerate).
blocking.collapse_chord(face, hinge, tol_vertex=1e-6, tol_curve=1e-3, tol_surface=1e-2)

# And the inverse: opening a chord puts a column back. The 2 named faces, from edge_faces(), are
# where the fan of blocks around the edge is cut — the fan comes apart in 2 arcs, the edge comes
# apart in 2, and a block is inserted in the gap. Where the edge is on the boundary of the blocking
# its fan is already open at both ends, so one of the 2 is a boundary face and cutting at it costs
# nothing. edge_corners() and face_corners() are how a caller holding positions finds the cells it
# means.
edge = blocking.edge_ids()[0]
fan = blocking.edge_faces(edge)

# How far the column runs is not said: the walk carries the 2 cuts from one edge to the next and
# stops when nothing carries them further. Returns False, changing nothing, when the 2 faces are the
# same or do not both carry the edge, when the chain runs back into itself, when it stops somewhere
# the cuts do not leave the blocks around a corner in exactly 2 groups, or when the walk finds more
# than one way to carry on — the structure offering several columns from that start, which is the
# caller's to choose between rather than the operation's to guess.
blocking.open_chord(edge, fan[0], fan[1], thickness=0.25, tol_vertex=1e-6, tol_curve=1e-3,
                    tol_surface=1e-2)

# The 2 are not inverses cell for cell, and cannot be: a fold takes with it the corners that belonged
# only to the column it closed, so an opening builds a column on the corners it now finds. It puts a
# column back, not that column.

# Undo. Every operation that changes the blocking snapshots it first, so taking an edit back is
# putting that snapshot in place — not a stack of inverse operations, because these operations have
# no inverses (collapsing the layer a cut just made does not restore the block). Ids survive an undo,
# so a block id from before the undone edit finds its block again. An operation the kernel *refuses*
# costs no undo step, and any new edit discards the redo history.
print(blocking.can_undo(), blocking.can_redo())
blocking.undo()
blocking.redo()

# A snapshot is a whole copy of the blocking, so the history costs depth times its size: a few
# thousand cells at degree 3 is a couple of megabytes each.
blocking.set_history_depth(20)

# What each block encloses, from its own stored geometry — exact at 1 subdivision when its faces are
# planar, converging as it grows otherwise. Negative means that block came out inverted.
print(blocking.block_volumes(subdivisions=1))

# Which block each generated cell came from, one entry per cell of mesh_hexes()/mesh_quads(). The
# mesh is emitted block by block, so this is what lets a per-block value — a colour, say — be spread
# onto every cell subdividing that block.
print(blocking.mesh_hex_owners(subdivisions=2))   # into the block_volumes()/delete_block() order
print(blocking.mesh_quad_owners(subdivisions=2))  # standalone quad blocks only; a hex emits no quads

# The halves keep the geometry they were cut out of exactly (De Casteljau subdivision, not a refit),
# so meshing after a cut traces the same points as before. Note classify() rebuilds faces and blocks
# from their boundaries and so discards that — cut first, classify after.

# The exported mesh's nodes carry 2 POINT_DATA scalars, "classification_dim"/"classification_tag"
# (-1/-1 when unclassified) — the same dim/tag pair node_classification_dims() reports, for every
# node the subdivision generates, not just the original block corners.
#
# The file carries more than the blocks. Every block edge lying on a model curve is written as a VTK
# line, every block face lying on a model surface as a VTK quad, so the export describes the boundary
# and not only the volume. CELL_DATA then gives, for every element: the same "classification_dim"/
# "classification_tag" pair (the curve id on a line, the surface id on a quad) and "block_id", which
# names the block a hexahedron came from and reads -1 on a line or a quad.
blocking.write_vtk(subdivisions=4, path="blocking.vtk")
```

## Running the tests

The Python test suite (`python/tst/`, using [pytest](https://pytest.org)) is registered as a
regular CTest test (`test_python`) alongside every C++ suite, once `GECKO_BUILD_TESTS` (the
default) and `GECKO_BUILD_PYTHON` are both on:

```bash
pip install -r python/requirements-test.txt
cmake --build build
cd build && ctest -R test_python --output-on-failure
```

## Code coverage

Because `gecko`'s C++ source (`python/src/bindings.cpp`, over the `Gecko::App` façades) is compiled by the very same CMake build
as the rest of the library, it picks up the project's `--coverage` instrumentation automatically
whenever `GECKO_CODE_COVERAGE`/`GECKO_CODE_COVERAGE_HTML_REPORT` is on (see
[Building & Testing](../developer-guide/testing.md#code-coverage-html-report)) — and because
`test_python` is a regular CTest test, running it exercises that instrumented code exactly like any
C++ test does. In other words: **Python tests contribute to the same C++ coverage report as the
C++ tests**, covering both the binding glue code and whatever C++ library code it calls into — no
separate coverage tooling needed for the Python side.
