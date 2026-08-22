# biy — Block It Yourself

`biy` is Gecko's interactive 3D viewer for building a block structure against a geometric model.
It loads a Gmsh `.msh` file as a [`GeomModel`](python.md#geommodel), shows it in a
[Polyscope](https://polyscope.run) window, and lets you build, deform and classify a
[`Blocking`](python.md#blocking) on top of it — from the on-screen panel, by dragging block corners
with the mouse, or by typing Python at the console, all acting on the same live objects.

![The biy viewer: a hex block fitted inside a translucent cylinder, its corners colored by what
they are classified on](../assets/biy.png)

!!! warning "Not covered by CI"
    Unlike the rest of Gecko, `biy` is **not built or tested in CI** — the hosted runners provide no
    GL/windowing stack for Polyscope. The library code underneath it (`Blocking::move_node`, the
    Python facades it displays through) *is* covered by the regular C++ and Python suites, but the
    viewer itself is only ever verified locally. Build and run it yourself after changing it.

## Building

`biy` is opt-in: on top of a Python development environment it needs a GL/windowing stack, so it
isn't built by default.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DGECKO_BUILD_BIY=ON
cmake --build build --target biy
```

## Running

```bash
./build/biy/biy test_data/two_cubes.msh        # cubic blocks (order 3, the default)
./build/biy/biy test_data/two_cubes.msh 1      # straight blocks
```

The optional second argument is the **block order** the session *starts* at: `1` for straight edges,
higher for Bezier edges of that degree. It is not fixed — the `order` field in the panel changes it
at any time, and the whole structure is refitted at the new order (see [Changing the order](#changing-the-order)).

Two things start together:

- a **3D window** showing the model's facets, with a panel of actions; and
- a **Python console** on the terminal you launched from, where `model` and `blocking` are already
  bound.

They are not two copies of the data kept in sync — the console's `blocking` *is* the object the
window draws, so a change made either way shows up immediately in the other.

## The panels

biy replaces Polyscope's own windows with two of its own, on opposite sides of the 3D view and
both the width Polyscope gives its panels (`panel_width` in the config), plus Polyscope's own
"Selection" box below Scene once something has been picked in the 3D view.

**BIY operations** (left) — Polyscope's *Reset view*, *Screenshot*, *View* and *Appearance*
controls, rebuilt here from its public builders, followed by biy's own:

| Button              | What it does                                                                                                 |
| ------------------- | ------------------------------------------------------------------------------------------------------------ |
| Create bounding box | Creates a hex block spanning the model's bounding box, with a 10% margin                                     |
| Build connectivity  | Sews coincident blocks together (see `Blocking::build_connectivity`)                                         |
| Classify            | Snaps the whole blocking onto the model, within the 3 tolerances below                                       |
| Export VTK          | Writes the generated mesh to `biy_blocking.vtk`, with `classification_dim`/`classification_tag` node scalars |

Polyscope's *Debug* section is dropped — it is an internal texture inspector, and unlike the others
it is not exposed as a callable builder.

**Scene** (right) — replaces Polyscope's "Structures" list with a tree that speaks biy's own
vocabulary.
Each row delegates to the structure's own Polyscope panel, so the visibility controls are the usual
ones; a row greys out when the model has nothing of that kind (a boundary-representation file has no
volume, a 2D blocking no blocks).

| Section  | Rows                              |
| -------- | --------------------------------- |
| Geometry | volumes, surfaces, curves, points |
| Blocking | blocks, faces, edges, vertices    |
| Mesh     | quads, hexes                      |

A **volume** is drawn as its true outer boundary, computed directly from the tetrahedral mesh (a
face belonging to exactly one tet); *not* what `polyscope::registerTetMesh` itself would show —
despite its name, a Polyscope volume mesh renders every tet face, interior ones included, with no
way to turn that off. *surfaces* is a separate thing entirely: every surface the file tags, which
on a model with an internal partition between two volumes includes that partition (`two_cubes.msh`:
286 tagged triangles against a 260-triangle outer boundary — the 26-triangle difference *is* the
partition). The two are shown independently for exactly that reason.

**Blocking** shows the block structure itself, **Mesh** the mesh it generates. `subdivisions`
therefore controls the Mesh section alone; the blocking is drawn at its own display resolution,
fine enough to show curvature whatever the mesh is set to. It is capped at **20**, because a 3D
blocking's mesh grows cubically with it — 20 is already 8000 hexes per block, and the window bogs
down before anything else does.

## Orientation gizmo

A small compass sits bottom-right of the 3D view: 6 dots, one per world axis, projected through the
current camera onto a flat 2D canvas — red for X, green for Y, blue for Z, a filled dot and letter
for the positive direction, a hollow ring for the negative one. It updates every frame, so it
reflects the scene's actual orientation continuously, including while the view is being dragged
with the mouse — not just after clicking one of its own dots.

Click any dot to fly the camera to look straight down that axis, centered on the scene, at the same
framing distance `Reset view` uses. It is a flat compass rather than a real 3D widget on purpose:
"bottom-right of the screen" is a screen-space request, and pinning actual 3D geometry to a screen
corner under every zoom level, without it ever clipping through nearby scene content, needs far more
care than projecting 6 directions through the camera's own basis each frame.

## Changing the order

The **order** field in the panel sets the degree of every edge curve, face surface and block volume.
Changing it rebuilds the whole structure at the new degree and refits it onto the model — topology
and classification are carried across untouched, only the representation changes.

That distinction matters: raising the order does not merely add control points, it *uses* them. An
edge lying on a curved model curve can only be its chord at order 1; at order 3 it comes back
following the curve. On a block fitted to `cylinder.msh`, the worst departure from straightness goes
from 0 at order 1 to 0.22 at order 3 — which is the geometry it was always supposed to have.

Lowering it works the same way and is equally lossless to the topology; the geometry simply has
fewer degrees of freedom to follow the model with.

The rebuild is a reclassification rather than a degree elevation, deliberately. Elevation is exact
and would keep each edge's current shape, but that shape is the best fit at the *old* degree, which
is precisely what asking for a higher order is asking to improve. An edge with nothing to fit onto
still lands where elevation would have put it — a straight edge rebuilt straight at any degree is
the same curve.

The panel caps the order at 10. Nothing in the kernel does: the degree is a plain number carried by
each curve. But a fit has only so much to say, and past that the extra control points crowd together
without following the model any closer, while costing on every face and block that touches the edge.

## Classification colors

Corners, block edges **and** block faces are all colored by what `classify()` put them on, so the
state of a blocking being fitted to its model is readable at a glance:

| Color  | Classified on                                                        |
| ------ | -------------------------------------------------------------------- |
| Violet | nothing yet — `classify()` hasn't run, or found nothing in tolerance |
| Yellow | a model vertex                                                       |
| Red    | a model curve                                                        |
| Blue   | a model surface                                                      |
| Green  | a model volume                                                       |

The corner being dragged is drawn white and larger, so it can't be confused with any of these.

Edges and faces are not classified by proximity like corners are, but from **their own boundary**:
an edge goes on the lowest-dimensional entity containing both its corners' classifications, a face
on the one containing all 4 of its edges'. So an edge between 2 corners of one model curve lands on
that curve, while a diagonal whose corners belong to different curves lands on the surface they
share — even though its midpoint may sit closer to a curve. Proximity is only the fallback, used
when the boundary can't decide (an unclassified corner, say).

Classifying a cell also **fits its geometry to what it landed on**, and that reaches inside the cell,
not just around its rim. An edge on a model curve is fitted through points taken on that curve; a
face on a model surface then has its interior pulled onto that surface, and the block behind it
follows its own faces. Without that last step a face would be nothing but a blend of its 4 boundary
edges — which on a sphere leaves the edges on the model and the middle of every face a fifth of the
radius inside it, no matter how high the order goes. Blue faces really are on the blue surface.

## Snapping

`classify()` takes **three** tolerances — vertex, curve, surface — editable in the panel. They are
separate because the scales genuinely differ: one loose enough to catch a surface would snap corners
onto the wrong vertex.

Those same tolerances drive the **snap on release**: let go of a dragged corner and it settles onto
whatever it landed near, with every edge and face touching it reclassified and refitted to match.
The status line reports what it snapped onto. That update is local rather than a full `classify()`
pass — sound precisely because edges and faces infer from their boundary, so no cell the corner
doesn't touch can have changed.

## Control points

**Control points** has one checkbox per kind, drawing the handles that actually drive each curved
cell together with the scaffold joining them:

| Checkbox | Shows                                                | Points per cell |
| -------- | ---------------------------------------------------- | --------------- |
| edges    | each edge's control polygon                          | `order + 1`     |
| faces    | each face's control net, along `u` and `v`           | `(order + 1)²`  |
| blocks   | each block's control lattice, along `u`, `v` and `w` | `(order + 1)³`  |

They are colored differently because they overlap: a face's net already contains its edges' control
points, and a block's lattice contains both. All three are disabled at order 1, where every control
point is simply a block corner.

Expect the lattice to stick out beyond the block it drives — a Bezier lies inside the convex hull of
its control points, not through them, so the handles sit outside the shape they bend.

That same property is why classification **fits** the edge rather than projecting its control points
onto the geometry. Moving a control point onto a geometric curve does not put the *curve* there: it
passes through its 2 endpoints only, staying strictly inside its handles, and so ends up bowed short
of the geometry it should follow. `classify()` instead samples points on the geometry and solves for
the control points that reproduce them — which is exactly why the handles you see do not lie on the
model.

An edge classified on a **curve** additionally leaves each of its ends along that curve's own
tangent there, since the first handle's direction *is* the curve's starting direction. Fitting
positions alone leaves those directions free, and they come out badly wrong — about 30° off on a
plain circular arc, which shows up as handles splaying away from the geometry instead of running
alongside it. A surface has no single such direction, so edges classified on one are fitted by
position only.

!!! note
    Green (volume) is not one of the four states you might expect, but `classify()` does produce it,
    and it shows up more than you'd think on tetrahedral models: `FacetedVolume::distance()` is
    currently a documented stub that returns `0.0` for every point (see
    `geom/inc/gecko/geom/FacetedEntities.h`). So on a model with volumes, any corner that doesn't
    match a vertex, curve or surface classifies onto a volume rather than staying unclassified —
    violet then only appears before `classify()` has run.

## Moving corners

The left mouse button does one of four things, picked with the **mouse mode** radio buttons at the
top of the panel or with a keypress:

| Mode   | Key | Left button                                                    |
| ------ | --- | -------------------------------------------------------------- |
| Camera | `C` | Polyscope's usual navigation: rotate, pan, zoom                |
| Edit   | `E` | Picks up a block corner and moves it                           |
| Cut    | `X` | Cuts the sheet under the cursor — click or `Space` (see below) |
| Delete | `D` | Deletes the block under the cursor — click or `Space`          |

These are genuine modes rather than a modifier like `Ctrl`+drag, because of how Polyscope is
built: it processes camera navigation at the top of each frame, *before* the per-frame user
callback runs. A drag therefore can't be intercepted after the fact — navigation has to be switched
off (`options::doDefaultMouseInteraction`) ahead of the frame the drag happens in. `Ctrl` is also
already taken: Polyscope uses `Shift`+`Ctrl`+drag for zooming.

In **Edit** mode, press the left button on a corner and it follows the mouse in the plane facing the
camera. While held, the corner is drawn larger and in the highlight color, returning to normal on
release. Every edge, face and block touching the corner is refitted live, so the block visibly
deforms as you drag. Releasing snaps it onto the model and reclassifies everything it touches — see
[Snapping](#snapping) above.

The scene itself stays put while you edit. Polyscope normally recomputes the scene's bounding box
and length scale whenever a structure changes, which drags the ground plane along with it — so biy
freezes both to the model once it's loaded (`options::automaticallyComputeSceneExtents`). Growing a
block, or pulling a corner far out, no longer shifts the ground or rescales the view.

## Cutting blocks

**Cut** mode splits blocks along a *sheet*: pick an edge and a point on it, and the cut runs through
every edge parallel to it, right across the blocking.

Point at a block edge and the whole sheet lights up, with a marker on each of its edges showing
exactly where the cut would land. The parameter follows the cursor along the edge, snapping to the
middle as you near it (cutting in half being the common case), and the panel's `Cut at` slider takes
an exact value when pointing is not precise enough. Click **or press `Space`** to cut.

The target is the whole highlighted sheet, not the hairline underneath it: pointing at an edge, at
the highlight drawn over it, or at any of the round markers all aim the same cut. That matters more
than it sounds, because the highlight is drawn thicker than the edges and each marker thicker still,
and a marker appears directly under the cursor that summoned it — so a moment after a sheet lights
up you are no longer pointing at the edge at all.

`Space` cuts as well as a click, because a trackpad tap is an unreliable way to say "here": it fires
a click and a small cursor jolt at the same moment. For the same reason the cut acts on the sheet
**currently highlighted** rather than re-testing what sits under the cursor at the instant of the
click — what you see highlighted is what gets cut — and a key press avoids nudging the pointer at
all. Every cut, and every refusal, also prints a line to the terminal the Python console runs in, so
a press that seems to do nothing says why there and not only in the panel's status line.

Seeing the sheet before committing is the point of the preview: a cut is never local. Two blocks
sewn together share the very same edges, so cutting one drags its neighbour in — otherwise the face
they share would come apart and the blocking would stop being conformal. On `two_cubes.msh`, cutting
across the shared face touches one block; cutting along it touches both. The status line reports how
many edges the cut split and how the block count moved.

A sheet that closes back onto itself has no single well-defined cut — some face would have to be cut
twice, or two blocks would disagree on which side the parameter runs from. biy says so and refuses,
rather than cutting somewhere arbitrary.

### Curvature is kept exactly

The blocks a cut produces are not refitted to their new boundaries; they are the *restriction* of
what they were cut out of. Every edge curve, face surface and block volume is subdivided by De
Casteljau, which is exact — so the cut re-parameterizes the blocking without moving it anywhere.
Meshing after a cut traces the very same points as meshing before: on a block fitted to
`cylinder.msh`, cutting at 1/3 leaves every generated node within 5e-16 of the uncut geometry.

Rebuilding the halves from their boundaries instead would move them, which is why the cut does not
do it: the restriction of a Coons patch to a sub-rectangle is not the Coons patch of that
sub-rectangle's own boundary curves. The difference is invisible on a straight-edged block, where
the two agree, and plain to see on a curved one.

One consequence worth knowing: `classify()` rebuilds every face and block from its boundary — that
is its documented job — and so discards the exact geometry a cut was careful to keep. **Cut first,
classify after** is the order that keeps both.

## Deleting blocks

**Delete** mode removes a block and everything that existed only because of it. Point at a block —
it lights up — then click or press `Space`, exactly as in Cut mode and for the same reasons.

What goes is what the block owned alone. A face, edge or corner it *shared* with a neighbouring
block stays, as that neighbour's own boundary: two blocks sewn along a face become one block with
that face on its boundary, which is what deleting one of them ought to mean. Delete the last block
and the blocking is simply empty, which is a state to be in rather than a broken one — it still
meshes, and still takes a new block.

There is no precondition, unlike deleting a *face* in a 2D blocking: a block is always removable,
and leaving the blocking in several disconnected pieces is a legitimate thing to be in the middle
of. There is also no undo — the status line and the terminal both report what was removed.

The block under the cursor is outlined by a shell drawn just outside it — a hair larger than the
block, not the same size, since a highlight sitting exactly on the faces it marks z-fights with them
and reads on screen as no highlight at all. The status line names the block as the cursor moves, and
says so when the cursor is over none, which a highlight alone cannot.

Aiming needs the **blocks** row ticked in the Scene panel; with it hidden there is nothing to point
at, and the terminal says so rather than the click doing nothing. Only blocks you can see are
reachable — the blocks are drawn as solids, so the pick lands on the outermost one along the line of
sight. Peeling a grid from the outside in works; reaching straight into the middle of one does not.

## Configuration

At startup biy reads `biy_config.json` from the current directory. The file is optional, and may
set only the keys it cares about; anything missing keeps its default. A malformed file is reported
and ignored rather than being fatal. [`biy/biy_config.json`](https://github.com/franck-ledoux/gecko/blob/main/biy/biy_config.json)
is a copy of the defaults:

```json
{
  "panel_width": 305,

  "geometry_volume_color": [0.55, 0.6, 0.7],
  "geometry_surface_color": [0.35, 0.5, 0.8],
  "geometry_curve_color": [0.15, 0.15, 0.2],
  "geometry_curve_radius": 0.004,
  "geometry_point_color": [0.1, 0.1, 0.15],
  "geometry_point_radius": 0.008,
  "model_transparency": 0.45,
  "show_geometry_volumes": true,
  "show_geometry_surfaces": true,
  "show_geometry_curves": true,
  "show_geometry_points": true,

  "corner_radius": 0.01,
  "corner_highlight_radius": 0.02,
  "corner_highlight_color": [1.0, 1.0, 1.0],
  "corner_color_unclassified": [0.6, 0.2, 0.85],
  "corner_color_on_vertex": [1.0, 0.9, 0.1],
  "corner_color_on_curve": [0.9, 0.15, 0.15],
  "corner_color_on_surface": [0.15, 0.4, 0.95],
  "corner_color_on_volume": [0.2, 0.75, 0.3],

  "show_block_edges": true,
  "block_edge_radius": 0.003,
  "block_edge_color": [0.15, 0.15, 0.15],

  "show_edge_control": false,
  "show_face_control": false,
  "show_block_control": false,
  "control_point_radius": 0.006,
  "control_polygon_radius": 0.002,
  "edge_control_color": [0.1, 0.8, 0.8],
  "face_control_color": [0.95, 0.55, 0.1],
  "block_control_color": [0.55, 0.85, 0.3],

  "mesh_color": [0.6, 0.6, 0.65],
  "show_mesh": false,

  "tol_vertex": 0.1,
  "tol_curve": 0.1,
  "tol_surface": 0.1
}
```

| Key                                                                                               | Meaning                                                                                  |
| ------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------- |
| `panel_width`                                                                                     | Width of biy's 2 side panels, before UI scaling                                          |
| `sheet_color`, `sheet_radius`                                                                     | Color and thickness of the highlighted sheet in Cut mode                                 |
| `cut_point_color`, `cut_point_radius`                                                             | Color and size of the markers showing where a cut would land                             |
| `cut_snap_tolerance`                                                                              | How close to an edge's middle the cursor snaps the cut to exactly 0.5                    |
| `delete_highlight_color`                                                                          | Color of the block about to be deleted, in Delete mode                                   |
| `delete_highlight_scale`                                                                          | How much larger than the block its highlight shell is drawn                              |
| `gizmo_size`, `gizmo_dot_radius`                                                                  | Size of the orientation gizmo's canvas, and of a dot at rest                             |
| `gizmo_color_x`, `gizmo_color_y`, `gizmo_color_z`                                                 | Color of each axis's dots and label                                                      |
| `geometry_volume_color`                                                                           | Color of the model's volumes                                                             |
| `geometry_surface_color`                                                                          | Color of the model's surfaces                                                            |
| `geometry_curve_color`, `geometry_curve_radius`                                                   | Color and thickness of the model's curves                                                |
| `geometry_point_color`, `geometry_point_radius`                                                   | Color and sphere radius of the model's points                                            |
| `show_geometry_volumes`, `show_geometry_surfaces`, `show_geometry_curves`, `show_geometry_points` | Which Geometry rows start visible                                                        |
| `mesh_color`, `show_mesh`                                                                         | Color of the generated mesh, and whether it starts visible                               |
| `corner_radius`                                                                                   | Size of a block corner at rest                                                           |
| `corner_highlight_radius`                                                                         | Size of the corner being dragged                                                         |
| `corner_highlight_color`                                                                          | Color of the corner being dragged, RGB in `[0,1]`                                        |
| `corner_color_unclassified`                                                                       | Color of a corner not classified onto anything                                           |
| `corner_color_on_vertex`                                                                          | Color of a corner classified on a model vertex                                           |
| `corner_color_on_curve`                                                                           | Color of a corner classified on a model curve                                            |
| `corner_color_on_surface`                                                                         | Color of a corner classified on a model surface                                          |
| `corner_color_on_volume`                                                                          | Color of a corner classified on a model volume                                           |
| `model_transparency`                                                                              | Opacity of the model surface in `[0,1]`; below 1 so the blocking inside it stays visible |
| `tol_vertex`, `tol_curve`, `tol_surface`                                                          | Starting snapping tolerances, one per entity dimension                                   |
| `show_block_edges`                                                                                | Whether block edges are drawn at startup                                                 |
| `show_control_points`                                                                             | Whether curved edges' control points are drawn at startup                                |
| `control_point_radius`                                                                            | Size of a control point                                                                  |
| `control_polygon_radius`                                                                          | Thickness of the polygon joining them                                                    |
| `control_point_color`                                                                             | Color of both, RGB in `[0,1]`                                                            |
| `block_edge_radius`                                                                               | Thickness of the block edges                                                             |
| `block_edge_color`                                                                                | Color of the block edges, RGB in `[0,1]`                                                 |

Radii are Polyscope *relative* values — a fraction of the scene's bounding box — so they stay
sensible whatever units the model uses.

## The Python console

Everything the panel does is also reachable from the console, which is the point: the panel is a
convenience, the console is the full API.

```python
>>> vs = model.mesh_vertices()
>>> lo = [min(v[k] for v in vs) for k in range(3)]
>>> hi = [max(v[k] for v in vs) for k in range(3)]
>>> blocking.create_hex_block([(lo[0],lo[1],lo[2]), (hi[0],lo[1],lo[2]), (hi[0],hi[1],lo[2]), (lo[0],hi[1],lo[2]),
...                            (lo[0],lo[1],hi[2]), (hi[0],lo[1],hi[2]), (hi[0],hi[1],hi[2]), (lo[0],hi[1],hi[2])])
0
>>> blocking.node_ids()
[0, 1, 2, 3, 4, 5, 6, 7]
>>> blocking.move_node(6, hi[0] + 1.5, hi[1] + 1.5, hi[2] + 1.0)
>>> blocking.classify(0.3)
>>> blocking.node_classification_dims()   # -1 free, 0 vertex, 1 curve, 2 surface, 3 volume
[0, 1, 2, 3, 0, 0, 0, 0]
>>> blocking.write_vtk(2, "blocking.vtk")
```

The console exposes exactly the API documented in [Python Bindings](python.md) — it registers the
same bindings as the standalone `gecko` module, so `import gecko` works here too.

It behaves like Python's own interactive prompt, `codeop` and all: multi-line blocks and bracket
continuations are accepted, showing `...` until the statement is complete, and a block ends with a
blank line. That also means a script piped in on stdin must leave a blank line after any trailing
block.

Press `Ctrl-D` in the console, or close the window, to quit.

## Implementation notes

Polyscope drives rendering from the main thread through its non-blocking `frameTick()`, rather than
its usual blocking `show()`; the Python interpreter is embedded (not a subprocess) and its REPL runs
on a second thread. A single mutex guards the facade objects, taken by the REPL around each
statement and by the render thread around each frame, so both can touch the same state safely. That
avoids the serialization layer a separate-process viewer would need, at the cost of the two sides
sharing a process.

`BIY_SCREENSHOT=<path>` renders a few frames, saves the window to that file and exits — a way to
check the rendering path without a human at the keyboard.
