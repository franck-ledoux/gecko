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
./build/biy/biy test_data/two_cubes.msh
```

Two things start together:

- a **3D window** showing the model's facets, with a panel of actions; and
- a **Python console** on the terminal you launched from, where `model` and `blocking` are already
  bound.

They are not two copies of the data kept in sync — the console's `blocking` *is* the object the
window draws, so a change made either way shows up immediately in the other.

## The panel

| Button              | What it does                                                                             |
| ------------------- | ---------------------------------------------------------------------------------------- |
| Create bounding box | Creates a hex block spanning the model's bounding box, with a 10% margin                 |
| Build connectivity  | Sews coincident blocks together (see `Blocking::build_connectivity`)                     |
| Classify            | Snaps the blocking onto the model's vertices/curves/surfaces, within the given tolerance |
| Export VTK          | Writes the generated mesh to `biy_blocking.vtk`                                          |

`subdivisions` controls how finely the blocking is meshed for display and export: `1` shows the raw
block structure, higher values show the mesh it generates. **Show block edges** draws the block
structure's own edges as a curve network — the edges of the blocks themselves, traced along their
curves, as opposed to the subdivision lines of the mesh.

## Corner colors

Block corners are colored by what `classify()` put them on, so the state of a blocking being fitted
to its model is readable at a glance:

| Corner | Classified on                                                        |
| ------ | -------------------------------------------------------------------- |
| Violet | nothing yet — `classify()` hasn't run, or found nothing in tolerance |
| Yellow | a model vertex                                                       |
| Red    | a model curve                                                        |
| Blue   | a model surface                                                      |
| Green  | a model volume                                                       |

The corner being dragged is drawn white and larger, so it can't be confused with any of these.

!!! note
    Green (volume) is not one of the four states you might expect, but `classify()` does produce it,
    and it shows up more than you'd think on tetrahedral models: `FacetedVolume::distance()` is
    currently a documented stub that returns `0.0` for every point (see
    `geom/inc/gecko/geom/FacetedEntities.h`). So on a model with volumes, any corner that doesn't
    match a vertex, curve or surface classifies onto a volume rather than staying unclassified —
    violet then only appears before `classify()` has run.

## Moving corners

The left mouse button does one of two things, picked with the **mouse mode** radio buttons at the
top of the panel or with a keypress:

| Mode   | Key | Left button                                     |
| ------ | --- | ----------------------------------------------- |
| Camera | `C` | Polyscope's usual navigation: rotate, pan, zoom |
| Edit   | `E` | Picks up a block corner and moves it            |

The two are a genuine mode rather than a modifier like `Ctrl`+drag, because of how Polyscope is
built: it processes camera navigation at the top of each frame, *before* the per-frame user
callback runs. A drag therefore can't be intercepted after the fact — navigation has to be switched
off (`options::doDefaultMouseInteraction`) ahead of the frame the drag happens in. `Ctrl` is also
already taken: Polyscope uses `Shift`+`Ctrl`+drag for zooming.

In **Edit** mode, press the left button on a corner and it follows the mouse in the plane facing the
camera. While held, the corner is drawn larger and in the highlight color, returning to normal on
release. Every edge, face and block touching the corner is refitted live, so the block visibly
deforms as you drag. Dragging does not re-classify — run **Classify** afterwards to snap the moved
corners back onto the geometry.

The scene itself stays put while you edit. Polyscope normally recomputes the scene's bounding box
and length scale whenever a structure changes, which drags the ground plane along with it — so biy
freezes both to the model once it's loaded (`options::automaticallyComputeSceneExtents`). Growing a
block, or pulling a corner far out, no longer shifts the ground or rescales the view.

## Configuration

At startup biy reads `biy_config.json` from the current directory. The file is optional, and may
set only the keys it cares about; anything missing keeps its default. A malformed file is reported
and ignored rather than being fatal. [`biy/biy_config.json`](https://github.com/franck-ledoux/gecko/blob/main/biy/biy_config.json)
is a copy of the defaults:

```json
{
  "corner_radius": 0.01,
  "corner_highlight_radius": 0.02,
  "corner_highlight_color": [1.0, 1.0, 1.0],

  "corner_color_unclassified": [0.6, 0.2, 0.85],
  "corner_color_on_vertex": [1.0, 0.9, 0.1],
  "corner_color_on_curve": [0.9, 0.15, 0.15],
  "corner_color_on_surface": [0.15, 0.4, 0.95],
  "corner_color_on_volume": [0.2, 0.75, 0.3],

  "model_transparency": 0.45,
  "show_block_edges": true,
  "block_edge_radius": 0.003,
  "block_edge_color": [0.15, 0.15, 0.15]
}
```

| Key                         | Meaning                                                                                  |
| --------------------------- | ---------------------------------------------------------------------------------------- |
| `corner_radius`             | Size of a block corner at rest                                                           |
| `corner_highlight_radius`   | Size of the corner being dragged                                                         |
| `corner_highlight_color`    | Color of the corner being dragged, RGB in `[0,1]`                                        |
| `corner_color_unclassified` | Color of a corner not classified onto anything                                           |
| `corner_color_on_vertex`    | Color of a corner classified on a model vertex                                           |
| `corner_color_on_curve`     | Color of a corner classified on a model curve                                            |
| `corner_color_on_surface`   | Color of a corner classified on a model surface                                          |
| `corner_color_on_volume`    | Color of a corner classified on a model volume                                           |
| `model_transparency`        | Opacity of the model surface in `[0,1]`; below 1 so the blocking inside it stays visible |
| `show_block_edges`          | Whether block edges are drawn at startup                                                 |
| `block_edge_radius`         | Thickness of the block edges                                                             |
| `block_edge_color`          | Color of the block edges, RGB in `[0,1]`                                                 |

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
