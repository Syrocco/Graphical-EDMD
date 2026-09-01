You can try it [here on github Pages](https://syrocco.github.io/Graphical-EDMD/)

# Graphical Event Driven Molecular Dynamics



Usually, to deal with collisions, one integrates the equation of motion using a small timestep. A notoriously hard task is to deal with object interpenetration. Tricks must be used to supress these overlaps and the resulting dynamics is **always** nonphysical. Artifacts due to these corrections often prevent high densities to be reached. Here, (almost) everything is solved exactly by computing the time before the next collision happening in the system and then advancing the simulation time to this next collision + processing it.
This type of algorithm is extremely fast but you pay the price of having to deal with an asynchronized evolution of events. Since every action is driven by the scheduling of an "event" (every collision in the future, every cell partitioning, every display of the screen in the game loop becomes an event), one has to implement a fast data structure able to quickly sort these events, retrieve them, delete them, reschedule them. This is achieved by a BST and a buffer/priority queue.

At full speed, and for ~2000 particles the algorithm can perform around $10^6$ collisions/sec.

## Compiling

The non-graphical simulation remains a C build:

```sh
make                    # 2D, no graphical interface
make THREE_D=1          # 3D, no graphical interface
```

The graphical build requires `make`, a C99 compiler, a C++17 compiler, and
[raylib](https://github.com/raysan5/raylib/wiki/Working-on-GNU-Linux). Dear
ImGui (docking branch) and rlImGui are pinned and vendored under `src/vendor`,
so they do not need to be installed separately.

```sh
make graphics           # 2D raylib + Dear ImGui interface
make graphics THREE_D=1 # 3D raylib + Dear ImGui interface
./a.out
```

Compile-time physics capabilities can be selected explicitly. They are shown
in the GUI's **Advanced > Build capabilities** table, and each feature set gets
its own object directory:

```sh
make graphics THREE_D=1 STRESS=1 HEAT=1 \
    EDMD_VERIFY_COLLISION_ACCOUNTING=1 \
    EDMD_VERIFY_QUEUE=1 EDMD_VERIFY_OVERLAPS=1
make graphics TANGENTIAL=1 # legacy 2D scalar-spin rough particles
make graphics INTRUDER=1   # legacy 2D polygonal intruder
```

`NONADDITIVE=1` (the historical macro name) selects the usual radius-sum
contact distance; `NONADDITIVE=0` selects the geometric-mean rule. Dimension,
transport, contact, and verification flags are all part of the build-cache
key, preventing objects from an incompatible feature set from being reused.

The Makefile rejects combinations whose implementation is not physically
complete, such as the 2D rough-particle or intruder models in 3D, and those
models combined with `STRESS`/`HEAT`. The fixed-axis tangential delta described
below is a separate runtime law and works in a normal `THREE_D=1` build.

The graphical interface is currently supported on macOS and Linux, not
Windows. On macOS the Makefile discovers a Homebrew raylib installation. The
2D and 3D graphical objects are kept in separate build directories, so the two
commands can be used successively without reusing objects from the wrong
dimension.

### Browser / WebAssembly build

Install and activate the [Emscripten SDK](https://emscripten.org/docs/getting_started/downloads.html)
so that `emcc`, `em++`, `emar`, and `emmake` are on `PATH`. The first build also
needs `curl` or `wget`, `tar`, and a SHA-256 utility. It downloads and verifies
the pinned raylib 5.5 source; Dear ImGui and rlImGui remain the vendored copies.

```sh
make web                    # complete 2D + 3D browser site
make web-2d                 # rebuild only the 2D app
make web-3d                 # rebuild only the 3D app
make web-serve              # build both, then serve http://127.0.0.1:8000/
```

Deploy the entire `docs/` directory. Its root launcher and the self-contained
`docs/2d/` and `docs/3d/` applications use only relative URLs, so the same tree
works at a domain root or a GitHub Pages project path. A small switch in the
bottom-right moves between dimensions without a toolbar or raylib header. The
root also accepts `?dimension=3d` for a direct link to the 3D application.

The application must be served over HTTP(S): browser security prevents a
WebAssembly sidecar from being fetched by a page opened directly with a
`file://` URL. Opening `docs/index.html` that way now explains the limitation
instead of failing silently. Use `make web-serve` for local testing; GitHub
Pages supplies HTTPS and the correct Wasm MIME type. `make web-clean` removes
generated browser outputs, web objects, and the workspace-local Emscripten
cache while retaining the verified raylib source for the next build.

Browser builds use WebGL 2 and a single-threaded animation loop, so they work
on ordinary static hosting such as GitHub Pages without SharedArrayBuffer
headers. As with native graphics builds, each dimension and compile-time
feature set has its own object directory.

The custom page is a full-window canvas with only a transient loading/error
overlay; the old raylib header, logo, controls, and console band are gone. It
runs the same dockable ImGui application as the native graphics target.
Snapshots become browser downloads, restart dumps can be selected from the
local machine, and structure-factor work is spread across animation frames so
it does not freeze the page. On HiDPI displays, the browser canvas keeps its
logical CSS size separate from its physical backing buffer, and scene render
targets use the same live pixel ratio through resizes and browser zoom changes.

## Running

Compile and execute `./a.out`. You can pass options on the command line, such
as `./a.out -N 1000 --phi 0.8`. In the graphical build these values seed the
initial configuration and can then be changed through the interface. The
default model is an elastic, ballistic hard-particle system with no thermostat,
forcing, or square-well potential. Graphical runs are unbounded and remain
under interactive control; headless runs retain the historical default horizon
of `30000` simulation-time units. `--time <t>` overrides either default.

Without a graphical interface the program can write restart/snapshot files
readable by [OVITO](https://www.ovito.org/about/) and the analysis streams
described below. `--load <path>` accepts full dump checkpoints. Partial
initial-state formats are intentionally rejected because they cannot restore
velocities or complete 3D state.

## Graphical interface

The raylib scene is embedded in a dockable Dear ImGui workspace. The default
layout provides the simulation viewport, a searchable Configuration window,
particle Inspector, Diagnostics, Live analysis, and an optional Structure
factor window.
Windows can be moved, resized, tabbed, or redocked; use **Layout > Reset default
workspace** to restore the initial arrangement.

The Configuration window is split into modular **Setup**, **Dynamics**,
**Boundaries**, **Appearance**, and **Advanced** sections. Its
search box filters settings by name and related terms. Options that cannot be
combined, are unavailable in the selected dimension, or are fixed by the build
remain visible but disabled; hover them to see the reason. Red validation
errors prevent an invalid restart, while yellow warnings are advisory.

Configuration changes are hot-reloaded at the smallest safe scope. Collision
response parameters such as restitution update in place. Flight-law changes
preserve the particles and rebuild their predictions; potential range changes
preserve the trajectory while rebuilding only the cell/well index; potential
species and charge edits rebuild only well memberships and affected events;
removing a wall preserves the current state and reschedules boundary events.
Settings which define a new initial state (particle count, initializer,
composition, box sizing, and similar setup fields) reconstruct the system at
`t = 0`. Adding a wall or membrane first tries the preserving path; if an
existing particle intersects the new constraint, the same requested
configuration is applied automatically by rebuilding a clean system at `t = 0`.

Numeric sliders and drags send these updates at a bounded live-preview rate
while the pointer is held, followed by a final update on release. Exact numeric
entry waits until editing finishes and then uses a short debounce. Appearance
settings take effect on the current frame without rebuilding the event queue.
Actually invalid or mutually incompatible configurations are not applied; the
previous valid state keeps running and the nonmodal status window explicitly
reports the rollback. Valid edits which are incompatible only with the current
trajectory instead trigger the automatic clean rebuild described above, with a
nonmodal explanation. There is no Apply/Revert step. The explicit **Restart**
action repeats the last successfully loaded configuration.

Run/pause, single-event stepping, snapshot capture, and optional physical-time
pacing are available in the top toolbar. **Pace** ties simulation-time progress
to a chosen number of simulation-time units per wall-clock second. Its
logarithmic slider covers interactive slow motion through fast playback and an
adjacent exact entry accepts a wider rate; neither changes the physical system.
The Diagnostics window reports accumulated lag when event processing cannot
keep up. Appearance controls are live and do not restart the simulation:
visibility of boundaries, axes, velocity arrows, selection and FPS overlays;
particle colors and palettes; target FPS; scene supersampling; color refresh;
3D field of view and sphere tessellation quality; and structure-factor sample
cap and refresh rate. Automatic sphere quality is useful for dense 3D systems.

Viewport controls:

- Mouse wheel: zoom.
- Right-drag: pan in 2D or orbit in 3D.
- `R`: restart the simulation; `Home` or `C`: refit the active camera.
- `Space`: run or pause; `F`: toggle fullscreen; `S`: save a snapshot.
- Click a particle to inspect it. Left-dragging applies a damped spring-like
  force toward the pointer in the simulation plane (or a camera-facing plane
  in 3D). It does not teleport the particle or automatically pause the run;
  affected collision and transport events are rescheduled continuously.
- While paused, a selected particle's velocity can be edited and the particle
  is rescheduled. Mass is also editable unless `STRESS` or `HEAT` is compiled,
  where it is disabled because transport accounting depends on the existing
  mass.

An overlap, invalid collision prediction, or exhausted event queue pauses the
graphical simulation instead of terminating the process. The status window is
nonmodal: configuration, inspection, plots, and restart remain usable, but Run,
Space, and single-step stay locked until a successful reconstruction clears the
fault.

The 3D renderer uses perspective, depth testing, true sphere meshes, and
complete velocity-arrow geometry. Its live lighting selector offers the
inexpensive original shading, balanced studio key/fill/hemisphere lighting,
or a single-pass depth-enhanced profile with subtle silhouette, rim, and
distance cues. The latter is deliberately lightweight rather than claiming
screen-space ambient occlusion or ray tracing. Conservative camera-frustum
culling skips off-screen sphere instances, velocity arrows, and selection
outlines; Diagnostics reports the visible and culled particle counts. Dynamic
clip planes keep large boxes visible without sacrificing ordinary-scene depth
precision. Random growth and hexagonal
close-packed (HCP) lattice starts are available, as are periodic boundaries,
individual X/Y planar walls, an X/Y wall box, a fully closed X/Y/Z wall box,
and spherical confinement. Lees-Edwards shear and fixed-axis tangential delta
kicks (`u_z x n`, where `n` is the contact normal) are dimension-complete.
Analyses that are inherently planar, such as the live 3D structure factor, are
labeled as XY projections. On HiDPI displays the window and scene target use
the physical framebuffer resolution; the Appearance supersampling value is an
additional quality multiplier, with `1` meaning native display-pixel
resolution.

The dockable Live analysis window keeps bounded in-memory histories for energy
per particle, temperature, pressure, and mean-squared displacement. Each curve
has an independent safe log-log view; nonpositive or nonfinite values are
skipped instead of producing invalid coordinates. MSD starts a fresh reference
when enabled or after a configuration generation changes. The window also
provides an online-averaged radial distribution `g(r)`, with linear or log-log
axes and a configurable particle cap, plus the existing asynchronously
computed structure factor. The 3D structure-factor view remains an XY
projection; `g(r)` uses 3D shell-volume normalization. These plots do not open
files or schedule simulation-time sampling events; particle caps keep their CPU
cost explicit and bounded. **Save Snapshot** remains an explicit manual action.

Headless runs retain persistent analysis. `--dtimeThermo <dt>` controls the
sampling interval; `--msd`, `--pair-correlation`,
`--bond-order-pair-correlation`, `--g6-pair-correlation`,
`--structure-factor density|velocity`, and `--velocity-distribution` enable the
corresponding output streams. Periodic full dumps remain available to the
headless pipeline and can be used as restart/checkpoint files.

The fixed-axis tangential drive adds angular momentum with one sign. With a
perfectly specular circular/spherical wall (`resT = 1`) there is no tangential
torque sink, so sufficiently strong sustained driving can spin the entire
system up instead of reaching a steady state. The validator and GUI warn about
this combination; use tangential wall dissipation (`resT < 1`), a thermal wall,
or damping when a stationary confined state is required.

Potential simulations use a cutoff-aware, radius-aware cell stencil and
ordered, duplicate-free in-well neighbor lists. A zero-depth potential takes
the hard-particle fast path, and pairs excluded by the selected species/charge
rule no longer generate redundant square-well entry/exit events. The event
calendar skips empty time buckets directly and tracks its occupancy without repeatedly scanning every particle.
The graphical host bounds event work per frame, reuses render buffers,
adaptively lowers sphere tessellation for large systems, and caps/cancels
structure-factor jobs. These choices prevent run-time data structures from
growing with simulation age and keep input and rendering responsive without
changing event ordering.

## Structure of the program

#### HUGE DISCLAIMER: I AM NOT --BY ANY STANDARD-- A PROGRAMMER. DO NOT USE MY PROGRAM AS A STANDARD FOR THIS KIND OF SIMULATION/GRAPHICAL DISPLAY
### Technical part
The simulation core is C and the graphical host adds raylib plus a small C++17
Dear ImGui layer.

- Particles are stored in cell-lists by an array based doubly linked list in order to efficiently compute the time before the next collision.
- Events (particles crossing cell-lists boundaries, particles colliding, thermostating, displaying, taking user inputs) are stored in a Binary Search Tree which receive events from [a partially sorted buffer/priority queue](https://arxiv.org/abs/physics/0606226).

### Graphical part
The interface communicates with the legacy engine through a dimension-neutral
configuration schema and a classified reconfiguration lifecycle. It separates
scalar, event-calendar, spatial-index, and genuinely state-defining changes so
the renderer does not discard a trajectory unnecessarily. GUI panels are
separate modules from the raylib renderer, making new settings and panels easier
to add without coupling ImGui widgets directly to the event-processing loop.

## Sources

If you want to understand the structure of an EDMD algorithm: [Github repo EDMD](https://github.com/FSmallenburg/EDMD) with [research article](https://arxiv.org/abs/2201.01100).
