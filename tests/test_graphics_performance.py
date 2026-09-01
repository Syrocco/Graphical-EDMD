#!/usr/bin/env python3
"""Deterministic regression checks for graphics/performance safeguards.

The graphical host needs a display, and throughput varies substantially across
drivers and machines.  These tests therefore verify the structural contracts
that bound per-frame work and prevent time-dependent data-structure growth;
they intentionally make no wall-clock performance assertions.
"""

from __future__ import annotations

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
EDMD_SOURCE = (ROOT / "src" / "EDMD.c").read_text(encoding="utf-8")
GRAPHICS_SOURCE = (ROOT / "src" / "graphics.c").read_text(encoding="utf-8")
SIM_CONFIG_SOURCE = (ROOT / "src" / "sim_config.c").read_text(
    encoding="utf-8"
)
GUI_SOURCE = (ROOT / "src" / "gui_imgui.cpp").read_text(encoding="utf-8")
GUI_HEADER_SOURCE = (ROOT / "src" / "gui_imgui.h").read_text(encoding="utf-8")
GUI_PANELS_SOURCE = (ROOT / "src" / "gui_panels.cpp").read_text(
    encoding="utf-8"
)
GUI_ANALYSIS_SOURCE = (ROOT / "src" / "gui_analysis.cpp").read_text(
    encoding="utf-8"
)
LIVE_ANALYSIS_SOURCE = (ROOT / "src" / "live_analysis.c").read_text(
    encoding="utf-8"
)
WEB_SHELL_SOURCE = (ROOT / "web" / "shell.html").read_text(encoding="utf-8")
WEB_LAUNCHER_SOURCE = (ROOT / "web" / "launcher.html").read_text(
    encoding="utf-8"
)
MAKEFILE = (ROOT / "Makefile").read_text(encoding="utf-8")


def compact(text: str) -> str:
    """Collapse formatting so checks survive harmless indentation changes."""

    return re.sub(r"\s+", " ", text).strip()


def function_body(source: str, name: str) -> str:
    """Return a C function definition, ignoring braces in comments/literals."""

    definition = re.compile(
        rf"\b{re.escape(name)}\s*\([^;{{}}]*\)\s*\{{", re.MULTILINE
    ).search(source)
    if definition is None:
        raise AssertionError(f"function {name}() was not found")

    opening = source.find("{", definition.start())
    depth = 0
    state = "code"
    escaped = False
    index = opening
    while index < len(source):
        char = source[index]
        following = source[index + 1] if index + 1 < len(source) else ""

        if state == "line-comment":
            if char == "\n":
                state = "code"
        elif state == "block-comment":
            if char == "*" and following == "/":
                state = "code"
                index += 1
        elif state in {"string", "character"}:
            delimiter = '"' if state == "string" else "'"
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == delimiter:
                state = "code"
        else:
            if char == "/" and following == "/":
                state = "line-comment"
                index += 1
            elif char == "/" and following == "*":
                state = "block-comment"
                index += 1
            elif char == '"':
                state = "string"
            elif char == "'":
                state = "character"
            elif char == "{":
                depth += 1
            elif char == "}":
                depth -= 1
                if depth == 0:
                    return source[opening + 1 : index]
        index += 1

    raise AssertionError(f"function {name}() has no matching closing brace")


class GraphicsPerformanceContracts(unittest.TestCase):
    def assert_in_order(self, text: str, *needles: str) -> None:
        cursor = 0
        for needle in needles:
            location = text.find(needle, cursor)
            self.assertNotEqual(
                location, -1, f"expected {needle!r} after offset {cursor}"
            )
            cursor = location + len(needle)

    def test_web_pipeline_is_a_relative_dual_dimension_static_site(self) -> None:
        web_rules = compact(
            MAKEFILE[MAKEFILE.index(".PHONY: web web-all") :]
        )
        self.assertIn("web: web-all", web_rules)
        self.assert_in_order(
            web_rules,
            '$(MAKE) web-2d WEB_OUTPUT_DIR="$(WEB_OUTPUT_DIR)"',
            '$(MAKE) web-3d WEB_OUTPUT_DIR="$(WEB_OUTPUT_DIR)"',
            '$(MAKE) web-launcher WEB_OUTPUT_DIR="$(WEB_OUTPUT_DIR)"',
        )
        self.assertIn('WEB_OUTPUT_DIR="$(WEB_OUTPUT_DIR)/2d"', web_rules)
        self.assertIn('WEB_OUTPUT_DIR="$(WEB_OUTPUT_DIR)/3d"', web_rules)
        self.assertIn('touch "$(WEB_OUTPUT_DIR)/.nojekyll"', web_rules)
        self.assertIn("web-serve: web", web_rules)

        # Static GitHub Pages deployments must not depend on a domain-root URL
        # or on cross-origin-isolation headers.
        self.assertIn('href="./2d/"', WEB_LAUNCHER_SOURCE)
        self.assertIn('href="./3d/"', WEB_LAUNCHER_SOURCE)
        self.assertNotRegex(WEB_LAUNCHER_SOURCE, r'(?:href|src)="/')
        self.assertIn('params.get("dimension")', WEB_LAUNCHER_SOURCE)
        self.assertIn("target.hash = window.location.hash", WEB_LAUNCHER_SOURCE)
        self.assertIn('window.location.protocol === "file:"', WEB_LAUNCHER_SOURCE)

        link_flags = compact(
            MAKEFILE[
                MAKEFILE.index("WEB_LINK_FLAGS =") : MAKEFILE.index("EDMD:")
            ]
        )
        self.assertNotIn("PTHREADS", link_flags)
        self.assertIn("-sENVIRONMENT=web", link_flags)

        self.assertIn('class="dimension-switch"', WEB_SHELL_SOURCE)
        self.assertIn("bottom: max(10px", WEB_SHELL_SOURCE)
        self.assertIn("opacity: 0.42", WEB_SHELL_SOURCE)
        self.assertIn('aria-label="Simulation dimension"', WEB_SHELL_SOURCE)
        self.assertIn("target.search = params.toString()", WEB_SHELL_SOURCE)
        self.assertIn("target.hash = window.location.hash", WEB_SHELL_SOURCE)
        self.assertIn('window.location.protocol === "file:"', WEB_SHELL_SOURCE)

    def test_event_processing_has_a_wall_clock_frame_budget(self) -> None:
        # main() contains mutually exclusive #if G branches whose raw braces do
        # not balance until preprocessing, so inspect the bounded-loop region
        # instead of trying to parse the whole unpreprocessed function.
        budget_start = EDMD_SOURCE.find("simulationFrameStart = GetTime()")
        self.assertNotEqual(budget_start, -1)
        draw_after_budget = EDMD_SOURCE.find("draw(argc", budget_start)
        self.assertNotEqual(draw_after_budget, -1)
        bounded_loop = compact(
            EDMD_SOURCE[budget_start : draw_after_budget + len("draw(argc")]
        )
        self.assertRegex(
            bounded_loop,
            r"GetTime\(\)\s*-\s*simulationFrameStart\s*<\s*1\.0\s*/\s*180\.0",
        )
        self.assert_in_order(
            bounded_loop, "simulationFrameStart = GetTime()", "draw(argc"
        )

    def test_optional_simulation_time_pacing_is_exact_and_reports_lag(self) -> None:
        # The host advances continuous event-free time only to the current
        # wall-clock target; a later queued event remains untouched.
        pacing_gate = compact(
            EDMD_SOURCE[
                EDMD_SOURCE.find("graphicsSimulationTimeLimit(t)") - 120 :
                EDMD_SOURCE.find("graphicsSimulationTimeLimit(t)") + 520
            ]
        )
        self.assert_in_order(
            pacing_gate,
            "graphicsSimulationTimeLimit(t)",
            "nextEvent->t > pacedTimeLimit",
            "t = pacedTimeLimit",
            "break;",
        )
        pacing = compact(
            function_body(GRAPHICS_SOURCE, "graphicsSimulationTimeLimit")
        )
        self.assertIn("synchronizeSimulationPacing(true, currentSimulationTime)", pacing)
        self.assertIn("return HUGE_VAL", pacing)
        self.assertIn("fmax(currentSimulationTime, pacingState.targetSimulationTime)", pacing)

        toolbar = compact(function_body(GUI_SOURCE, "draw_toolbar"))
        self.assertIn("pace_simulation_time", toolbar)
        self.assertIn("target_simulation_time_rate", toolbar)
        self.assertIn("simulation_time_lag", toolbar)
        self.assertIn('ImGui::SliderScalar("##pace-rate-slider"', toolbar)
        self.assertIn("ImGuiSliderFlags_Logarithmic", toolbar)
        self.assertIn("ImGuiSliderFlags_AlwaysClamp", toolbar)
        self.assertIn('ImGui::InputDouble("##pace-rate-input"', toolbar)
        self.assertIn("minimum_slider_rate = 1e-3", toolbar)
        self.assertIn("maximum_slider_rate = 1e3", toolbar)

        # Pacing is a live host concern: changing either rate widget must not
        # enqueue a hot reload or restart of the physical simulation.
        pacing_controls = toolbar[
            toolbar.index('ImGui::Checkbox("Pace"') :
            toolbar.index('ImGui::Button("Step event"')
        ]
        self.assertNotIn("assign_command", pacing_controls)
        self.assertNotIn("assign_restart", pacing_controls)

    def test_live_plots_offer_safe_per_series_log_log_axes(self) -> None:
        plottable = compact(
            function_body(GUI_ANALYSIS_SOURCE, "series_sample_is_plottable")
        )
        self.assertIn("sample_time > 0.0 && sample_value > 0.0", plottable)
        coordinate = compact(
            function_body(GUI_ANALYSIS_SOURCE, "plot_coordinate")
        )
        self.assertIn("std::log10(value)", coordinate)
        plot = compact(function_body(GUI_ANALYSIS_SOURCE, "draw_series_plot"))
        self.assertIn('ImGui::Checkbox("log-log", log_log)', plot)
        self.assertIn("previous_valid = false", plot)
        self.assertIn("series_sample_is_plottable", plot)
        live = compact(function_body(GUI_ANALYSIS_SOURCE, "draw_live_analysis"))
        for field in (
            "energy_log_log",
            "temperature_log_log",
            "pressure_log_log",
            "msd_log_log",
            "pair_correlation_log_log",
        ):
            self.assertIn(field, live)
        pair_plot = compact(function_body(GUI_ANALYSIS_SOURCE, "draw_pair_plot"))
        self.assertIn('ImGui::Checkbox("log-log", log_log)', pair_plot)
        self.assertIn("pair_sample_is_plottable", pair_plot)
        self.assertIn("previous_valid = false", pair_plot)
        self.assertIn("logarithmic color scale", GUI_ANALYSIS_SOURCE)

    def test_gui_has_no_sampling_tab_or_cadence_tuning_controls(self) -> None:
        self.assertNotIn('BeginTabItem("Sampling")', GUI_PANELS_SOURCE)
        self.assertNotIn('collapsing_header(context, "Sampling"', GUI_PANELS_SOURCE)
        self.assertNotIn("Plot refresh", GUI_ANALYSIS_SOURCE)
        self.assertNotIn("Pair refresh", GUI_ANALYSIS_SOURCE)
        self.assertNotIn("Structure-factor refresh", GUI_ANALYSIS_SOURCE)
        self.assertNotIn('"Refresh interval", "structure factor', GUI_PANELS_SOURCE)

    def test_failed_live_update_explicitly_reports_transactional_rollback(self) -> None:
        handler = compact(function_body(GRAPHICS_SOURCE, "handleGuiCommand"))
        self.assertIn("simulation_reconfigure_core", handler)
        self.assertIn("SIM_RECONFIGURE_FAILED", handler)
        self.assertIn("Live update rolled back", handler)
        self.assertIn("GUI_CONFIG_APPLY_REJECTED", handler)
        self.assertIn("The previous valid configuration and simulation", handler)
        self.assertIn("state were retained", handler)
        self.assertIn("no partial changes were applied", handler)

    def test_live_update_preserves_renderer_state_unless_restart_is_required(self) -> None:
        handler = compact(function_body(GRAPHICS_SOURCE, "handleGuiCommand"))
        in_place = handler[
            handler.index("if (action != SIM_RECONFIGURE_RESTART_REQUIRED)") :
            handler.index("/* State-defining edits and runtime constraints")
        ]
        self.assertIn("GUI_CONFIG_APPLY_SUCCEEDED", in_place)
        self.assertNotIn("clearParticleSelection", in_place)
        self.assertNotIn("resetCamera", in_place)
        self.assertNotIn("simulation_restart_core", in_place)

        restart_fallback = handler[
            handler.index("/* State-defining edits and runtime constraints") :
            handler.index("case GUI_COMMAND_RESTART_APPLIED")
        ]
        self.assertIn("simulation_restart_core", restart_fallback)
        self.assertIn("GUI_CONFIG_APPLY_RESTARTED", restart_fallback)
        self.assertIn("Configuration applied by rebuilding", restart_fallback)
        self.assertIn("restartNotice", restart_fallback)

    def test_nonnegative_physical_drags_are_bounded(self) -> None:
        drag = compact(function_body(GUI_PANELS_SOURCE, "input_double_row"))
        self.assertIn("&drag_minimum, &drag_maximum", drag)
        self.assertIn("ImGuiSliderFlags_AlwaysClamp", drag)
        self.assertIn("std::clamp(*value, drag_minimum, drag_maximum)", drag)
        nonnegative = compact(
            function_body(GUI_PANELS_SOURCE, "input_nonnegative_double_row")
        )
        self.assertIn("0.0, DBL_MAX", nonnegative)
        slider_int = compact(function_body(GUI_PANELS_SOURCE, "slider_int_row"))
        self.assertIn("std::clamp(*value, minimum, maximum)", slider_int)
        slider_double = compact(
            function_body(GUI_PANELS_SOURCE, "slider_double_row")
        )
        self.assertIn("std::clamp(*value, minimum, maximum)", slider_double)
        advanced = compact(function_body(GUI_PANELS_SOURCE, "draw_advanced"))
        for label in (
            "Side length",
            "Inner radius",
            "Mass",
            "Collision delta",
            "Bath particle radius",
            "Kinetic-model burn time",
        ):
            self.assertRegex(
                advanced,
                rf"input_nonnegative_double_row\([^;]*{re.escape(label)}",
            )

    def test_gui_exposes_the_full_validated_numeric_ranges(self) -> None:
        setup = compact(function_body(GUI_PANELS_SOURCE, "draw_setup"))
        self.assertIn("triangular_close_packing = 0.90689968211710892530", setup)
        self.assertIn("hcp_close_packing = 0.74048048969306104117", setup)
        self.assertRegex(
            setup,
            r"packing_slider_max\s*=\s*init\.mode\s*==\s*"
            r"SIM_INIT_HEX_LATTICE\s*\?\s*hex_close_packing\s*:\s*1\.0",
        )

        appearance = compact(
            function_body(GUI_PANELS_SOURCE, "draw_appearance")
        )
        self.assertIn("&live.target_fps, 1, 1000", appearance)
        self.assertIn(
            "&live.structure_factor_sample_cap, 1, 8192", appearance
        )

    def test_potential_cell_stencil_is_cutoff_aware_and_deduplicated(self) -> None:
        box = compact(function_body(EDMD_SOURCE, "boxConstantHelper"))
        self.assertRegex(
            box,
            r"targetCellSize\s*=\s*2\s*\*\s*sig\s*\*\s*"
            r"collisionGridMaxRadius\s*/\s*collisionCellNeighborRange",
        )
        self.assertIn("Nxcells = (int)(Lx/targetCellSize)", box)
        self.assertIn("Nycells = (int)(Ly/targetCellSize)", box)

        neighbors = compact(
            function_body(EDMD_SOURCE, "neighboringCellCoordinates")
        )
        self.assertIn("if (!duplicate) coordinates[number++] = coordinate", neighbors)
        self.assertIn("coordinates[i] == coordinate", neighbors)

    def test_potential_membership_cannot_grow_from_duplicate_entries(self) -> None:
        add = compact(function_body(EDMD_SOURCE, "addParticleInWellList"))
        self.assert_in_order(
            add,
            "particleWellLowerBound(p1, num)",
            "p1->particlesInWell[lower] == num",
            "return;",
            "p1->numberOfParticlesInWell == p1->particlesInWellCapacity",
            "p1->numberOfParticlesInWell++",
        )

        remove = compact(function_body(EDMD_SOURCE, "removeParticleInWellList"))
        self.assertIn("particleWellLowerBound(p1, num)", remove)
        self.assertIn("p1->numberOfParticlesInWell--", remove)
        self.assertIn("memmove", remove)

        lookup = compact(function_body(EDMD_SOURCE, "isParticleInWellList"))
        self.assertIn("numberOfParticlesInWell <= 16", lookup)
        self.assertIn("particleWellLowerBound(p1, num)", lookup)

    def test_noninteracting_species_skip_square_well_events(self) -> None:
        pairs = compact(function_body(EDMD_SOURCE, "pairsInit"))
        self.assertIn("potentialPairStepEnergy(p1, p2) != 0.0", pairs)

        scheduler = compact(
            function_body(EDMD_SOURCE, "scheduledPotentialEventTime")
        )
        self.assertIn(
            "U == 0.0 || potentialPairStepEnergy(p1, p2) == 0.0",
            scheduler,
        )
        self.assert_in_order(
            scheduler,
            "potentialPairStepEnergy(p1, p2) == 0.0",
            "*eventType = COLLISION",
            "scheduledCollisionTimeNormal(p1, p2)",
        )

    def test_particle_and_gpu_buffers_are_reused_at_stable_capacity(self) -> None:
        predictions = compact(
            function_body(GRAPHICS_SOURCE, "refreshRenderPredictions")
        )
        self.assert_in_order(
            predictions,
            "if (renderPredictionCapacity != N)",
            "realloc(renderPredictions",
            "renderPredictionCapacity = N",
        )

        circles = compact(function_body(GRAPHICS_SOURCE, "ensureCircleBatch"))
        self.assert_in_order(
            circles,
            "batch->loaded && (batch->capacity == capacity)",
            "return true;",
            "unloadCircleBatch(batch)",
            "loadCircleBatch(batch, capacity)",
        )

        arrows = compact(
            function_body(GRAPHICS_SOURCE, "ensureVelocityArrowBatch")
        )
        self.assert_in_order(
            arrows,
            "batch->loaded && batch->capacity == capacity",
            "return true;",
            "unloadVelocityArrowBatch(batch)",
            "loadVelocityArrowBatch(batch, capacity)",
        )

        spheres = compact(function_body(GRAPHICS_SOURCE, "ensureSphereBatch"))
        self.assert_in_order(
            spheres,
            "batch->loaded && (batch->capacity == capacity)",
            "batch->rings == rings",
            "batch->slices == slices",
            "return true;",
            "unloadSphereBatch(batch)",
            "loadSphereBatch(batch, capacity, rings, slices)",
        )

        target = compact(function_body(GRAPHICS_SOURCE, "ensureSceneTarget"))
        self.assert_in_order(
            target,
            "sceneTargetWidth == width",
            "sceneTargetHeight == height",
            "return true;",
            "UnloadRenderTexture(sceneTarget)",
            "LoadRenderTexture(width, height)",
        )

    def test_2d_particles_are_procedural_batched_circles(self) -> None:
        loader = compact(
            function_body(GRAPHICS_SOURCE, "loadProceduralCircleMaterial")
        )
        self.assertIn("LoadShaderFromMemory", loader)
        self.assertIn("circleEdgeScale", loader)
        self.assertIn("immediate circles", loader)
        self.assertIn("fwidth(distanceSquared)", GRAPHICS_SOURCE)
        self.assertIn("smoothstep", GRAPHICS_SOURCE)

        draw = compact(function_body(GRAPHICS_SOURCE, "drawParticlesAndBox"))
        self.assertIn("appendCircleQuad", draw)
        self.assertIn("drawCircleBatch(&particleCircleBatch)", draw)
        self.assertIn("DrawCircleV", draw)
        self.assertEqual(
            compact(function_body(GRAPHICS_SOURCE, "drawCircleBatch")).count(
                "DrawMesh("
            ),
            1,
        )
        self.assertNotIn("tex_data.h", GRAPHICS_SOURCE)
        self.assertNotIn("circleTexture", GRAPHICS_SOURCE)
        self.assertNotIn("LoadImageFromMemory", GRAPHICS_SOURCE)

    def test_default_scene_palette_matches_dark_imgui_canvas(self) -> None:
        expected = ("4296fa", "f29d49", "181a1f", "8b949e")
        for rgb in expected:
            self.assertIn(rgb, GRAPHICS_SOURCE.lower())
            self.assertIn(rgb, SIM_CONFIG_SOURCE.lower())

        selection = compact(function_body(GRAPHICS_SOURCE, "drawParticlesAndBox"))
        self.assertIn("selectionColor", selection)
        arrows = compact(function_body(GRAPHICS_SOURCE, "appendVelocityArrow2D"))
        self.assertIn("236, 239, 244", arrows)

    def test_3d_spheres_are_instanced_with_adaptive_mesh_complexity(self) -> None:
        resolution = function_body(GRAPHICS_SOURCE, "sphereResolution3D")
        automatic_levels = [
            tuple(map(int, match))
            for match in re.findall(
                r"particleCount\s*<=\s*(\d+)\s*\)\s*\{\s*"
                r"\*rings\s*=\s*(\d+)\s*;\s*\*slices\s*=\s*(\d+)",
                resolution,
                re.DOTALL,
            )
        ]
        self.assertEqual(
            automatic_levels,
            [(2500, 16, 24), (5000, 12, 16), (10000, 8, 12), (20000, 6, 8)],
        )
        triangle_proxies = [rings * slices for _, rings, slices in automatic_levels]
        self.assertEqual(triangle_proxies, sorted(triangle_proxies, reverse=True))
        self.assertRegex(resolution, r"else\s*\{\s*\*rings\s*=\s*4\s*;\s*\*slices\s*=\s*6")
        self.assertIn("GenMeshSphere", function_body(GRAPHICS_SOURCE, "loadSphereBatch"))
        self.assertIn("DrawMeshInstanced", function_body(GRAPHICS_SOURCE, "drawSphereBatch"))

    def test_velocity_arrows_are_directional_speed_aware_geometry(self) -> None:
        shape_2d = compact(function_body(GRAPHICS_SOURCE, "velocityArrowShape"))
        self.assertIn("log1p(speed)", shape_2d)
        self.assertIn("direction", shape_2d)
        self.assertIn("outerShaftHalfWidth", shape_2d)
        self.assertIn("innerHeadHalfWidth", shape_2d)

        append_2d = compact(
            function_body(GRAPHICS_SOURCE, "appendVelocityArrow2D")
        )
        self.assertEqual(append_2d.count("appendVelocityArrowLayer"), 2)
        self.assertIn("outline", append_2d)
        self.assertIn("fill", append_2d)

        draw_3d = compact(function_body(GRAPHICS_SOURCE, "drawVelocityArrow3D"))
        self.assertIn("log1p(speed)", draw_3d)
        self.assertEqual(draw_3d.count("DrawCylinderEx"), 2)
        self.assertIn("shaftEnd", draw_3d)
        self.assertIn("tip", draw_3d)
        self.assertNotIn("arrowTexture", GRAPHICS_SOURCE)

    def test_structure_factor_work_is_capped_reused_and_cancellable(self) -> None:
        hard_cap = re.search(
            r"#define\s+STRUCTURE_FACTOR_HARD_SAMPLE_CAP\s+(\d+)",
            GRAPHICS_SOURCE,
        )
        self.assertIsNotNone(hard_cap)
        self.assertLessEqual(int(hard_cap.group(1)), 8192)

        desired = compact(
            function_body(GRAPHICS_SOURCE, "desiredStructureFactorSampleCount")
        )
        self.assertIn("(population < sampleCap) ? population : sampleCap", desired)

        allocate = compact(
            function_body(GRAPHICS_SOURCE, "allocateStructureFactorPositions")
        )
        self.assert_in_order(
            allocate,
            "desiredCount <= structSampleCapacity",
            "return true;",
            "realloc(positions",
        )

        release = compact(
            function_body(GRAPHICS_SOURCE, "releaseStructureFactorPositions")
        )
        self.assert_in_order(
            release,
            "structJobCancelRequested = true",
            "awaitStructFactor()",
            "free(positions)",
        )
        worker = compact(
            function_body(GRAPHICS_SOURCE, "computeStructureFactorGraphic")
        )
        row_kernel = compact(
            function_body(GRAPHICS_SOURCE, "computeStructureFactorRow")
        )
        self.assertIn("computeStructureFactorRow(i, completedGeneration)", worker)
        self.assertIn("(n & 63) == 0", row_kernel)
        self.assertIn("structureFactorJobWasCancelled(generation)", row_kernel)

    def test_msd_only_does_not_walk_potential_memberships(self) -> None:
        update = compact(
            function_body(LIVE_ANALYSIS_SOURCE, "live_analysis_update")
        )
        self.assertIn(
            "energetics_requested = state.settings.energy || "
            "state.settings.temperature || state.settings.pressure",
            update,
        )
        self.assertRegex(
            update,
            r"if \(energetics_requested\) sample_energetics\(",
        )

    def test_render_target_resizing_is_quantized_and_bounded(self) -> None:
        sizing = compact(
            function_body(GRAPHICS_SOURCE, "requestedSceneTargetSize")
        )
        self.assertIn("logicalWidth*dpiX*supersampling", sizing)
        self.assertIn("logicalHeight*dpiY*supersampling", sizing)
        self.assertIn("/16.0f", sizing)
        self.assertIn("*width > 6144", sizing)
        self.assertIn("*height > 6144", sizing)

    def test_global_shortcuts_survive_docking_keyboard_capture(self) -> None:
        frame = compact(function_body(GUI_SOURCE, "gui_imgui_frame"))
        self.assertIn("!ImGui::GetIO().WantTextInput", frame)
        self.assertIn("!ImGui::IsAnyItemActive()", frame)
        self.assertRegex(
            frame,
            r"ImGui::IsKeyPressed\(ImGuiKey_Space, false\).*?"
            r"live->running = !live->running",
        )
        self.assertRegex(
            frame,
            r"ImGui::IsKeyPressed\(ImGuiKey_R, false\).*?"
            r"assign_restart\(result.command, g_context->applied, live->running\)",
        )
        self.assertRegex(
            frame,
            r"ImGuiKey_Home.*?ImGuiKey_C.*?GUI_COMMAND_RESET_VIEW",
        )

    def test_numeric_sliders_live_update_while_held_and_exact_entry_debounces(self) -> None:
        frame = compact(function_body(GUI_SOURCE, "gui_imgui_frame"))
        self.assertIn("automatic_apply_pending", frame)
        self.assertIn("reactive_numeric_active", frame)
        self.assertIn("reactive_numeric_deactivated_after_edit", frame)
        self.assertIn("kReactiveApplyIntervalSeconds", frame)
        # A release after the last held preview must not reapply the same
        # already-applied value a second time.
        self.assertIn("const bool reactive_numeric_apply = dirty &&", frame)
        self.assertIn("reactive_numeric_apply || settled_apply", frame)
        self.assertIn("result.command.type == GUI_COMMAND_NONE", frame)
        self.assertIn("now - g_context->last_physics_edit_time >= 0.30", frame)
        self.assertIn("!ImGui::IsAnyItemActive()", frame)
        self.assert_in_order(
            frame,
            "live->running = !live->running",
            "assign_config_apply(result.command, g_context->draft, live->running)",
        )
        apply = compact(function_body(GUI_SOURCE, "assign_config_apply"))
        self.assertIn("GUI_COMMAND_APPLY_CONFIG", apply)
        self.assertNotIn("GUI_COMMAND_RESTART_APPLIED", apply)

        for helper, pointer_widget, exact_widget in (
            ("slider_int_row", "ImGui::SliderInt", "ImGui::InputInt"),
            ("slider_double_row", "ImGui::SliderScalar", "ImGui::InputDouble"),
            ("input_int_row", "ImGui::DragInt", "ImGui::InputInt"),
            ("input_double_row", "ImGui::DragScalar", "ImGui::InputDouble"),
        ):
            body = compact(function_body(GUI_PANELS_SOURCE, helper))
            self.assert_in_order(
                body,
                pointer_widget,
                "track_reactive_numeric_edit(context, reactive)",
                exact_widget,
            )

    def test_apply_acknowledgement_distinguishes_live_update_from_restart(self) -> None:
        self.assertIn("GUI_COMMAND_APPLY_CONFIG", GUI_HEADER_SOURCE)
        self.assertNotIn("GUI_COMMAND_APPLY_RESTART", GUI_HEADER_SOURCE)
        acknowledge = compact(
            function_body(GUI_SOURCE, "gui_imgui_acknowledge_applied_config")
        )
        succeeded = acknowledge[
            acknowledge.index("case GUI_CONFIG_APPLY_SUCCEEDED:") :
            acknowledge.index("case GUI_CONFIG_APPLY_REJECTED:")
        ]
        rejected = acknowledge[
            acknowledge.index("case GUI_CONFIG_APPLY_REJECTED:") :
            acknowledge.index("case GUI_CONFIG_APPLY_RESTARTED:")
        ]
        restarted = acknowledge[
            acknowledge.index("case GUI_CONFIG_APPLY_RESTARTED:") :
            acknowledge.index("default:")
        ]

        self.assertIn("g_context->applied = *applied_config", succeeded)
        self.assertIn("g_context->draft = *applied_config", succeeded)
        self.assertNotIn("simulation_faulted = false", succeeded)
        self.assertNotIn("live_analysis_reset_generation", succeeded)
        self.assertNotIn("particle_edit", succeeded)

        self.assertIn("automatic_apply_blocked = true", rejected)
        self.assertNotIn("g_context->draft = *applied_config", rejected)
        self.assertIn("simulation_faulted = false", restarted)
        self.assertIn("live_analysis_reset_generation", restarted)

        legacy = compact(function_body(GUI_SOURCE, "gui_imgui_set_applied_config"))
        self.assertIn("GUI_CONFIG_APPLY_RESTARTED", legacy)

        help_window = compact(function_body(GUI_SOURCE, "draw_help_window"))
        self.assertIn("sent as live updates", help_window)
        self.assertIn("Sliders and drags send live updates while held", help_window)

    def test_successful_rebuild_notice_is_a_short_fading_toast(self) -> None:
        acknowledge = compact(
            function_body(GUI_SOURCE, "gui_imgui_acknowledge_applied_config")
        )
        self.assertIn("outcome == GUI_CONFIG_APPLY_RESTARTED", acknowledge)
        self.assertIn(
            "Simulation rebuilt for the new configuration", acknowledge
        )
        self.assertIn("toast_started_at = ImGui::GetTime()", acknowledge)
        self.assertIn("toast_visible = true", acknowledge)

        toast = compact(function_body(GUI_SOURCE, "draw_rebuild_toast"))
        self.assertIn("hold_seconds = 2.5", toast)
        self.assertIn("fade_seconds = 3.5", toast)
        self.assertIn("viewport->WorkPos.x + 16.0f", toast)
        self.assertIn("ImGuiWindowFlags_NoInputs", toast)
        self.assertIn("context.toast_visible = false", toast)
        self.assertNotIn("context.runtime_error", toast)

        # Rollbacks and engine faults retain the persistent status window.
        status = compact(function_body(GUI_SOURCE, "draw_runtime_status"))
        self.assertIn("context.runtime_message_visible", status)
        self.assertIn("Restart paused", status)
        report_fault = compact(function_body(GUI_SOURCE, "gui_imgui_report_fault"))
        self.assertIn("runtime_message_visible = true", report_fault)

    def test_dimension_impossible_controls_are_omitted_with_safe_recovery(self) -> None:
        setup = compact(function_body(GUI_PANELS_SOURCE, "draw_setup"))
        setup_3d = setup[
            setup.index("const EnumOption modes_3d[]") :
            setup.index("const EnumOption modes_2d[]")
        ]
        self.assertNotIn("SIM_INIT_SQUARE_LATTICE", setup_3d)
        self.assertNotIn("SIM_INIT_MEMBRANE_OSMOSIS", setup_3d)
        self.assertIn("SIM_INIT_HEX_LATTICE", setup_3d)
        self.assertIn("box_modes_2d", setup)
        self.assertNotIn("box_modes_3d", setup)
        self.assertIn(
            "init.box_length_mode = SIM_BOX_FROM_PACKING_FRACTION", setup
        )

        boundaries = compact(
            function_body(GUI_PANELS_SOURCE, "draw_boundaries")
        )
        boundaries_2d = boundaries[
            boundaries.index("const EnumOption boundaries_2d[]") :
            boundaries.index("if (boundary_changed)")
        ]
        self.assertNotIn("SIM_BOUNDARY_WALL_XYZ", boundaries_2d)
        self.assertIn("SIM_BOUNDARY_CIRCULAR_SPHERICAL", boundaries_2d)
        self.assertIn("if (!is3d)", boundaries)
        self.assertIn("geometry.center_membrane = false", boundaries)

        appearance = compact(
            function_body(GUI_PANELS_SOURCE, "draw_appearance")
        )
        self.assertIn("if (is3d)", appearance)
        self.assertIn("sphere_quality_3d", appearance)
        color_modes_3d = appearance[
            appearance.index("const EnumOption color_modes_3d[]") :
            appearance.index("const EnumOption color_modes_2d[]")
        ]
        self.assertNotIn("SIM_COLOR_BY_HEX_ORDER", color_modes_3d)
        self.assertNotIn("SIM_COLOR_BY_SQUARE_ORDER", color_modes_3d)
        self.assertIn("live.color_mode = SIM_COLOR_BY_TYPE", appearance)

    def test_particle_drag_applies_a_rescheduled_force_without_teleporting(self) -> None:
        begin = compact(function_body(GRAPHICS_SOURCE, "beginParticleDrag"))
        self.assertIn("gui_imgui_simulation_faulted()", begin)
        self.assertIn("pointOnDragPlane(ray, planePoint, planeNormal, &hit)", begin)
        self.assertIn("lastWallTime", begin)
        self.assertNotIn("addWell && U != 0.0", begin)
        self.assertNotIn("graphicsPauseSimulation", begin)
        self.assertNotIn("removeFromCell", begin)
        self.assertNotIn("freeFly(p)", begin)

        get_input = compact(function_body(GRAPHICS_SOURCE, "getInput"))
        self.assertIn("pendingParticlePress.pending", get_input)
        moved = compact(
            function_body(GRAPHICS_SOURCE, "pendingParticlePressMoved")
        )
        self.assertIn("PARTICLE_DRAG_THRESHOLD_PIXELS", moved)
        self.assertIn("dx*dx + dy*dy > thresholdSquared", moved)
        start_pending = compact(
            function_body(GRAPHICS_SOURCE, "startPendingParticleDrag")
        )
        self.assert_in_order(
            start_pending,
            "beginParticleDrag(screenWindow, screenState, dragMouse)",
            "updateParticleDrag(screenWindow, dragMouse)",
            "if (finishImmediately)",
            "finishParticleDrag(screenState)",
        )
        self.assertRegex(
            get_input,
            r"IsMouseButtonReleased\(MOUSE_BUTTON_LEFT\).*?if \(moved\).*?"
            r"startPendingParticleDrag\(screenWindow, screenState, currentScreen, true\)",
        )
        self.assertRegex(
            get_input,
            r"IsMouseButtonDown\(MOUSE_BUTTON_LEFT\) && moved.*?"
            r"startPendingParticleDrag\(screenWindow, screenState, currentScreen, false\)",
        )
        press_branch = get_input[get_input.rfind(
            "else if (mouseInScene && IsMouseButtonPressed"
        ):]
        self.assertIn("pendingParticlePress = (particlePress)", press_branch)
        self.assertNotIn("beginParticleDrag", press_branch)

        hit_test = compact(
            function_body(GRAPHICS_SOURCE, "getWhatsUnderClick")
        )
        self.assertNotIn("freeFly(", hit_test)
        self.assertNotIn("graphicsPauseSimulation", hit_test)

        update = compact(function_body(GRAPHICS_SOURCE, "updateParticleDrag"))
        self.assert_in_order(
            update,
            "freeFly(p)",
            "const double oldVx = p->vx",
            "PBC(&dx, &dy",
            "const double response = 1.0 - pow",
            "p->vx += response*(dx - p->vx)",
            "p->coll++",
            "removeEventFromQueue(eventList[p->num])",
            "removeEventFromQueue(eventList[N + p->num])",
            "crossingEvent(p->num)",
            "collisionEvent(p->num)",
            "physicalQ()",
        )
        self.assertIn("beginTransportVelocityBatch()", update)
        self.assertIn("recordTransportVelocityChange(p, oldVx, oldVy, oldVz)", update)
        self.assertNotIn("p->x =", update)
        self.assertNotIn("p->y =", update)
        self.assertNotIn("p->z =", update)
        self.assertNotIn("removeFromCell", update)
        self.assertNotIn("addToCell", update)

        finish = compact(function_body(GRAPHICS_SOURCE, "finishParticleDrag"))
        self.assertIn("screenState->particleUnderClick = p", finish)
        self.assertNotIn("draggedParticlePlacementValid", finish)
        self.assertNotIn("live_analysis_reset_generation", finish)
        self.assertNotIn("crossingEvent", finish)

    def test_binary_mass_and_narrow_wall_initialization_are_stable(self) -> None:
        start = EDMD_SOURCE.index("void particlesInit(){")
        end = EDMD_SOURCE.index("void pairsInit(){", start)
        particles_init = compact(EDMD_SOURCE[start:end])
        self.assertIn("p->m = pow(sizeratio, 3)", particles_init)
        self.assertNotRegex(
            particles_init,
            r"p->m = pow\(sizeratio, 3\); p->m = 1;",
        )
        self.assertNotIn("drand(2.5, Lx - 2.5)", particles_init)
        self.assertNotIn("drand(2.5, Ly - 2.5)", particles_init)
        self.assertNotIn("drand(2.5, Lz - 2.5)", particles_init)
        self.assertIn("addWallx || addMidWall", particles_init)
        self.assertIn("p->y = addWally", particles_init)
        self.assertIn("p->z = addWallz", particles_init)

        bounded = compact(
            function_body(EDMD_SOURCE, "randomCoordinateWithMargin")
        )
        self.assertIn("largestSafeMargin", bounded)
        self.assertIn("maximum > minimum ? drand(minimum, maximum)", bounded)

    def test_fault_status_is_nonmodal_and_blocks_progress_controls(self) -> None:
        toolbar = compact(function_body(GUI_SOURCE, "draw_toolbar"))
        self.assertIn("BeginDisabled(context.simulation_faulted)", toolbar)
        self.assertIn("live.running || context.simulation_faulted", toolbar)
        status = compact(function_body(GUI_SOURCE, "draw_runtime_status"))
        self.assertNotIn("OpenPopup", status)
        self.assertIn("Restart paused", status)
        self.assertIn("assign_restart", status)

        inspector = compact(function_body(GUI_SOURCE, "draw_inspector"))
        self.assertIn(
            "particle_edit_disabled = live.running || context.simulation_faulted",
            inspector,
        )
        update = compact(
            function_body(GRAPHICS_SOURCE, "updateSelectedParticle")
        )
        self.assert_in_order(
            update,
            "gui_imgui_simulation_faulted()",
            "return;",
            "freeFly(p)",
        )

    def test_explicit_restart_preserves_run_pause_state(self) -> None:
        assign = compact(
            function_body(GUI_SOURCE, "assign_restart_command")
        )
        self.assertIn("restart.rendering.start_running = running", assign)

        toolbar = compact(function_body(GUI_SOURCE, "draw_toolbar"))
        self.assertIn(
            "assign_restart(command, context.applied, live.running)", toolbar
        )
        status = compact(function_body(GUI_SOURCE, "draw_runtime_status"))
        self.assertIn(
            "assign_restart(command, context.applied, live.running)", status
        )
        self.assertIn("Recovery stays paused", status)

        handler = compact(function_body(GRAPHICS_SOURCE, "handleGuiCommand"))
        self.assert_in_order(
            handler,
            "simulation_restart_core(&command->config",
            "screenState->running = command->config.rendering.start_running",
            "guiLive.running = screenState->running",
        )

    def test_make_test_artifacts_are_configuration_scoped(self) -> None:
        prediction_dir = re.search(
            r"^PREDICTION_TEST_BUILD_DIR\s*=\s*(.+)$", MAKEFILE, re.MULTILINE
        )
        self.assertIsNotNone(prediction_dir)
        self.assertIn("$(GRAPHICS_DIMENSION)", prediction_dir.group(1))
        self.assertIn("$(GRAPHICS_FEATURE_SET)", prediction_dir.group(1))
        self.assertIn("$(PREDICTION_TEST_CORE_OBJECT)", MAKEFILE)
        self.assertIn("$(PREDICTION_TEST_BINARY)", MAKEFILE)


if __name__ == "__main__":
    unittest.main(verbosity=2)
