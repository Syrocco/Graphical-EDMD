# Vendored GUI dependencies

The graphical build vendors two source dependencies so that it remains
reproducible without system-wide ImGui packages:

- **Dear ImGui 1.92.7 docking:** upstream tag `v1.92.7-docking`, commit
  `b1bcb12a624af7509894c8e77dd47416997777fa`. Dear ImGui is distributed under
  the MIT License; see [`imgui/LICENSE.txt`](imgui/LICENSE.txt).
- **rlImGui for raylib 5.5:** upstream release tag `Raylib_5_5`, release commit
  `118221c`. rlImGui is distributed under the zlib License; see
  [`rlimgui/LICENSE`](rlimgui/LICENSE).

These revisions are intentionally paired because rlImGui's renderer uses the
Dear ImGui 1.92 texture API. Update and test them together. Font Awesome is
disabled with `NO_FONT_AWESOME` because its optional generated assets are not
vendored.
