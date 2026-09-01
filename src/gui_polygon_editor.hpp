#ifndef EDMD_GUI_POLYGON_EDITOR_HPP
#define EDMD_GUI_POLYGON_EDITOR_HPP

#include "sim_config.h"

namespace edmd_gui {

/* Draw the 2D convex-particle shape editor.

   Geometry manipulated with the mouse is kept in an editor-local staging
   buffer.  The supplied configuration is changed, and true is returned, only
   when the interaction is committed (normally on mouse release). */
bool draw_polygon_editor(SimParticleShapeConfig &shape);

} // namespace edmd_gui

#endif
