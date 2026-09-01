#ifndef EDMD_SIM_LIFECYCLE_H
#define EDMD_SIM_LIFECYCLE_H

#include "sim_config.h"

#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Replace the currently allocated graphical simulation with config.
 *
 * Validation and input-file preflight happen before the current simulation is
 * destroyed.  Physics/configuration changes are therefore transactional up to
 * the point at which the legacy allocation routines are entered.  Those
 * routines still terminate the process on an out-of-memory error.
 *
 * On failure, false is returned and error receives a user-facing explanation
 * when error_capacity is non-zero.  On success, error is set to an empty
 * string.  The renderer owns its caches and selection state; callers should
 * invalidate those around a successful restart.  This function is not thread
 * safe and must run on the simulation thread between event batches.
 */
bool simulation_restart_core(const SimulationConfig *config,
                             char *error,
                             size_t error_capacity);

/* Result of applying a validated graphical configuration to the simulation
 * which is already running.  The successful in-place results are ordered by
 * the amount of derived state which had to be refreshed; none of them replace
 * the particle array or reset simulation time/counters. */
typedef enum SimulationReconfigureAction {
    SIM_RECONFIGURE_FAILED = -1,
    SIM_RECONFIGURE_UNCHANGED = 0,
    SIM_RECONFIGURE_GLOBALS_ONLY,
    SIM_RECONFIGURE_EVENTS_REBUILT,
    SIM_RECONFIGURE_SPATIAL_REBUILT,
    SIM_RECONFIGURE_RESTART_REQUIRED
} SimulationReconfigureAction;

/* Apply next relative to the configuration which owns the current particle
 * state.  Validation and every geometry compatibility check run before any
 * mutation.  State-defining edits and valid edits which cannot preserve the
 * current trajectory return RESTART_REQUIRED and do not change the running
 * system, allowing the caller to rebuild it from next.  FAILED is reserved
 * for an invalid request/API use and likewise leaves the system untouched.
 *
 * Rendering fields are deliberately ignored here because their owner applies
 * them directly.  Call only on the simulation thread between event batches. */
SimulationReconfigureAction simulation_reconfigure_core(
    const SimulationConfig *current,
    const SimulationConfig *next,
    char *message,
    size_t message_capacity);

#ifdef __cplusplus
}
#endif

#endif
