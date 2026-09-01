#include "sim_lifecycle.h"

#include "EDMD.h"

#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef G
#define G 0
#endif

#ifndef STRESS
#define STRESS 0
#endif

#ifndef HEAT
#define HEAT 0
#endif

#ifndef THREE_D
#define THREE_D 0
#endif

enum {
    RESTART_VALIDATION_ISSUE_CAPACITY = 128,
    TEMPERATURE_PROFILE_CAPACITY = 40,
    ENGINE_INITIAL_PATH_CAPACITY = 350,
    EVENT_LIST_AUXILIARY_SLOTS = 10
};

#define MAX_ENGINE_PARTICLE_COUNT \
    ((INT_MAX - EVENT_LIST_AUXILIARY_SLOTS) / 2)

static void restart_error(char *error, size_t error_capacity,
                          const char *format, ...)
{
    va_list arguments;

    if (error == NULL || error_capacity == 0) {
        return;
    }

    va_start(arguments, format);
    (void)vsnprintf(error, error_capacity, format, arguments);
    va_end(arguments);
    error[error_capacity - 1] = '\0';
}

static bool restart_config_is_valid(const SimulationConfig *config,
                                    char *error, size_t error_capacity)
{
    BuildCapabilities capabilities;
    ConfigIssue issues[RESTART_VALIDATION_ISSUE_CAPACITY];
    size_t issue_count;
    size_t stored_count;
    size_t error_count = 0;
    const ConfigIssue *first_error = NULL;

    capabilities = simulation_build_capabilities();
    issue_count = simulation_config_validate(
        config, &capabilities, issues, RESTART_VALIDATION_ISSUE_CAPACITY);
    stored_count = issue_count;
    if (stored_count > RESTART_VALIDATION_ISSUE_CAPACITY) {
        stored_count = RESTART_VALIDATION_ISSUE_CAPACITY;
    }

    for (size_t index = 0; index < stored_count; ++index) {
        if (issues[index].severity == CONFIG_ISSUE_ERROR) {
            if (first_error == NULL) {
                first_error = &issues[index];
            }
            ++error_count;
        }
    }

    if (first_error != NULL) {
        if (error_count == 1) {
            restart_error(error, error_capacity, "%s: %s",
                          first_error->field, first_error->message);
        } else {
            restart_error(error, error_capacity, "%s: %s (%zu errors total)",
                          first_error->field, first_error->message,
                          error_count);
        }
        return false;
    }

    /* A validator should keep its error set within the generous stack buffer.
     * Refuse an unexpectedly truncated result instead of applying an unchecked
     * tail of issues after the old simulation has been freed. */
    if (issue_count > RESTART_VALIDATION_ISSUE_CAPACITY) {
        restart_error(error, error_capacity,
                      "Configuration validation produced too many issues");
        return false;
    }

    return true;
}

static bool path_is_terminated(const char path[SIM_CONFIG_PATH_CAPACITY])
{
    return memchr(path, '\0', SIM_CONFIG_PATH_CAPACITY) != NULL;
}

static bool token_is_int(const char *token, int *value)
{
    char *end = NULL;
    long parsed;

    if (token == NULL || token[0] == '\0') return false;
    errno = 0;
    parsed = strtol(token, &end, 10);
    if (errno == ERANGE || parsed < INT_MIN || parsed > INT_MAX ||
        end == token || *end != '\0') {
        return false;
    }
    if (value != NULL) *value = (int)parsed;
    return true;
}

static bool token_is_double(const char *token, double *value)
{
    char *end = NULL;
    double parsed;

    if (token == NULL || token[0] == '\0') return false;
    errno = 0;
    parsed = strtod(token, &end);
    if (errno == ERANGE || end == token || *end != '\0' ||
        !isfinite(parsed)) {
        return false;
    }
    if (value != NULL) *value = parsed;
    return true;
}

static int split_tokens(char *line, char **tokens, int capacity)
{
    char *save = NULL;
    char *token = strtok_r(line, " \t\r\n", &save);
    int count = 0;

    while (token != NULL) {
        if (count < capacity) tokens[count] = token;
        ++count;
        token = strtok_r(NULL, " \t\r\n", &save);
    }
    return count;
}

static bool line_was_complete(const char *line, FILE *input)
{
    return strchr(line, '\n') != NULL || feof(input);
}

static bool particle_count_fits_engine(int particle_count)
{
    return particle_count > 0 &&
           particle_count <= MAX_ENGINE_PARTICLE_COUNT;
}

static bool config_has_x_walls(const SimulationConfig *config)
{
    return config->geometry.boundary == SIM_BOUNDARY_WALL_X ||
           config->geometry.boundary == SIM_BOUNDARY_WALL_XY ||
           config->geometry.boundary == SIM_BOUNDARY_WALL_XYZ;
}

static bool config_has_y_walls(const SimulationConfig *config)
{
    return config->geometry.boundary == SIM_BOUNDARY_WALL_Y ||
           config->geometry.boundary == SIM_BOUNDARY_WALL_XY ||
           config->geometry.boundary == SIM_BOUNDARY_WALL_XYZ;
}

#if THREE_D
static bool config_has_z_walls(const SimulationConfig *config)
{
    return config->geometry.boundary == SIM_BOUNDARY_WALL_XYZ;
}
#endif

typedef enum DumpPreflightSection {
    DUMP_SECTION_NONE,
    DUMP_SECTION_TIMESTEP,
    DUMP_SECTION_ATOM_COUNT,
    DUMP_SECTION_BOX,
    DUMP_SECTION_ATOMS,
    DUMP_SECTION_OTHER
} DumpPreflightSection;

typedef struct DumpFramePreflight {
    bool present;
    bool coherent;
    bool timestep_header;
    bool timestep_value;
    bool atom_count_header;
    bool atom_count_value;
    bool box_header;
    bool atoms_header;
    bool atom_properties;
    bool atom_rows_valid;
    int atom_count;
    int box_rows;
    double box_lower[3];
    double box_upper[3];
    double minimum_particle_extent[3];
    double maximum_particle_extent[3];
    size_t atom_rows;
    int atom_columns;
    int type_column;
    int x_column;
    int y_column;
    int z_column;
    int vx_column;
    int vy_column;
    int vz_column;
    int radius_column;
    int mass_column;
    DumpPreflightSection section;
} DumpFramePreflight;

static void reset_dump_frame(DumpFramePreflight *frame)
{
    memset(frame, 0, sizeof(*frame));
    frame->present = true;
    frame->coherent = true;
    frame->atom_rows_valid = true;
    frame->type_column = -1;
    frame->x_column = -1;
    frame->y_column = -1;
    frame->z_column = -1;
    frame->vx_column = -1;
    frame->vy_column = -1;
    frame->vz_column = -1;
    frame->radius_column = -1;
    frame->mass_column = -1;
    frame->minimum_particle_extent[0] = INFINITY;
    frame->minimum_particle_extent[1] = INFINITY;
    frame->minimum_particle_extent[2] = INFINITY;
    frame->maximum_particle_extent[0] = -INFINITY;
    frame->maximum_particle_extent[1] = -INFINITY;
    frame->maximum_particle_extent[2] = -INFINITY;
    frame->section = DUMP_SECTION_TIMESTEP;
    frame->timestep_header = true;
}

static bool dump_item_is(const char *line, const char *name)
{
    const size_t length = strlen(name);
    return strncmp(line, name, length) == 0 &&
           (line[length] == '\0' || line[length] == '\r' ||
            line[length] == '\n' || isspace((unsigned char)line[length]));
}

static bool parse_single_dump_double(char *line)
{
    char *tokens[2];
    double value = 0.0;
    return split_tokens(line, tokens, 2) == 1 &&
           token_is_double(tokens[0], &value);
}

static bool parse_single_dump_int(char *line, int *value)
{
    char *tokens[2];
    return split_tokens(line, tokens, 2) == 1 &&
           token_is_int(tokens[0], value);
}

static bool parse_dump_bounds(char *line, double *lower_result,
                              double *upper_result)
{
    char *tokens[3];
    double lower = 0.0;
    double upper = 0.0;
    if (split_tokens(line, tokens, 3) != 2 ||
        !token_is_double(tokens[0], &lower) ||
        !token_is_double(tokens[1], &upper) || upper <= lower) {
        return false;
    }
    if (lower_result != NULL) *lower_result = lower;
    if (upper_result != NULL) *upper_result = upper;
    return true;
}

static void parse_dump_atom_header(char *line, DumpFramePreflight *frame)
{
    char *save = NULL;
    char *token = strtok_r(line, " \t\r\n", &save);
    int column = -2;

    while (token != NULL) {
        if (column >= 0) {
            if (strcmp(token, "type") == 0) frame->type_column = column;
            else if (strcmp(token, "x") == 0) frame->x_column = column;
            else if (strcmp(token, "y") == 0) frame->y_column = column;
            else if (strcmp(token, "z") == 0) frame->z_column = column;
            else if (strcmp(token, "vx") == 0) frame->vx_column = column;
            else if (strcmp(token, "vy") == 0) frame->vy_column = column;
            else if (strcmp(token, "vz") == 0) frame->vz_column = column;
            else if (strcmp(token, "radius") == 0) frame->radius_column = column;
            else if (strcmp(token, "m") == 0) frame->mass_column = column;
        }
        ++column;
        token = strtok_r(NULL, " \t\r\n", &save);
    }

    frame->atom_columns = column;
    frame->atom_properties = frame->type_column >= 0 &&
        frame->x_column >= 0 && frame->y_column >= 0 &&
        frame->vx_column >= 0 && frame->vy_column >= 0 &&
        frame->radius_column >= 0 && frame->mass_column >= 0
#if THREE_D
        && frame->z_column >= 0 && frame->vz_column >= 0
#endif
        ;
}

static bool parse_dump_atom_row(char *line,
                                const DumpFramePreflight *frame,
                                double *x_result, double *y_result,
                                double *z_result,
                                double *radius_result)
{
    char *save = NULL;
    char *token = strtok_r(line, " \t\r\n", &save);
    int column = 0;
    int parsed_type = 0;
    bool type_ok = false;
    bool x_ok = false;
    bool y_ok = false;
    bool z_ok = THREE_D == 0;
    bool vx_ok = false;
    bool vy_ok = false;
    bool vz_ok = THREE_D == 0;
    bool radius_ok = false;
    bool mass_ok = false;
    double parsed_x = 0.0;
    double parsed_y = 0.0;
    double parsed_z = 0.0;
    double parsed_radius = 0.0;

    while (token != NULL) {
        double value = 0.0;
        if (column == frame->type_column) {
            type_ok = token_is_int(token, &parsed_type);
        } else if (column == frame->x_column) {
            x_ok = token_is_double(token, &parsed_x);
        } else if (column == frame->y_column) {
            y_ok = token_is_double(token, &parsed_y);
        } else if (column == frame->z_column) {
            z_ok = token_is_double(token, &parsed_z);
        } else if (column == frame->vx_column) {
            vx_ok = token_is_double(token, &value);
        } else if (column == frame->vy_column) {
            vy_ok = token_is_double(token, &value);
        } else if (column == frame->vz_column) {
            vz_ok = token_is_double(token, &value);
        } else if (column == frame->radius_column) {
            radius_ok = token_is_double(token, &parsed_radius) &&
                        parsed_radius > 0.0;
        } else if (column == frame->mass_column) {
            mass_ok = token_is_double(token, &value) && value > 0.0;
        }
        ++column;
        token = strtok_r(NULL, " \t\r\n", &save);
    }

    const bool valid = column >= frame->atom_columns && type_ok && x_ok &&
        y_ok && z_ok && vx_ok && vy_ok && vz_ok && radius_ok && mass_ok;
    if (valid) {
        if (x_result != NULL) *x_result = parsed_x;
        if (y_result != NULL) *y_result = parsed_y;
        if (z_result != NULL) *z_result = parsed_z;
        if (radius_result != NULL) *radius_result = parsed_radius;
    }
    return valid;
}

static bool preflight_dump_configuration(FILE *input, const char *path,
                                         const SimulationConfig *config,
                                         char *error, size_t error_capacity)
{
    char line[1024];
    DumpFramePreflight frame;

    memset(&frame, 0, sizeof(frame));

    /* The legacy parser exits the process when a required dump section is
     * absent.  Track only the final frame because particlesInit() jumps to it;
     * sections from earlier frames must not accidentally satisfy its checks. */
    while (fgets(line, sizeof(line), input) != NULL) {
        if (!line_was_complete(line, input)) {
            restart_error(error, error_capacity,
                          "Dump file '%s' contains a line longer than %zu bytes",
                          path, sizeof(line) - 1);
            return false;
        }

        if (dump_item_is(line, "ITEM: TIMESTEP")) {
            reset_dump_frame(&frame);
            continue;
        }
        if (!frame.present) continue;

        if (dump_item_is(line, "ITEM: NUMBER OF ATOMS")) {
            if (frame.atom_count_header) frame.coherent = false;
            frame.atom_count_header = true;
            frame.section = DUMP_SECTION_ATOM_COUNT;
            continue;
        }
        if (dump_item_is(line, "ITEM: BOX BOUNDS")) {
            if (frame.box_header) frame.coherent = false;
            frame.box_header = true;
            frame.section = DUMP_SECTION_BOX;
            continue;
        }
        if (dump_item_is(line, "ITEM: ATOMS")) {
            if (frame.atoms_header) frame.coherent = false;
            frame.atoms_header = true;
            frame.section = DUMP_SECTION_ATOMS;
            parse_dump_atom_header(line, &frame);
            continue;
        }
        if (strncmp(line, "ITEM:", 5) == 0) {
            frame.section = DUMP_SECTION_OTHER;
            continue;
        }

        switch (frame.section) {
            case DUMP_SECTION_TIMESTEP:
                if (frame.timestep_value || !parse_single_dump_double(line)) {
                    frame.coherent = false;
                } else {
                    frame.timestep_value = true;
                }
                break;
            case DUMP_SECTION_ATOM_COUNT:
                if (frame.atom_count_value ||
                    !parse_single_dump_int(line, &frame.atom_count)) {
                    frame.coherent = false;
                } else {
                    frame.atom_count_value = true;
                }
                break;
            case DUMP_SECTION_BOX:
                if (line[0] != '\n' && line[0] != '\r') {
                    /* 2D LAMMPS dumps conventionally retain a third dummy
                     * z row (often "0 0"), which the engine does not read. */
                    if (frame.box_rows < (THREE_D ? 3 : 2)) {
                        if (!parse_dump_bounds(
                                line, &frame.box_lower[frame.box_rows],
                                &frame.box_upper[frame.box_rows])) {
                            frame.coherent = false;
                        }
                    }
                    ++frame.box_rows;
                }
                break;
            case DUMP_SECTION_ATOMS:
                if (line[0] != '\n' && line[0] != '\r') {
                    if (frame.atom_properties) {
                        double x = 0.0;
                        double y = 0.0;
                        double z = 0.0;
                        double radius = 0.0;
                        if (!parse_dump_atom_row(line, &frame, &x, &y,
                                                 &z,
                                                 &radius)) {
                            frame.atom_rows_valid = false;
                        } else {
                            frame.minimum_particle_extent[0] = fmin(
                                frame.minimum_particle_extent[0], x - radius);
                            frame.maximum_particle_extent[0] = fmax(
                                frame.maximum_particle_extent[0], x + radius);
                            frame.minimum_particle_extent[1] = fmin(
                                frame.minimum_particle_extent[1], y - radius);
                            frame.maximum_particle_extent[1] = fmax(
                                frame.maximum_particle_extent[1], y + radius);
#if THREE_D
                            frame.minimum_particle_extent[2] = fmin(
                                frame.minimum_particle_extent[2], z - radius);
                            frame.maximum_particle_extent[2] = fmax(
                                frame.maximum_particle_extent[2], z + radius);
#endif
                        }
                    }
                    ++frame.atom_rows;
                }
                break;
            case DUMP_SECTION_NONE:
            case DUMP_SECTION_OTHER:
                break;
            }
    }

    if (ferror(input)) {
        restart_error(error, error_capacity,
                      "Could not read dump file '%s': %s",
                      path, strerror(errno));
        return false;
    }
    if (!frame.present) {
        restart_error(error, error_capacity,
                      "Dump file '%s' contains no timestep frame", path);
        return false;
    }
    if (!frame.coherent || !frame.timestep_header || !frame.timestep_value ||
        !frame.atom_count_header || !frame.atom_count_value ||
        !frame.box_header || !frame.atoms_header) {
        restart_error(error, error_capacity,
                      "The final frame of dump file '%s' has missing, duplicate, or malformed required sections",
                      path);
        return false;
    }
    if (!particle_count_fits_engine(frame.atom_count)) {
        restart_error(error, error_capacity,
                      "The final dump frame requests %d particles; the signed 2*N event index supports at most %d",
                      frame.atom_count, MAX_ENGINE_PARTICLE_COUNT);
        return false;
    }
    if (frame.box_rows < (THREE_D ? 3 : 2)) {
        restart_error(error, error_capacity,
                      "The final frame of dump file '%s' has only %d box-bound rows; this %dD build requires %d",
                      path, frame.box_rows, THREE_D ? 3 : 2,
                      THREE_D ? 3 : 2);
        return false;
    }
    for (int axis = 0; axis < (THREE_D ? 3 : 2); ++axis) {
        if (frame.box_lower[axis] != 0.0) {
            const char axis_name = axis == 0 ? 'x' : (axis == 1 ? 'y' : 'z');
            restart_error(error, error_capacity,
                          "Dump file '%s' uses unsupported nonzero %c-box origin %.17g; the current loader requires a zero lower bound because it does not translate particle coordinates",
                          path, axis_name, frame.box_lower[axis]);
            return false;
        }
    }
    if (!frame.atom_properties) {
        restart_error(error, error_capacity,
                      "The final frame of dump file '%s' does not provide every particle column required by this %dD build",
                      path, THREE_D ? 3 : 2);
        return false;
    }
    if (frame.atom_rows < (size_t)frame.atom_count) {
        restart_error(error, error_capacity,
                      "The final frame of dump file '%s' declares %d particles but provides only %zu atom rows",
                      path, frame.atom_count, frame.atom_rows);
        return false;
    }
    if (!frame.atom_rows_valid) {
        restart_error(error, error_capacity,
                      "The final frame of dump file '%s' has a malformed atom row or a non-finite/non-positive required value",
                      path);
        return false;
    }
    if (config_has_x_walls(config) &&
        (frame.minimum_particle_extent[0] < frame.box_lower[0] ||
         frame.maximum_particle_extent[0] > frame.box_upper[0])) {
        restart_error(error, error_capacity,
                      "The final dump frame in '%s' places a particle outside the selected X walls (particle extent [%.17g, %.17g], box [%.17g, %.17g])",
                      path, frame.minimum_particle_extent[0],
                      frame.maximum_particle_extent[0], frame.box_lower[0],
                      frame.box_upper[0]);
        return false;
    }
    if (config_has_y_walls(config) &&
        (frame.minimum_particle_extent[1] < frame.box_lower[1] ||
         frame.maximum_particle_extent[1] > frame.box_upper[1])) {
        restart_error(error, error_capacity,
                      "The final dump frame in '%s' places a particle outside the selected Y walls (particle extent [%.17g, %.17g], box [%.17g, %.17g])",
                      path, frame.minimum_particle_extent[1],
                      frame.maximum_particle_extent[1], frame.box_lower[1],
                      frame.box_upper[1]);
        return false;
    }
#if THREE_D
    if (config_has_z_walls(config) &&
        (frame.minimum_particle_extent[2] < frame.box_lower[2] ||
         frame.maximum_particle_extent[2] > frame.box_upper[2])) {
        restart_error(error, error_capacity,
                      "The final dump frame in '%s' places a particle outside the selected Z walls (particle extent [%.17g, %.17g], box [%.17g, %.17g])",
                      path, frame.minimum_particle_extent[2],
                      frame.maximum_particle_extent[2], frame.box_lower[2],
                      frame.box_upper[2]);
        return false;
    }
#endif
    return true;
}

static bool preflight_initial_state(const SimulationConfig *config,
                                    char *error, size_t error_capacity)
{
    const SimInitMode mode = config->initialization.mode;
    const char *path = config->initialization.load_path;
    FILE *input;

    if (mode != SIM_INIT_DUMP_FILE) {
        return true;
    }

    if (!path_is_terminated(config->initialization.load_path)) {
        restart_error(error, error_capacity,
                      "The initial-state path is not null-terminated");
        return false;
    }
    if (path[0] == '\0') {
        restart_error(error, error_capacity,
                      "Choose an initial-state file before restarting");
        return false;
    }
    if (strlen(path) >= ENGINE_INITIAL_PATH_CAPACITY) {
        restart_error(error, error_capacity,
                      "Initial-state path is too long (maximum %d characters)",
                      ENGINE_INITIAL_PATH_CAPACITY - 1);
        return false;
    }

    errno = 0;
    input = fopen(path, "r");
    if (input == NULL) {
        restart_error(error, error_capacity,
                      "Unable to open initial-state file '%s': %s",
                      path, strerror(errno));
        return false;
    }

    const bool valid = preflight_dump_configuration(
        input, path, config, error, error_capacity);
    (void)fclose(input);
    return valid;
}

#if G == 1

/* Mutable engine state.  These declarations intentionally live here instead
 * of broadening EDMD.h: the lifecycle bridge is the only new C module which
 * should manipulate restart internals directly. */
extern int N;
extern int Nbig;
extern int Nsmall;
extern bool addWell;

extern double t;
extern double U;
extern double Lx;
extern double Ly;
#if THREE_D
extern double Lz;
#endif
extern unsigned long int ncol;
extern unsigned long int ncross;
extern double paulTime;
extern int actualPaulList;
extern double t1;
extern double nextScreen;
extern double dtimeThermo;
extern int logSpacing;

extern double px;
extern double py;
extern double E;
extern double pressure;
extern double dxr;
extern double dyr;
extern double dist;
extern double active;
extern double Ep;
extern double MSD;
extern double COM_x;
extern double COM_y;
extern double lastScreen;
extern double lastCollNum;

#if THREE_D
extern double pz;
extern double dzr;
extern double COM_z;
#endif

extern double collTermX;
extern double collTermY;
extern double collTermXY;
extern double collTermYX;
#if THREE_D
extern double collTermZ;
extern double collTermXZ;
extern double collTermZX;
extern double collTermYZ;
extern double collTermZY;
#endif

extern double dp;
extern double dpLeft;
extern double dpRight;
extern double dpMid;
extern double JxProfile[TEMPERATURE_PROFILE_CAPACITY];
extern double JyProfile[TEMPERATURE_PROFILE_CAPACITY];
#if THREE_D
extern double JzProfile[TEMPERATURE_PROFILE_CAPACITY];
#endif
extern double energyProfile[TEMPERATURE_PROFILE_CAPACITY];
extern double countTempProfile[TEMPERATURE_PROFILE_CAPACITY];

extern particle *particles;
extern particle **cellList;
extern node **eventList;
extern node *root;
extern node *nextEvent;
extern void (*freeFly)(particle *);
extern FILE *file;
extern char initialFilePath[ENGINE_INITIAL_PATH_CAPACITY];

void init_genrand(unsigned long seed);

static void free_running_simulation(void)
{
    outputLifecycleClose();
    freeArrays();
    root = NULL;
    nextEvent = NULL;
    file = NULL;
}

static void reset_runtime_accumulators(void)
{
    t = 0.0;
    ncol = 0;
    ncross = 0;
    paulTime = 0.0;
    actualPaulList = 0;
    t1 = 0.0;

    Nbig = 0;
    Nsmall = 0;

    px = 0.0;
    py = 0.0;
    E = 0.0;
    pressure = 0.0;
    dxr = 0.0;
    dyr = 0.0;
    dist = 0.0;
    active = 0.0;
    Ep = 0.0;
    MSD = 0.0;
    COM_x = 0.0;
    COM_y = 0.0;
    lastScreen = 0.0;
    lastCollNum = 0.0;

#if THREE_D
    pz = 0.0;
    dzr = 0.0;
    COM_z = 0.0;
#endif

    collTermX = 0.0;
    collTermY = 0.0;
    collTermXY = 0.0;
    collTermYX = 0.0;
#if THREE_D
    collTermZ = 0.0;
    collTermXZ = 0.0;
    collTermZX = 0.0;
    collTermYZ = 0.0;
    collTermZY = 0.0;
#endif

    dp = 0.0;
    cumulativeWallImpulse = 0.0;
    dpLeft = 0.0;
    dpRight = 0.0;
    dpMid = 0.0;

    (void)memset(JxProfile, 0, sizeof(double) * TEMPERATURE_PROFILE_CAPACITY);
    (void)memset(JyProfile, 0, sizeof(double) * TEMPERATURE_PROFILE_CAPACITY);
#if THREE_D
    (void)memset(JzProfile, 0, sizeof(double) * TEMPERATURE_PROFILE_CAPACITY);
#endif
    (void)memset(energyProfile, 0,
                 sizeof(double) * TEMPERATURE_PROFILE_CAPACITY);
    (void)memset(countTempProfile, 0,
                 sizeof(double) * TEMPERATURE_PROFILE_CAPACITY);
}

static void refresh_population_counts(void)
{
    Nbig = 0;
    for (int index = 0; index < N; ++index) {
        Nbig += particles[index].type;
    }
    Nsmall = N - Nbig;
}

enum ReconfigurePlanBit {
    RECONFIGURE_PLAN_GLOBALS = 1u << 0,
    RECONFIGURE_PLAN_SYNCHRONIZE = 1u << 1,
    RECONFIGURE_PLAN_PARTICLE_EVENTS = 1u << 2,
    RECONFIGURE_PLAN_SPATIAL_INDEX = 1u << 3,
    RECONFIGURE_PLAN_THERMOSTAT_EVENT = 1u << 4,
    RECONFIGURE_PLAN_PHYSICAL_TOTALS = 1u << 5,
    RECONFIGURE_PLAN_SHEAR_REBASE = 1u << 6,
    RECONFIGURE_PLAN_CHARGES = 1u << 7,
    RECONFIGURE_PLAN_HEX_TYPES = 1u << 8,
    RECONFIGURE_PLAN_POTENTIAL_MEMBERSHIPS = 1u << 9,
    RECONFIGURE_PLAN_FULL_RESTART = 1u << 10
};

typedef unsigned int ReconfigurePlan;

static bool changed_double(double left, double right)
{
    return left != right;
}

static bool initialization_fields_changed(
    const SimInitializationConfig *left,
    const SimInitializationConfig *right)
{
    return left->mode != right->mode ||
        strcmp(left->load_path, right->load_path) != 0 ||
        left->particle_count != right->particle_count ||
        changed_double(left->packing_fraction, right->packing_fraction) ||
        left->composition != right->composition ||
        changed_double(left->size_ratio, right->size_ratio) ||
        changed_double(left->small_fraction, right->small_fraction) ||
        changed_double(left->polydisperse_a, right->polydisperse_a) ||
        changed_double(left->polydisperse_b, right->polydisperse_b) ||
        changed_double(left->initial_energy, right->initial_energy) ||
        left->random_seed != right->random_seed ||
        left->box_length_mode != right->box_length_mode ||
        changed_double(left->aspect_ratio, right->aspect_ratio) ||
        changed_double(left->box_lx, right->box_lx) ||
        changed_double(left->box_ly, right->box_ly) ||
        changed_double(left->box_lz, right->box_lz) ||
        changed_double(left->liquid_packing_fraction,
                       right->liquid_packing_fraction) ||
        changed_double(left->gas_packing_fraction,
                       right->gas_packing_fraction) ||
        changed_double(left->hex_hole_radius, right->hex_hole_radius) ||
        left->automatic_growth_rate != right->automatic_growth_rate ||
        changed_double(left->growth_rate, right->growth_rate);
}

static bool geometry_response_changed(const SimGeometryConfig *left,
                                      const SimGeometryConfig *right)
{
    return changed_double(left->wall_normal_restitution,
                          right->wall_normal_restitution) ||
        changed_double(left->wall_tangent_restitution,
                       right->wall_tangent_restitution) ||
        left->thermal_wall != right->thermal_wall ||
        left->shear_wall != right->shear_wall ||
        changed_double(left->shear_wall_speed, right->shear_wall_speed);
}

static bool dynamics_response_changed(const SimDynamicsConfig *left,
                                      const SimDynamicsConfig *right)
{
    return left->restitution_model != right->restitution_model ||
        changed_double(left->particle_restitution,
                       right->particle_restitution) ||
        changed_double(left->exponential_velocity_scale,
                       right->exponential_velocity_scale) ||
        changed_double(left->exponential_amplitude,
                       right->exponential_amplitude) ||
        left->injection != right->injection ||
        changed_double(left->delta, right->delta) ||
        changed_double(left->delta_small, right->delta_small) ||
        changed_double(left->delta_large, right->delta_large) ||
        changed_double(left->collision_memory_time,
                       right->collision_memory_time) ||
        changed_double(left->evolving_delta_tau,
                       right->evolving_delta_tau) ||
        changed_double(left->evolving_delta_max,
                       right->evolving_delta_max) ||
        changed_double(left->injected_energy, right->injected_energy) ||
        changed_double(left->background_injected_energy,
                       right->background_injected_energy) ||
        changed_double(left->injection_recovery_time,
                       right->injection_recovery_time) ||
        changed_double(left->injection_exponent,
                       right->injection_exponent) ||
        left->randomize_injected_energy !=
            right->randomize_injected_energy ||
        left->simplified_energy_injection !=
            right->simplified_energy_injection ||
        left->rough_tangential_delta != right->rough_tangential_delta ||
        changed_double(left->rough_tangential_restitution,
                       right->rough_tangential_restitution) ||
        changed_double(left->rough_angular_kick,
                       right->rough_angular_kick);
}

static bool driving_response_changed(const SimDrivingConfig *left,
                                     const SimDrivingConfig *right)
{
    return changed_double(left->temperature, right->temperature) ||
        changed_double(left->gamma, right->gamma) ||
        left->euler_noise != right->euler_noise;
}

static bool potential_is_effective(const SimPotentialConfig *potential)
{
    return potential->enabled && potential->energy != 0.0;
}

static bool trajectory_changed(const SimDynamicsConfig *left,
                               const SimDynamicsConfig *right)
{
    if (left->flight != right->flight) return true;
    if ((left->flight == SIM_FLIGHT_CONSTANT_FIELD ||
         right->flight == SIM_FLIGHT_CONSTANT_FIELD) &&
        changed_double(left->constant_field, right->constant_field)) {
        return true;
    }
    if ((left->flight == SIM_FLIGHT_LINEAR_DAMPING ||
         right->flight == SIM_FLIGHT_LINEAR_DAMPING) &&
        changed_double(left->damping_gamma, right->damping_gamma)) {
        return true;
    }
    return (left->flight == SIM_FLIGHT_LEES_EDWARDS_SHEAR ||
            right->flight == SIM_FLIGHT_LEES_EDWARDS_SHEAR) &&
        changed_double(left->shear_strain_rate,
                       right->shear_strain_rate);
}

static bool dynamics_known_fields_changed(const SimDynamicsConfig *left,
                                          const SimDynamicsConfig *right)
{
    return trajectory_changed(left, right) ||
        /* Inactive trajectory parameters are still remembered live. */
        changed_double(left->constant_field, right->constant_field) ||
        changed_double(left->damping_gamma, right->damping_gamma) ||
        changed_double(left->shear_strain_rate,
                       right->shear_strain_rate) ||
        dynamics_response_changed(left, right);
}

static ReconfigurePlan classify_reconfiguration(
    const SimulationConfig *current,
    const SimulationConfig *next)
{
    ReconfigurePlan plan = 0;
    const bool old_potential = potential_is_effective(&current->potential);
    const bool new_potential = potential_is_effective(&next->potential);
    const bool species_changed = current->potential.species_mode !=
        next->potential.species_mode;
    const bool charge_fractions_changed =
        changed_double(current->potential.positive_charge_fraction,
                       next->potential.positive_charge_fraction) ||
        changed_double(current->potential.neutral_charge_fraction,
                       next->potential.neutral_charge_fraction);

    if (current->schema_version != next->schema_version ||
        initialization_fields_changed(&current->initialization,
                                      &next->initialization)) {
        return RECONFIGURE_PLAN_FULL_RESTART;
    }

    if (current->geometry.center_membrane !=
        next->geometry.center_membrane) {
        plan |= RECONFIGURE_PLAN_GLOBALS |
                RECONFIGURE_PLAN_SYNCHRONIZE |
                RECONFIGURE_PLAN_PARTICLE_EVENTS;
    }
    if (current->geometry.boundary != next->geometry.boundary) {
        plan |= RECONFIGURE_PLAN_GLOBALS |
                RECONFIGURE_PLAN_SYNCHRONIZE |
                RECONFIGURE_PLAN_PARTICLE_EVENTS;
        /* Changing periodic images changes existing well memberships, but the
         * linked-cell dimensions and particle-to-cell ownership stay valid. */
        if (old_potential || new_potential)
            plan |= RECONFIGURE_PLAN_POTENTIAL_MEMBERSHIPS;
    }
    if (geometry_response_changed(&current->geometry, &next->geometry))
        plan |= RECONFIGURE_PLAN_GLOBALS;

    if (trajectory_changed(&current->dynamics, &next->dynamics)) {
        plan |= RECONFIGURE_PLAN_GLOBALS |
                RECONFIGURE_PLAN_SYNCHRONIZE |
                RECONFIGURE_PLAN_PARTICLE_EVENTS |
                RECONFIGURE_PLAN_PHYSICAL_TOTALS;
        if (current->dynamics.flight == SIM_FLIGHT_LEES_EDWARDS_SHEAR ||
            next->dynamics.flight == SIM_FLIGHT_LEES_EDWARDS_SHEAR) {
            plan |= RECONFIGURE_PLAN_SHEAR_REBASE;
        }
    }
    if (dynamics_response_changed(&current->dynamics, &next->dynamics) ||
        changed_double(current->dynamics.constant_field,
                       next->dynamics.constant_field) ||
        changed_double(current->dynamics.damping_gamma,
                       next->dynamics.damping_gamma) ||
        changed_double(current->dynamics.shear_strain_rate,
                       next->dynamics.shear_strain_rate)) {
        plan |= RECONFIGURE_PLAN_GLOBALS;
    }

    if (current->driving.thermostat != next->driving.thermostat ||
        ((current->driving.thermostat != SIM_THERMOSTAT_NONE ||
          next->driving.thermostat != SIM_THERMOSTAT_NONE) &&
         changed_double(current->driving.thermostat_interval,
                        next->driving.thermostat_interval))) {
        plan |= RECONFIGURE_PLAN_GLOBALS |
                RECONFIGURE_PLAN_THERMOSTAT_EVENT;
    }
    if (changed_double(current->driving.thermostat_interval,
                       next->driving.thermostat_interval))
        plan |= RECONFIGURE_PLAN_GLOBALS;
    if (driving_response_changed(&current->driving, &next->driving))
        plan |= RECONFIGURE_PLAN_GLOBALS;

    if (old_potential != new_potential ||
        ((old_potential || new_potential) &&
         changed_double(current->potential.range,
                        next->potential.range))) {
        plan |= RECONFIGURE_PLAN_GLOBALS |
                RECONFIGURE_PLAN_SYNCHRONIZE |
                RECONFIGURE_PLAN_PARTICLE_EVENTS |
                RECONFIGURE_PLAN_SPATIAL_INDEX;
    }
    if (current->potential.enabled != next->potential.enabled ||
        changed_double(current->potential.energy,
                       next->potential.energy) ||
        changed_double(current->potential.range,
                       next->potential.range))
        plan |= RECONFIGURE_PLAN_GLOBALS;
    if (species_changed) {
        plan |= RECONFIGURE_PLAN_GLOBALS;
        if (old_potential || new_potential) {
            plan |= RECONFIGURE_PLAN_SYNCHRONIZE |
                    RECONFIGURE_PLAN_PARTICLE_EVENTS |
                    RECONFIGURE_PLAN_POTENTIAL_MEMBERSHIPS;
        }
    }
    if (charge_fractions_changed) {
        plan |= RECONFIGURE_PLAN_GLOBALS;
    }
    /* Charge labels are runtime metadata only while a charged potential is
     * effective.  Remember inactive fraction edits without touching particles,
     * then realize them when that model is enabled. */
    if (new_potential &&
        next->potential.species_mode == SIM_POTENTIAL_CHARGED &&
        (!old_potential ||
         current->potential.species_mode != SIM_POTENTIAL_CHARGED ||
         charge_fractions_changed)) {
        plan |= RECONFIGURE_PLAN_CHARGES;
        if (old_potential && new_potential) {
            plan |= RECONFIGURE_PLAN_SYNCHRONIZE |
                    RECONFIGURE_PLAN_PARTICLE_EVENTS |
                    RECONFIGURE_PLAN_POTENTIAL_MEMBERSHIPS;
        }
    }
    if (current->initialization.mode == SIM_INIT_HEX_LATTICE) {
        const bool old_split = current->potential.enabled &&
            current->potential.species_mode == SIM_POTENTIAL_SAME_SPECIES;
        const bool new_split = next->potential.enabled &&
            next->potential.species_mode == SIM_POTENTIAL_SAME_SPECIES;
        if (old_split != new_split) {
            plan |= RECONFIGURE_PLAN_HEX_TYPES;
            if (current->geometry.center_membrane ||
                next->geometry.center_membrane) {
                plan |= RECONFIGURE_PLAN_SYNCHRONIZE |
                        RECONFIGURE_PLAN_PARTICLE_EVENTS;
            }
        }
    }

    /* Synchronizing lazy constant-field/damped trajectories updates velocities;
     * keep the legacy kinetic totals coherent without charging scalar-only
     * response edits for an unnecessary O(N) sweep. */
    if ((plan & RECONFIGURE_PLAN_SYNCHRONIZE) != 0 &&
        (current->dynamics.flight == SIM_FLIGHT_CONSTANT_FIELD ||
         current->dynamics.flight == SIM_FLIGHT_LINEAR_DAMPING)) {
        plan |= RECONFIGURE_PLAN_PHYSICAL_TOTALS;
    }

    /* These owners allocate files, baselines, or an intruder geometry.  Until
     * they have dedicated incremental bridges, explicitly route their edits
     * through the established transactional restart. */
    if (memcmp(&current->observables, &next->observables,
               sizeof(current->observables)) != 0 ||
        memcmp(&current->output, &next->output,
               sizeof(current->output)) != 0 ||
        memcmp(&current->intruder, &next->intruder,
               sizeof(current->intruder)) != 0) {
        plan |= RECONFIGURE_PLAN_FULL_RESTART;
    }

    /* Rendering is owned by the frontend and intentionally contributes no
     * core work.  Detect an unknown future field in a currently handled group
     * by comparing its representation when none of its known fields changed. */
    if (memcmp(&current->geometry, &next->geometry,
               sizeof(current->geometry)) != 0 &&
        current->geometry.boundary == next->geometry.boundary &&
        current->geometry.center_membrane == next->geometry.center_membrane &&
        !geometry_response_changed(&current->geometry, &next->geometry)) {
        plan |= RECONFIGURE_PLAN_FULL_RESTART;
    }
    if (memcmp(&current->dynamics, &next->dynamics,
               sizeof(current->dynamics)) != 0 &&
        !dynamics_known_fields_changed(&current->dynamics, &next->dynamics)) {
        plan |= RECONFIGURE_PLAN_FULL_RESTART;
    }
    if (memcmp(&current->driving, &next->driving,
               sizeof(current->driving)) != 0 &&
        current->driving.thermostat == next->driving.thermostat &&
        !changed_double(current->driving.thermostat_interval,
                        next->driving.thermostat_interval) &&
        !driving_response_changed(&current->driving, &next->driving)) {
        plan |= RECONFIGURE_PLAN_FULL_RESTART;
    }

    return plan;
}

static double canonical_coordinate(double coordinate, double length)
{
    if (!isfinite(coordinate) || !isfinite(length) || length <= 0.0)
        return NAN;
    coordinate = fmod(coordinate, length);
    if (coordinate < 0.0) coordinate += length;
    /* Keep the canonical interval half open even if rounding returns L. */
    return coordinate < length ? coordinate : 0.0;
}

static particlePrediction particle_at_current_time(const particle *value)
{
    return freeFly == &freeFlyGrow
        ? predictFreeFlyGrow(value, t)
        : predictFreeFlyNormal(value, t);
}

static bool preflight_boundary_change(const SimulationConfig *current,
                                      const SimulationConfig *next,
                                      char *message,
                                      size_t message_capacity)
{
    const bool outer_changed = current->geometry.boundary !=
        next->geometry.boundary;
    const bool membrane_added = !current->geometry.center_membrane &&
        next->geometry.center_membrane;
    const bool old_spherical = current->geometry.boundary ==
        SIM_BOUNDARY_CIRCULAR_SPHERICAL;
    const bool spherical = next->geometry.boundary ==
        SIM_BOUNDARY_CIRCULAR_SPHERICAL;
    bool outer_constraint_added = false;
    if (outer_changed) {
        if (spherical) {
            outer_constraint_added = !old_spherical;
        } else if (!old_spherical) {
            outer_constraint_added =
                (config_has_x_walls(next) &&
                 !config_has_x_walls(current)) ||
                (config_has_y_walls(next) &&
                 !config_has_y_walls(current))
#if THREE_D
                || (config_has_z_walls(next) &&
                    !config_has_z_walls(current))
#endif
                ;
        }
    }
    /* Pure removal cannot make the already valid state intersect a newly
     * introduced surface, so avoid a redundant full-particle prediction pass. */
    if (!outer_constraint_added && !membrane_added) {
        return true;
    }

    const double length_scale = fmax(1.0, fmax(Lx, Ly));
    const double tolerance = 128.0*DBL_EPSILON*length_scale;
    if (spherical) {
        if (fabs(Lx - Ly) > tolerance
#if THREE_D
            || fabs(Lx - Lz) > tolerance
#endif
        ) {
            restart_error(message, message_capacity,
                          "An in-place update cannot preserve the current "
                          "non-cubic box (%.9g x %.9g"
#if THREE_D
                          " x %.9g"
#endif
                          ") for spherical confinement; a clean rebuild from "
                          "the requested configuration is required.", Lx, Ly
#if THREE_D
                          , Lz
#endif
            );
            return false;
        }
    }

    for (int index = 0; index < N; ++index) {
        const particle *value = particles + index;
        const particlePrediction prediction =
            particle_at_current_time(value);
        const double x = canonical_coordinate(prediction.x, Lx);
        const double y = canonical_coordinate(prediction.y, Ly);
#if THREE_D
        const double z = canonical_coordinate(prediction.z, Lz);
#endif
        const double radius = prediction.rad;
        if (!isfinite(x) || !isfinite(y) || !isfinite(radius) ||
            radius < 0.0
#if THREE_D
            || !isfinite(z)
#endif
        ) {
            restart_error(message, message_capacity,
                          "An in-place update cannot preserve particle %d's "
                          "non-finite current state for the requested boundary; "
                          "a clean rebuild from the requested configuration is "
                          "required.",
                          index);
            return false;
        }

        bool contained = true;
        if (outer_constraint_added && spherical) {
            const double accessible = 0.5*Lx - radius;
            const double dx = x - 0.5*Lx;
            const double dy = y - 0.5*Ly;
            double distance_squared = dx*dx + dy*dy;
#if THREE_D
            const double dz = z - 0.5*Lz;
            distance_squared += dz*dz;
#endif
            contained = accessible >= -tolerance &&
                distance_squared <= accessible*accessible +
                    tolerance*length_scale;
        } else if (outer_constraint_added &&
                   next->geometry.boundary != SIM_BOUNDARY_PERIODIC) {
            if (config_has_x_walls(next))
                contained = contained && x + tolerance >= radius &&
                    x <= Lx - radius + tolerance;
            if (config_has_y_walls(next))
                contained = contained && y + tolerance >= radius &&
                    y <= Ly - radius + tolerance;
#if THREE_D
            if (config_has_z_walls(next))
                contained = contained && z + tolerance >= radius &&
                    z <= Lz - radius + tolerance;
#endif
        }
        if (membrane_added && value->type == 1)
            contained = contained &&
                fabs(x - 0.5*Lx) + tolerance >= radius;
        if (!contained) {
            restart_error(message, message_capacity,
                          "An in-place update cannot preserve particle %d "
                          "because it intersects the requested boundary at "
                          "t=%.9g; a clean rebuild from the requested "
                          "configuration is required.",
                          index, t);
            return false;
        }
    }
    return true;
}

static int charge_category(int charge)
{
    if (charge > 0) return 2;
    if (charge == 0) return 1;
    return 0;
}

/* Meet the requested charge counts while changing the fewest existing
 * particle labels.  No random-number stream is consumed by a GUI edit. */
static void rebalance_particle_charges(const SimPotentialConfig *potential)
{
    int desired[3];
    int counts[3] = {0, 0, 0};
    desired[2] = (int)floor(potential->positive_charge_fraction*N + 0.5);
    desired[1] = (int)floor(potential->neutral_charge_fraction*N + 0.5);
    if (desired[2] < 0) desired[2] = 0;
    if (desired[2] > N) desired[2] = N;
    if (desired[1] < 0) desired[1] = 0;
    if (desired[1] > N - desired[2]) desired[1] = N - desired[2];
    desired[0] = N - desired[1] - desired[2];

    for (int index = 0; index < N; ++index)
        counts[charge_category(particles[index].charge)]++;
    for (int index = 0; index < N; ++index) {
        const int category = charge_category(particles[index].charge);
        if (counts[category] <= desired[category]) continue;
        int replacement = -1;
        for (int candidate = 0; candidate < 3; ++candidate) {
            if (counts[candidate] < desired[candidate]) {
                replacement = candidate;
                break;
            }
        }
        if (replacement < 0) break;
        particles[index].charge = replacement - 1;
        counts[category]--;
        counts[replacement]++;
    }
}

static void apply_hex_potential_types(const SimulationConfig *next)
{
    const bool split = next->potential.enabled &&
        next->potential.species_mode == SIM_POTENTIAL_SAME_SPECIES;
    for (int index = 0; index < N; ++index)
        particles[index].type = split && index >= N/2 ? 0 : 1;
    refresh_population_counts();
}

#endif

bool simulation_restart_core(const SimulationConfig *config,
                             char *error,
                             size_t error_capacity)
{
    if (error != NULL && error_capacity > 0) {
        error[0] = '\0';
    }
    if (config == NULL) {
        restart_error(error, error_capacity,
                      "Cannot restart from a null configuration");
        return false;
    }
    if (!restart_config_is_valid(config, error, error_capacity)) {
        return false;
    }
    if (!preflight_initial_state(config, error, error_capacity)) {
        return false;
    }

#if G != 1
    restart_error(error, error_capacity,
                  "Simulation restart is only available in a graphical build");
    return false;
#else
    free_running_simulation();
    simulationClearFault();
    simulation_config_apply_globals(config);
    configureTransverseVelocityCorrelationSampling(&dtimeThermo, &logSpacing);
    if (config->initialization.mode == SIM_INIT_DUMP_FILE) {
        const size_t path_length = strlen(config->initialization.load_path);
        memcpy(initialFilePath, config->initialization.load_path,
               path_length + 1);
    }
    reset_runtime_accumulators();
    resetRuntimeShearPhase();

    init_genrand((unsigned long)(unsigned int)
                 config->initialization.random_seed);
    srand((unsigned int)config->initialization.random_seed);

    /* argc == 0 is the documented graphical-restart path: globals have
     * already been applied, and constantInit only reconstructs derived box and
     * cell-list constants. */
    constantInit(0, NULL);
    particlesInit();
    refresh_population_counts();
    adaptCellGridToParticles();
    cellListInit();
    if (addWell && U != 0.0) {
        pairsInit();
    }
    physicalQ();
#if (STRESS || HEAT)
    initializeTransportAccumulators();
#endif
    format();

    /* Output failures are non-fatal to the freshly constructed simulation.
     * outputLifecycleOpen disables only a failed stream and leaves a warning
     * in error; the GUI captures the resulting effective configuration. */
    (void)outputLifecycleOpen(error, error_capacity);
    eventListInit();

    return true;
#endif
}

SimulationReconfigureAction simulation_reconfigure_core(
    const SimulationConfig *current,
    const SimulationConfig *next,
    char *message,
    size_t message_capacity)
{
    if (message != NULL && message_capacity > 0) message[0] = '\0';
    if (current == NULL || next == NULL) {
        restart_error(message, message_capacity,
                      "Cannot update from a null configuration");
        return SIM_RECONFIGURE_FAILED;
    }
    if (!restart_config_is_valid(next, message, message_capacity))
        return SIM_RECONFIGURE_FAILED;

#if G != 1
    restart_error(message, message_capacity,
                  "In-place reconfiguration is only available in a graphical build");
    return SIM_RECONFIGURE_FAILED;
#else
    const ReconfigurePlan plan = classify_reconfiguration(current, next);
    if (plan & RECONFIGURE_PLAN_FULL_RESTART) {
        restart_error(message, message_capacity,
                      "This edit changes the initial state, box construction, "
                      "output lifecycle, observable allocation, or intruder "
                      "geometry and requires a transactional restart.");
        return SIM_RECONFIGURE_RESTART_REQUIRED;
    }
    if (plan == 0)
        return SIM_RECONFIGURE_UNCHANGED;

    const bool needs_particles = (plan &
        (RECONFIGURE_PLAN_SYNCHRONIZE |
         RECONFIGURE_PLAN_PARTICLE_EVENTS |
         RECONFIGURE_PLAN_SPATIAL_INDEX |
         RECONFIGURE_PLAN_POTENTIAL_MEMBERSHIPS |
         RECONFIGURE_PLAN_CHARGES |
         RECONFIGURE_PLAN_HEX_TYPES)) != 0;
    if (needs_particles && (particles == NULL || N <= 0)) {
        restart_error(message, message_capacity,
                      "An in-place update cannot use the uninitialized running "
                      "particle state; a clean rebuild from the requested "
                      "configuration is required.");
        return SIM_RECONFIGURE_RESTART_REQUIRED;
    }
    if ((plan & (RECONFIGURE_PLAN_PARTICLE_EVENTS |
                 RECONFIGURE_PLAN_THERMOSTAT_EVENT)) != 0 &&
        eventList == NULL) {
        restart_error(message, message_capacity,
                      "An in-place update cannot use the uninitialized running "
                      "event calendar; a clean rebuild from the requested "
                      "configuration is required.");
        return SIM_RECONFIGURE_RESTART_REQUIRED;
    }
    if ((plan & (RECONFIGURE_PLAN_SPATIAL_INDEX |
                 RECONFIGURE_PLAN_POTENTIAL_MEMBERSHIPS)) != 0 &&
        cellList == NULL) {
        restart_error(message, message_capacity,
                      "An in-place update cannot use the uninitialized running "
                      "spatial index; a clean rebuild from the requested "
                      "configuration is required.");
        return SIM_RECONFIGURE_RESTART_REQUIRED;
    }
    if ((plan & (RECONFIGURE_PLAN_SYNCHRONIZE |
                 RECONFIGURE_PLAN_PARTICLE_EVENTS |
                 RECONFIGURE_PLAN_SPATIAL_INDEX |
                 RECONFIGURE_PLAN_POTENTIAL_MEMBERSHIPS)) != 0 &&
        simulationHasFault()) {
        restart_error(message, message_capacity,
                      "An in-place structural update cannot use the fault-locked "
                      "event queue; a clean rebuild from the requested "
                      "configuration is required.");
        return SIM_RECONFIGURE_RESTART_REQUIRED;
    }
    if (!preflight_boundary_change(current, next,
                                   message, message_capacity)) {
        /* The requested configuration is valid; only the state of the old
         * trajectory prevents an in-place update.  Let the host take the
         * same preflighted restart path used for box, composition, and
         * initializer edits instead of presenting this as a rejected option. */
        return SIM_RECONFIGURE_RESTART_REQUIRED;
    }

    const bool trajectory_update =
        (plan & RECONFIGURE_PLAN_SHEAR_REBASE) != 0 ||
        trajectory_changed(&current->dynamics, &next->dynamics);
#if (STRESS || HEAT)
    if (trajectory_update)
        beginTransportVelocityBatch();
#else
    (void)trajectory_update;
#endif
    if (plan & RECONFIGURE_PLAN_SYNCHRONIZE)
        synchronizeParticlesAtCurrentTime();
    if (plan & RECONFIGURE_PLAN_SHEAR_REBASE) {
        rebaseRuntimeShear(
            current->dynamics.flight == SIM_FLIGHT_LEES_EDWARDS_SHEAR,
            current->dynamics.shear_strain_rate,
            next->dynamics.flight == SIM_FLIGHT_LEES_EDWARDS_SHEAR,
            next->dynamics.shear_strain_rate);
    }

    simulation_config_apply_globals(next);

#if (STRESS || HEAT)
    if (trajectory_update)
        rebaseTransportAccumulatorsAfterModelChange();
#endif
    if (plan & RECONFIGURE_PLAN_CHARGES)
        rebalance_particle_charges(&next->potential);
    if (plan & RECONFIGURE_PLAN_HEX_TYPES)
        apply_hex_potential_types(next);
    if (plan & RECONFIGURE_PLAN_SPATIAL_INDEX)
        rebuildSpatialIndexAtCurrentTime();
    else if (plan & RECONFIGURE_PLAN_POTENTIAL_MEMBERSHIPS)
        rebuildPotentialMembershipsAtCurrentTime();
    if (plan & RECONFIGURE_PLAN_PARTICLE_EVENTS)
        rescheduleParticleEventsAtCurrentTime();
    if (plan & RECONFIGURE_PLAN_THERMOSTAT_EVENT)
        rescheduleThermostatEventAtCurrentTime();
    if (plan & RECONFIGURE_PLAN_PHYSICAL_TOTALS)
        physicalQ();

    if (plan & RECONFIGURE_PLAN_SPATIAL_INDEX) {
        restart_error(message, message_capacity,
                      "Applied in place: particle state/time were preserved; "
                      "the cell grid, well memberships, and particle events "
                      "were refreshed.");
        return SIM_RECONFIGURE_SPATIAL_REBUILT;
    }
    if (plan & (RECONFIGURE_PLAN_PARTICLE_EVENTS |
                RECONFIGURE_PLAN_THERMOSTAT_EVENT)) {
        if (plan & RECONFIGURE_PLAN_POTENTIAL_MEMBERSHIPS) {
            restart_error(message, message_capacity,
                          "Applied in place: particle state/time and the cell "
                          "grid were preserved; only potential memberships "
                          "and affected event predictions were refreshed.");
        } else {
            restart_error(message, message_capacity,
                          "Applied in place: particle state/time were preserved "
                          "and only affected event predictions were refreshed.");
        }
        return SIM_RECONFIGURE_EVENTS_REBUILT;
    }
    restart_error(message, message_capacity,
                  "Applied immediately without rebuilding particles, cells, "
                  "or the event calendar.");
    return SIM_RECONFIGURE_GLOBALS_ONLY;
#endif
}
