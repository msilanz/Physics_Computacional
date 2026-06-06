// ============================================================================
// Von Kármán Vortex Street – Turbulent Inflow (Synthetic Eddy Method)
// 16:9 domain, 720p, 60 fps video, Multi‑Reynolds runner
// ============================================================================

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// ----------------------------------------------------------------------------
// Global parameters
// ----------------------------------------------------------------------------
double Reynolds;               // set by command line in slave processes
const double D = 0.125;        // cylinder diameter
const double U0 = 1.0;          // mean inflow velocity

// Turbulence parameters (physical values, dimensioned)
const double I = 0.1;           // turbulence intensity (u_rms / U0)
const double L_T = 0.45;        // turbulent length scale (meters) → 10% of channel height (4.5 m)
#define N_EDDIES 100            // number of synthetic eddies (compile-time constant)

int maxlevel = 9;
face vector muv[];

scalar f[];
scalar * tracers = {f};

// Video parameters
#define VIDEO_DURATION 30.0
#define VIDEO_FPS 120
#define TOTAL_FRAMES (VIDEO_DURATION * VIDEO_FPS)
#define SIMULATION_END 30.0

// ----------------------------------------------------------------------------
// Synthetic Eddy Method (SEM) data and functions
// ----------------------------------------------------------------------------
typedef struct {
    double x, y;       // physical position (m)
    double sigma;      // physical size (m)
    double epsilon;    // sign (+1 or -1)
} Eddy;

Eddy eddies[N_EDDIES];   // now size is compile-time constant

// Initialize eddies with physical positions and sizes
void init_eddies(void) {
    // Domain dimensions
    double xmin = -0.5, xmax = 7.5;
    double ymin = -2.25, ymax = 2.25;
    double y_range = ymax - ymin;   // = 4.5
    double x_range = xmax - xmin;   // = 8.0

    for (int i = 0; i < N_EDDIES; i++) {
        // Random physical position inside the domain
        eddies[i].x = xmin + ((double)rand() / RAND_MAX) * x_range;
        eddies[i].y = ymin + ((double)rand() / RAND_MAX) * y_range;
        // Random size around L_T (0.45 m) with ±50% variation
        eddies[i].sigma = L_T * (0.5 + (double)rand() / RAND_MAX);
        // Random sign (+1 or -1)
        eddies[i].epsilon = (rand() % 2) * 2 - 1;
    }
}

// Compute fluctuating velocity component at a given point (x,y) and time t
double compute_fluctuation(double x, double y, double t) {
    double u_fluct = 0.0;
    double norm = 0.0;
    double u_rms = I * U0;

    for (int i = 0; i < N_EDDIES; i++) {
        // Advect eddy with mean flow (U0)
        double x_adv = eddies[i].x - U0 * t;
        double xi = (x - x_adv) / eddies[i].sigma;
        double yi = (y - eddies[i].y) / eddies[i].sigma;
        double r2 = xi*xi + yi*yi;
        if (r2 < 1.0) {
            double f = 1.0 - r2;   // shape function
            u_fluct += eddies[i].epsilon * f;
            norm += f;
        }
    }
    if (norm > 0.0) {
        u_fluct = u_rms * u_fluct / sqrt(norm);
    }
    // Clip to ±30% of U0 to avoid unphysical bursts
    if (u_fluct > 0.3 * U0) u_fluct = 0.3 * U0;
    if (u_fluct < -0.3 * U0) u_fluct = -0.3 * U0;
    return u_fluct;
}
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Function that runs a single simulation with current global Reynolds
// ----------------------------------------------------------------------------
void run_one_simulation(void) {
    // Domain: width 8.0, height 4.5 → 16:9 aspect ratio
    L0 = 8.0;
    origin(-0.5, -L0 * (9.0 / 16.0) / 2.0);
    N = 512;
    mu = muv;

    srand(42);                // reproducible turbulence for each Reynolds
    init_eddies();            // set up eddies with physical coordinates

    display_control(Reynolds, 10, 1000);
    display_control(maxlevel, 6, 12);
    run();
}

// ----------------------------------------------------------------------------
// Basilisk events (physics and outputs)
// ----------------------------------------------------------------------------

// Viscosity (depends on Reynolds)
event properties(i++) {
    foreach_face()
        muv.x[] = fm.x[] * D * U0 / Reynolds;
}

// Inlet boundary condition with turbulent fluctuations
// The dirichlet condition is evaluated at every boundary face at each timestep.
u.n[left] = dirichlet( U0 + compute_fluctuation(x, y, t) );

// Other boundary conditions (unchanged)
u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);
u.n[top]   = neumann(0.);
u.t[top]   = neumann(0.);
u.n[bottom]= neumann(0.);
u.t[bottom]= neumann(0.);

// Inlet pressure condition
p[left]    = neumann(0.);
pf[left]   = neumann(0.);

// Tracer injection at bottom half of inlet
f[left]    = dirichlet(y < 0);

// Initial condition
event init(t = 0) {
    solid(cs, fs, sqrt(sq(x) + sq(y)) - D/2.);
    foreach()
        u.x[] = cs[] ? U0 : 0.;
}

// Logfile
event logfile(i++) {
    fprintf(stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}

// Movie output (720p, 60 fps, 30 seconds)
event movies(t += VIDEO_DURATION / TOTAL_FRAMES; t <= VIDEO_DURATION) {
    scalar omega[], mask_field[];
    vorticity(u, omega);
    foreach()
        mask_field[] = cs[] - 0.5;

    char *vort_name, *f_name;
    asprintf(&vort_name, "vort_%dfps_Re=%g_turb.mp4", VIDEO_FPS, Reynolds);
    asprintf(&f_name,   "f_%dfps_Re=%g_turb.mp4",     VIDEO_FPS, Reynolds);

    output_ppm(omega, file = vort_name,
        box = {{-0.5, -2.25}, {7.5, 2.25}},
        n = 1280,
        min = -10, max = 10,
        linear = true,
        mask = mask_field,
        opt = "crf=18 -r 60"
    );

    output_ppm(f, file = f_name,
        box = {{-0.5, -2.25}, {7.5, 2.25}},
        n = 1280,
        linear = true,
        min = 0, max = 1,
        mask = mask_field,
        opt = "crf=18 -r 60"
    );

    static int frame = 0;
    printf("Frame %d / %d\n", ++frame, (int)TOTAL_FRAMES);

    free(vort_name);
    free(f_name);
}

// Adaptive mesh refinement
event adapt(i++) {
    adapt_wavelet({cs, u, f}, {1e-2, 3e-2, 3e-2, 3e-2}, maxlevel, 4);
}

// End of simulation
event end(t = SIMULATION_END) {
    printf("Simulation finished for Re = %g. Video contains %.1f seconds.\n",
           Reynolds, VIDEO_DURATION);
    return 1;
}

// ----------------------------------------------------------------------------
// Main: Master / Slave process for multiple Reynolds numbers
// ----------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    if (argc == 2) {
        // Slave process: run one simulation with given Reynolds
        Reynolds = atof(argv[1]);
        printf("\n=== Slave: simulating Re = %g ===\n", Reynolds);
        run_one_simulation();
        return 0;
    }

    // Master process: define the Reynolds numbers to simulate
    double rey_num[] = {0.9, 20.0, 40.0, 100.0, 200.0};
    int length = sizeof(rey_num) / sizeof(rey_num[0]);

    // Build the command to call ourselves
    char *self_path;
    if (asprintf(&self_path, "./%s", argv[0]) == -1) {
        fprintf(stderr, "Memory allocation failed\n");
        return 1;
    }

    for (int i = 0; i < length; i++) {
        double Re = rey_num[i];
        printf("\n========================================\n");
        printf("Master: launching slave for Re = %g\n", Re);
        printf("========================================\n");

        char *command;
        if (asprintf(&command, "%s %g", self_path, Re) == -1) {
            fprintf(stderr, "Memory allocation failed for command\n");
            free(self_path);
            return 1;
        }

        int ret = system(command);
        if (ret != 0)
            fprintf(stderr, "Warning: slave for Re=%g failed (code %d)\n", Re, ret);

        free(command);
    }

    free(self_path);
    printf("\nAll %d simulations completed.\n", length);
    return 0;
}
