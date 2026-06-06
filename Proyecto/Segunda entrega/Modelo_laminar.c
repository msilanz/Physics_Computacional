// ============================================================================
// Von Kármán Vortex Street – 16:9 domain, 720p, 60 fps, 30 sec video
// Multi‑Reynolds runner: if called with an argument, run one simulation.
// If called without arguments, loop over the rey_num[] array and spawn itself.
// ============================================================================

#include "embed.h"
#include "navier-stokes/centered.h"
#include "tracer.h"
#include <stdio.h>
#include <stdlib.h>   // for atof, system, snprintf

// ----------------------------------------------------------------------------
// Global parameters (these will be set before each simulation)
// ----------------------------------------------------------------------------
double Reynolds;             // default; will be overridden by command‑line arg
const double D = 0.125;      // cylinder diameter
const double U0 = 1.0;       // inflow velocity

int maxlevel = 9;
face vector muv[];

scalar f[];
scalar * tracers = {f};

// Video parameters – do not change
#define VIDEO_DURATION 30.0
#define VIDEO_FPS 120
#define TOTAL_FRAMES (VIDEO_DURATION * VIDEO_FPS)
#define SIMULATION_END 30.0

// ----------------------------------------------------------------------------
// Function that runs a SINGLE simulation with the current global 'Reynolds'
// ----------------------------------------------------------------------------
void run_one_simulation() {
    // Domain size: width = 8.0, height = 4.5  → aspect ratio 8:4.5 = 16:9
    L0 = 8.0;
    // Origin: x from -0.5 to 7.5, y from -2.25 to 2.25
    origin(-0.5, -L0 * (9.0 / 16.0) / 2.0); // -L0*(9/16)/2 = -2.25
    N = 512;
    mu = muv;

    // Optional interactive controls
    display_control(Reynolds, 10, 1000);
    display_control(maxlevel, 6, 12);

    run();   // starts the event loop
}

// ----------------------------------------------------------------------------
// Event definitions (these are global – they will be used by run())
// ----------------------------------------------------------------------------

// Viscosity (constant, based on Reynolds)
event properties(i++) {
    foreach_face()
        muv.x[] = fm.x[] * D * U0 / Reynolds;
}

// Boundary conditions

// Inlet (left)
u.n[left]  = dirichlet(U0);
p[left]    = neumann(0.);
pf[left]   = neumann(0.);
f[left]    = dirichlet(y < 0); // tracer injected in bottom half

// Outlet (right)
u.n[right] = neumann(0.);
p[right]   = dirichlet(0.);
pf[right]  = dirichlet(0.);

// Cylinder surface (no-slip)
u.n[embed] = dirichlet(0.);
u.t[embed] = dirichlet(0.);

// Top and bottom walls (free-slip – adjust as needed)
u.n[top]    = neumann(0.);
u.t[top]    = neumann(0.);
u.n[bottom] = neumann(0.);
u.t[bottom] = neumann(0.);

event init(t = 0) {
    solid(cs, fs, sqrt(sq(x) + sq(y)) - D/2.);
    foreach()
        u.x[] = cs[] ? U0 : 0.;
}

event logfile(i++) {
    fprintf(stderr, "%d %g %d %d\n", i, t, mgp.i, mgu.i);
}

event movies(t += VIDEO_DURATION / TOTAL_FRAMES; t <= VIDEO_DURATION) {
    scalar omega[], mask_field[];
    vorticity(u, omega);
    foreach()
        mask_field[] = cs[] - 0.5;

    // File names depend on Reynolds and FPS – constructed with snprintf
    char vort_name[2048], f_name[2048];
    snprintf(vort_name, sizeof(vort_name), "vort_%dfps_Re=%g.mp4", VIDEO_FPS, Reynolds);
    snprintf(f_name,   sizeof(f_name),   "f_%dfps_Re=%g.mp4",      VIDEO_FPS, Reynolds);

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
}

// Adaptive mesh refinement
event adapt(i++) {
    adapt_wavelet({cs, u, f}, {1e-2, 3e-2, 3e-2, 3e-2}, maxlevel, 4);
}

event end(t = SIMULATION_END) {
    printf("Simulation finished for Re = %g. Video contains %.1f seconds.\n",
           Reynolds, VIDEO_DURATION);
    return 1;
}

// Main: master process loops over the array, launches slave processes
// ----------------------------------------------------------------------------
int main(int argc, char *argv[]) {
    if (argc == 2) {
        // slave process: set Reynolds from argument and run simulation
        Reynolds = atof(argv[1]);
        printf("\n=== slave: simulating Re = %g ===\n", Reynolds);
        run_one_simulation();
        return 0;
    }

    // Master process: define the array of Reynolds numbers (only here!)
    double rey_num[] = {0.9, 20.0, 40.0, 100.0, 200.0};
    int length = sizeof(rey_num) / sizeof(rey_num[0]);

    // Get current executable path
    char self_path[2048];
    snprintf(self_path, sizeof(self_path), "./%s", argv[0]);

    for (int i = 0; i < length; i++) {
        double Re = rey_num[i];
        printf("\n========================================\n");
        printf("Master: launching slave for Re = %g\n", Re);
        printf("========================================\n");

        char command[2048];
        snprintf(command, sizeof(command), "%s %g", self_path, Re);
        int ret = system(command);
        if (ret != 0)
            fprintf(stderr, "Warning: slave for Re=%g failed (code %d)\n", Re, ret);
    }

    printf("\nAll %d simulations completed.\n", length);
    return 0;
}
