/*
 * blackhole.c - Schwarzschild black hole ray tracer
 *
 * Renders a non-rotating black hole with a thin accretion disk by tracing
 * light rays backwards from the camera along curved (null geodesic) paths.
 *
 * Physics (units: Schwarzschild radius r_s = 1):
 *   A photon's path in Schwarzschild spacetime obeys, in Cartesian form,
 *       d^2x/dl^2 = -1.5 * h^2 * x / r^5
 *   where h = |x cross v| is the (conserved) specific angular momentum.
 *   This reproduces the exact orbit equation u'' + u = 1.5 u^2 (u = 1/r).
 *
 *   Accretion disk: Shakura-Sunyaev temperature profile from the ISCO
 *   (r = 3) outward, Keplerian orbital velocity, relativistic Doppler
 *   shift + beaming and gravitational redshift (g^4 intensity scaling).
 *
 * Output: binary PPM (P6) image, viewable with most image tools.
 *
 * Build:  gcc -O3 -fopenmp -o blackhole blackhole.c -lm
 *         (omit -fopenmp for single-threaded)
 * Usage:  ./blackhole [-w width] [-h height] [-i incl_deg] [-d distance]
 *                     [-f fov_deg] [-n frames] [-o output.ppm]
 *   With -n N > 1 the camera orbits the hole, writing frame_000.ppm ...
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ---------- small vector library ---------- */
typedef struct { double x, y, z; } vec3;

static inline vec3 v3(double x, double y, double z) { vec3 r = {x, y, z}; return r; }
static inline vec3 add(vec3 a, vec3 b) { return v3(a.x + b.x, a.y + b.y, a.z + b.z); }
static inline vec3 sub(vec3 a, vec3 b) { return v3(a.x - b.x, a.y - b.y, a.z - b.z); }
static inline vec3 mul(vec3 a, double s) { return v3(a.x * s, a.y * s, a.z * s); }
static inline double dot(vec3 a, vec3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static inline vec3 cross(vec3 a, vec3 b) {
    return v3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}
static inline double len(vec3 a) { return sqrt(dot(a, a)); }
static inline vec3 norm(vec3 a) { return mul(a, 1.0 / len(a)); }

/* ---------- scene parameters ---------- */
#define DISK_IN   3.0     /* ISCO = 3 r_s */
#define DISK_OUT  14.0
#define ESCAPE_R  60.0
#define MAX_STEPS 4000
#define DISK_T0   4500.0  /* peak-ish disk temperature in Kelvin (artistic) */

/* ---------- hashing / noise ---------- */
static inline uint32_t hash3(int x, int y, int z) {
    uint32_t h = (uint32_t)x * 374761393u + (uint32_t)y * 668265263u + (uint32_t)z * 2147483647u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return h ^ (h >> 16);
}
static inline double hashf(int x, int y, int z) { return (hash3(x, y, z) & 0xffffff) / 16777216.0; }

static double vnoise2(double x, double y, int seed) {
    int xi = (int)floor(x), yi = (int)floor(y);
    double fx = x - xi, fy = y - yi;
    double sx = fx * fx * (3 - 2 * fx), sy = fy * fy * (3 - 2 * fy);
    double a = hashf(xi, yi, seed), b = hashf(xi + 1, yi, seed);
    double c = hashf(xi, yi + 1, seed), d = hashf(xi + 1, yi + 1, seed);
    return (a + (b - a) * sx) + ((c + (d - c) * sx) - (a + (b - a) * sx)) * sy;
}

/* ---------- blackbody temperature -> linear RGB (normalised) ---------- */
static vec3 blackbody(double T) {
    /* Tanner Helland's fit, returns sRGB-ish in [0,1]; we treat as linear */
    double t = T / 100.0, r, g, b;
    if (t <= 66) {
        r = 255;
        g = 99.4708025861 * log(t) - 161.1195681661;
        b = t <= 19 ? 0 : 138.5177312231 * log(t - 10) - 305.0447927307;
    } else {
        r = 329.698727446 * pow(t - 60, -0.1332047592);
        g = 288.1221695283 * pow(t - 60, -0.0755148492);
        b = 255;
    }
#define CL(v) ((v) < 0 ? 0 : (v) > 255 ? 255 : (v))
    return v3(CL(r) / 255.0, CL(g) / 255.0, CL(b) / 255.0);
#undef CL
}

/* ---------- background sky: procedural stars + faint galaxy band ---------- */
static vec3 sky(vec3 d) {
    vec3 col = v3(0, 0, 0);
    /* stars: hash a direction grid at a few scales */
    for (int layer = 0; layer < 3; layer++) {
        double scale = 120.0 + layer * 110.0;
        vec3 p = mul(d, scale);
        int ix = (int)floor(p.x), iy = (int)floor(p.y), iz = (int)floor(p.z);
        double h = hashf(ix, iy, iz + layer * 1000);
        if (h > 0.993) {
            vec3 c = v3(ix + 0.5, iy + 0.5, iz + 0.5);
            double dist = len(sub(p, c));
            double bright = pow(h - 0.993, 1.0) / 0.007;
            double fall = exp(-dist * dist * 18.0);
            double temp = 3000 + 9000 * hashf(iz, ix, iy + 7);
            col = add(col, mul(blackbody(temp), bright * fall * 2.5));
        }
    }
    /* milky-way-like band along a tilted great circle */
    vec3 bandN = norm(v3(0.3, 1.0, 0.45));
    double bd = dot(d, bandN);
    double band = exp(-bd * bd * 12.0);
    double neb = vnoise2(atan2(d.z, d.x) * 6.0, d.y * 6.0, 3) * 0.6
               + vnoise2(atan2(d.z, d.x) * 14.0, d.y * 14.0, 5) * 0.4;
    col = add(col, mul(v3(0.35, 0.28, 0.45), band * neb * 0.05));
    return col;
}

/* ---------- accretion disk emission at a hit point ---------- */
/* returns emitted (observed) colour and sets *alpha (opacity) */
static vec3 disk_color(vec3 p, vec3 ray_dir, double *alpha) {
    double r = sqrt(p.x * p.x + p.z * p.z);
    /* Shakura-Sunyaev: T ~ r^-3/4 (1 - sqrt(r_in/r))^1/4 */
    double prof = pow(r / DISK_IN, -0.75) * pow(1.0 - sqrt(DISK_IN / r), 0.25);
    double T = DISK_T0 * prof / 0.488; /* normalise peak (~0.488) to T0 */

    /* Keplerian orbit (counter-clockwise viewed from +y), speed measured by
       a static observer: beta = sqrt(M/(r-2M)) = sqrt(0.5/(r-1)) */
    double beta = sqrt(0.5 / (r - 1.0));
    vec3 vdir = norm(v3(-p.z, 0, p.x));
    double gamma = 1.0 / sqrt(1.0 - beta * beta);
    /* photon travels opposite to our backward ray */
    vec3 photon = mul(ray_dir, -1.0);
    double doppler = 1.0 / (gamma * (1.0 - beta * dot(vdir, photon)));
    double grav = sqrt(1.0 - 1.0 / r);
    double g = doppler * grav;

    double Tobs = T * g;
    double intensity = pow(g, 4.0) * prof * prof;

    /* turbulent streaks: noise in (log r, angle) co-rotating space */
    double ang = atan2(p.z, p.x);
    double n = 0.55 * vnoise2(log(r) * 18.0, ang * 3.0 / M_PI * 4.0, 11)
             + 0.30 * vnoise2(log(r) * 45.0, ang * 3.0 / M_PI * 10.0, 17)
             + 0.15 * vnoise2(log(r) * 110.0, ang * 3.0 / M_PI * 22.0, 23);
    double streak = 0.35 + 1.1 * n * n;

    /* soft edges */
    double edge = fmin(1.0, (r - DISK_IN) / 0.4) * fmin(1.0, (DISK_OUT - r) / 4.0);
    if (edge < 0) edge = 0;

    *alpha = fmin(1.0, (0.55 + 0.6 * n) * edge);
    /* boost saturation a little so the disk reads as hot plasma, not grey */
    vec3 bb = blackbody(Tobs);
    bb = v3(pow(bb.x, 1.8), pow(bb.y, 1.8), pow(bb.z, 1.8));
    return mul(bb, intensity * streak * 1.1);
}

/* ---------- geodesic integration ---------- */
static inline vec3 accel(vec3 p, double h2) {
    double r2 = dot(p, p);
    double r5 = r2 * r2 * sqrt(r2);
    return mul(p, -1.5 * h2 / r5);
}

static vec3 trace(vec3 pos, vec3 dir) {
    vec3 col = v3(0, 0, 0);
    double transmit = 1.0;
    double h2 = dot(cross(pos, dir), cross(pos, dir));
    vec3 p = pos, v = dir;

    for (int i = 0; i < MAX_STEPS; i++) {
        double r = len(p);
        if (r < 1.0) return col;                 /* fell past the horizon */
        if (r > ESCAPE_R && dot(p, v) > 0)       /* escaped to infinity */
            return add(col, mul(sky(norm(v)), transmit));

        /* adaptive step: small near the hole, larger far away */
        double dt = fmin(0.08 * r, 0.02 + 0.05 * r * r);
        if (dt > 2.0) dt = 2.0;

        /* RK4 on (p, v) */
        vec3 k1v = accel(p, h2),                       k1p = v;
        vec3 k2v = accel(add(p, mul(k1p, dt / 2)), h2), k2p = add(v, mul(k1v, dt / 2));
        vec3 k3v = accel(add(p, mul(k2p, dt / 2)), h2), k3p = add(v, mul(k2v, dt / 2));
        vec3 k4v = accel(add(p, mul(k3p, dt)), h2),     k4p = add(v, mul(k3v, dt));
        vec3 np = add(p, mul(add(add(k1p, mul(k2p, 2)), add(mul(k3p, 2), k4p)), dt / 6));
        vec3 nv = add(v, mul(add(add(k1v, mul(k2v, 2)), add(mul(k3v, 2), k4v)), dt / 6));

        /* disk plane crossing (y = 0) */
        if ((p.y > 0) != (np.y > 0)) {
            double t = p.y / (p.y - np.y);
            vec3 hit = add(p, mul(sub(np, p), t));
            double rh = sqrt(hit.x * hit.x + hit.z * hit.z);
            if (rh > DISK_IN && rh < DISK_OUT) {
                double a;
                vec3 c = disk_color(hit, norm(v), &a);
                col = add(col, mul(c, transmit * a));
                transmit *= (1.0 - a);
                if (transmit < 0.01) return col;
            }
        }
        p = np;
        v = nv;
    }
    return col;
}

/* ---------- tone mapping (ACES approximation) + sRGB ---------- */
static double aces(double x) {
    const double a = 2.51, b = 0.03, c = 2.43, d = 0.59, e = 0.14;
    double y = (x * (a * x + b)) / (x * (c * x + d) + e);
    return y < 0 ? 0 : y > 1 ? 1 : y;
}
static uint8_t to_srgb8(double x) {
    x = aces(x);
    x = x <= 0.0031308 ? 12.92 * x : 1.055 * pow(x, 1.0 / 2.4) - 0.055;
    int v = (int)(x * 255.0 + 0.5);
    return (uint8_t)(v < 0 ? 0 : v > 255 ? 255 : v);
}

/* ---------- simple bloom: blur bright parts and add back ---------- */
static void bloom(vec3 *img, int W, int H) {
    vec3 *bright = malloc(sizeof(vec3) * W * H), *tmp = malloc(sizeof(vec3) * W * H);
    for (int i = 0; i < W * H; i++) {
        vec3 c = img[i];
        double l = 0.2126 * c.x + 0.7152 * c.y + 0.0722 * c.z;
        double k = l > 1.0 ? (l - 1.0) / l : 0.0;
        bright[i] = mul(c, k);
    }
    int R = W / 60 + 2;
    double sigma = R / 2.5;
    double *wts = malloc(sizeof(double) * (2 * R + 1)), wsum = 0;
    for (int k = -R; k <= R; k++) wsum += (wts[k + R] = exp(-k * k / (2 * sigma * sigma)));
    for (int k = 0; k <= 2 * R; k++) wts[k] /= wsum;
    #pragma omp parallel for
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            vec3 s = v3(0, 0, 0);
            for (int k = -R; k <= R; k++) {
                int xx = x + k; xx = xx < 0 ? 0 : xx >= W ? W - 1 : xx;
                s = add(s, mul(bright[y * W + xx], wts[k + R]));
            }
            tmp[y * W + x] = s;
        }
    #pragma omp parallel for
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            vec3 s = v3(0, 0, 0);
            for (int k = -R; k <= R; k++) {
                int yy = y + k; yy = yy < 0 ? 0 : yy >= H ? H - 1 : yy;
                s = add(s, mul(tmp[yy * W + x], wts[k + R]));
            }
            img[y * W + x] = add(img[y * W + x], mul(s, 0.6));
        }
    free(bright); free(tmp); free(wts);
}

/* ---------- render one frame ---------- */
static int render(const char *path, int W, int H, double incl_deg, double azim_deg,
                  double dist, double fov_deg, int spp) {
    double inc = incl_deg * M_PI / 180.0, az = azim_deg * M_PI / 180.0;
    /* camera above the disk plane by `incl` degrees, looking at the origin */
    vec3 cam = v3(dist * cos(inc) * sin(az), dist * sin(inc), -dist * cos(inc) * cos(az));
    vec3 fwd = norm(mul(cam, -1.0));
    vec3 right = norm(cross(v3(0, 1, 0), fwd));
    vec3 up = cross(fwd, right);
    double scale = tan(fov_deg * M_PI / 360.0);

    vec3 *img = malloc(sizeof(vec3) * W * H);
    if (!img) return -1;

    #pragma omp parallel for schedule(dynamic, 1)
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            vec3 acc = v3(0, 0, 0);
            for (int s = 0; s < spp; s++) {
                /* stratified 2x2 / NxN jitter for anti-aliasing */
                int side = (int)sqrt((double)spp);
                double jx = side > 1 ? ((s % side) + 0.5) / side : 0.5;
                double jy = side > 1 ? ((s / side) + 0.5) / side : 0.5;
                double u = (2.0 * (x + jx) / W - 1.0) * scale * W / H;
                double v = (1.0 - 2.0 * (y + jy) / H) * scale;
                vec3 dir = norm(add(fwd, add(mul(right, u), mul(up, v))));
                acc = add(acc, trace(cam, dir));
            }
            img[y * W + x] = mul(acc, 1.0 / spp);
        }
    }

    bloom(img, W, H);

    FILE *f = fopen(path, "wb");
    if (!f) { perror(path); free(img); return -1; }
    fprintf(f, "P6\n%d %d\n255\n", W, H);
    uint8_t *row = malloc(3 * W);
    for (int y = 0; y < H; y++) {
        for (int x = 0; x < W; x++) {
            vec3 c = mul(img[y * W + x], 1.2); /* exposure */
            row[3 * x] = to_srgb8(c.x);
            row[3 * x + 1] = to_srgb8(c.y);
            row[3 * x + 2] = to_srgb8(c.z);
        }
        fwrite(row, 3, W, f);
    }
    free(row);
    fclose(f);
    free(img);
    return 0;
}

static void usage(const char *prog) {
    fprintf(stderr,
        "Usage: %s [options]\n"
        "  -w N     image width        (default 1280)\n"
        "  -h N     image height       (default 720)\n"
        "  -i DEG   camera inclination above disk (default 8)\n"
        "  -a DEG   camera azimuth     (default 0)\n"
        "  -d R     camera distance in r_s (default 30)\n"
        "  -f DEG   vertical field of view (default 40)\n"
        "  -s N     samples per pixel, square number (default 4)\n"
        "  -n N     number of frames for an orbit animation (default 1)\n"
        "  -o FILE  output file (default blackhole.ppm)\n", prog);
}

int main(int argc, char **argv) {
    int W = 1280, H = 720, spp = 4, frames = 1;
    double incl = 8.0, azim = 0.0, dist = 30.0, fov = 40.0;
    const char *out = "blackhole.ppm";
    int opt;
    while ((opt = getopt(argc, argv, "w:h:i:a:d:f:s:n:o:")) != -1) {
        switch (opt) {
        case 'w': W = atoi(optarg); break;
        case 'h': H = atoi(optarg); break;
        case 'i': incl = atof(optarg); break;
        case 'a': azim = atof(optarg); break;
        case 'd': dist = atof(optarg); break;
        case 'f': fov = atof(optarg); break;
        case 's': spp = atoi(optarg); break;
        case 'n': frames = atoi(optarg); break;
        case 'o': out = optarg; break;
        default: usage(argv[0]); return 1;
        }
    }
    if (W <= 0 || H <= 0 || spp <= 0 || frames <= 0 || dist <= 1.5) { usage(argv[0]); return 1; }

    if (frames == 1) {
        fprintf(stderr, "Rendering %dx%d (%d spp) -> %s\n", W, H, spp, out);
        return render(out, W, H, incl, azim, dist, fov, spp) ? 1 : 0;
    }
    for (int i = 0; i < frames; i++) {
        char name[64];
        snprintf(name, sizeof name, "frame_%03d.ppm", i);
        double a = azim + 360.0 * i / frames;
        fprintf(stderr, "\rFrame %d/%d", i + 1, frames);
        if (render(name, W, H, incl, a, dist, fov, spp)) return 1;
    }
    fprintf(stderr, "\nDone. Make a GIF with e.g.: convert -delay 4 frame_*.ppm orbit.gif\n");
    return 0;
}
