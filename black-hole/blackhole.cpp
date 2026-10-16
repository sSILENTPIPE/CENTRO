// blackhole.cpp - a ray-traced Schwarzschild black hole with an accretion disk.
//
// Every pixel sends a photon backwards from the camera and follows it along a
// null geodesic of the Schwarzschild metric until it falls through the event
// horizon, hits the accretion disk or escapes to the sky. Nothing in the
// picture is painted in by hand; it all comes out of the geodesics:
//
//   * the black "shadow", which is 2.6 times wider than the horizon itself,
//   * the thin photon ring of light that orbited the hole before escaping,
//   * the far side of the disk, lensed up over the top and under the bottom,
//   * Doppler beaming: the side of the disk coming towards us is bluer and
//     far brighter than the side moving away,
//   * gravitational redshift, which dims and reddens the inner disk,
//   * an Einstein ring and multiple images of the background stars.
//
// Units: G = c = 1 and the mass of the hole M = 1, so lengths are measured in
// GM/c^2 and times in GM/c^3. Horizon r = 2, photon sphere r = 3, ISCO r = 6.
//
// Build:  g++ -O3 -std=c++17 -pthread blackhole.cpp -o blackhole
// Run:    ./blackhole --help

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <thread>
#include <vector>

namespace {

// ---------------------------------------------------------------------------
// Constants

constexpr double kPi = 3.14159265358979323846;
constexpr double kM = 1.0;                  // mass of the black hole
constexpr double kHorizon = 2.0 * kM;       // Schwarzschild radius r_s
constexpr double kPhotonSphere = 3.0 * kM;  // unstable circular orbits of light
constexpr double kIsco = 6.0 * kM;          // innermost stable circular orbit
const double kCriticalB = 3.0 * std::sqrt(3.0) * kM;  // impact parameter of the shadow edge

constexpr double kEscapeRadius = 1000.0;  // a photon this far out and receding has reached the sky
constexpr int kMaxSteps = 20000;

// Physical size of the unit of length and time for a hole of one solar mass.
constexpr double kMetresPerMsun = 1476.625;      // G M_sun / c^2
constexpr double kSecondsPerMsun = 4.925491e-6;  // G M_sun / c^3

// ---------------------------------------------------------------------------
// Small vector type, also used for linear RGB colours

struct Vec3 {
  double x = 0, y = 0, z = 0;
  Vec3() = default;
  constexpr Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
  Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
  Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
  Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
  Vec3 operator/(double s) const { return {x / s, y / s, z / s}; }
  Vec3& operator+=(const Vec3& o) {
    x += o.x;
    y += o.y;
    z += o.z;
    return *this;
  }
};

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double length(const Vec3& v) { return std::sqrt(dot(v, v)); }
inline Vec3 normalize(const Vec3& v) { return v / length(v); }
inline Vec3 mix(const Vec3& a, const Vec3& b, double t) { return a * (1.0 - t) + b * t; }
inline double luminance(const Vec3& c) { return 0.2126 * c.x + 0.7152 * c.y + 0.0722 * c.z; }

inline double smoothstep(double e0, double e1, double x) {
  double t = std::clamp((x - e0) / (e1 - e0), 0.0, 1.0);
  return t * t * (3.0 - 2.0 * t);
}

// ---------------------------------------------------------------------------
// Light rays
//
// A photon in Schwarzschild spacetime moves in a plane through the centre,
// and its orbit u(phi) = 1/r obeys
//
//     d^2u/dphi^2 + u = 3 M u^2.
//
// The Cartesian equation of motion
//
//     d^2x/dl^2 = -3 M |h|^2 x / |x|^5,     h = x cross dx/dl  (conserved),
//
// is a central force whose orbits satisfy exactly that equation (by Binet's
// formula), so rays can be integrated in plain 3-D Cartesian coordinates with
// none of the pole singularities of spherical ones. The parameter l is not
// the affine parameter, but only the shape of the path matters for imaging.
// Along the path |v|^2 - 2M|h|^2/r^3 stays constant and equals (|h|/b)^2,
// where b is the photon's impact parameter.

struct Photon {
  Vec3 pos, vel;
  Vec3 h;          // x cross v, normal to the plane of the orbit
  double k = 0.0;  // 3 M |h|^2: strength of the bending (0 in flat space)

  Photon(const Vec3& p, const Vec3& v, bool gravity)
      : pos(p), vel(v), h(cross(p, v)), k(gravity ? 3.0 * kM * dot(h, h) : 0.0) {}

  Vec3 accel(const Vec3& p) const {
    double r2 = dot(p, p);
    return p * (-k / (r2 * r2 * std::sqrt(r2)));
  }

  // One RK4 step of length about frac * r. Out in nearly flat space a ray
  // heading almost straight in or out (small angle psi to the radial
  // direction) barely bends, since the force is central, so there the step
  // stretches by up to 1/sin(psi); near the hole, and at every ray's closest
  // approach, it stays a fixed fraction of r.
  void step(double frac) {
    double r = length(pos), speed = length(vel);
    double sin_psi = length(h) / (r * speed);
    double stretch = std::min(1.0 / std::max(sin_psi, 0.01), 1.0 + r / 50.0);
    double ds = std::min(frac * r * stretch, 0.5 * r);
    double dt = ds / speed;
    Vec3 a1 = accel(pos);
    Vec3 p2 = pos + vel * (0.5 * dt), v2 = vel + a1 * (0.5 * dt);
    Vec3 a2 = accel(p2);
    Vec3 p3 = pos + v2 * (0.5 * dt), v3 = vel + a2 * (0.5 * dt);
    Vec3 a3 = accel(p3);
    Vec3 p4 = pos + v3 * dt, v4 = vel + a3 * dt;
    Vec3 a4 = accel(p4);
    pos += (vel + v2 * 2.0 + v3 * 2.0 + v4) * (dt / 6.0);
    vel += (a1 + a2 * 2.0 + a3 * 2.0 + a4) * (dt / 6.0);
  }
};

// A static observer at radius r sees a photon arrive at angle alpha from the
// radial direction. Its impact parameter is b = r sin(alpha) / sqrt(1 - 2M/r),
// and the same geodesic comes out of the Cartesian model by shrinking the
// radial part of the direction by sqrt(1 - 2M/r). With that, every ray leaving
// the observer has |v|^2 - 2M|h|^2/r^3 = 1 - 2M/r.
Vec3 observer_to_ray(const Vec3& pos, const Vec3& dir, bool gravity) {
  if (!gravity) return dir;
  double r = length(pos);
  Vec3 rhat = pos / r;
  return dir + rhat * (dot(dir, rhat) * (std::sqrt(1.0 - 2.0 * kM / r) - 1.0));
}

enum class Fate { Captured, Escaped, Absorbed };

// Steps a photon until it crosses the horizon, escapes past escape_radius or
// on_step(previous position, previous velocity) returns false.
template <class OnStep>
Fate propagate(Photon& ph, double step, double escape_radius, OnStep&& on_step) {
  for (int i = 0; i < kMaxSteps; ++i) {
    Vec3 p0 = ph.pos, v0 = ph.vel;
    ph.step(step);
    if (!on_step(p0, v0)) return Fate::Absorbed;
    double r2 = dot(ph.pos, ph.pos);
    if (r2 < kHorizon * kHorizon) return Fate::Captured;
    if (r2 > escape_radius * escape_radius && dot(ph.pos, ph.vel) > 0.0) return Fate::Escaped;
  }
  return Fate::Captured;  // still circling the photon sphere: it is as good as gone
}

// ---------------------------------------------------------------------------
// Colour of thermal light
//
// A blackbody seen with redshift g = nu_obs / nu_emit still looks exactly like
// a blackbody, at temperature g T (the invariance of I_nu / nu^3 makes all the
// powers of g cancel). So the observed colour and brightness of a piece of the
// disk is the Planck spectrum at g T integrated against the CIE 1931
// colour-matching functions; beaming (the g^4 in total brightness) is built in.

double cie_lobe(double nm, double mu, double s1, double s2) {
  double t = (nm - mu) / (nm < mu ? s1 : s2);
  return std::exp(-0.5 * t * t);
}

// Multi-lobe Gaussian fit to the CIE 1931 2-degree observer
// (Wyman, Sloan & Shirley 2013).
Vec3 cie_xyz(double nm) {
  return {1.056 * cie_lobe(nm, 599.8, 37.9, 31.0) + 0.362 * cie_lobe(nm, 442.0, 16.0, 26.7) -
              0.065 * cie_lobe(nm, 501.1, 20.4, 26.2),
          0.821 * cie_lobe(nm, 568.8, 46.9, 40.5) + 0.286 * cie_lobe(nm, 530.9, 16.3, 31.1),
          1.217 * cie_lobe(nm, 437.0, 11.8, 36.0) + 0.681 * cie_lobe(nm, 459.0, 26.0, 13.8)};
}

Vec3 xyz_to_linear_srgb(const Vec3& c) {
  return {3.2406 * c.x - 1.5372 * c.y - 0.4986 * c.z, -0.9689 * c.x + 1.8758 * c.y + 0.0415 * c.z,
          0.0557 * c.x - 0.2040 * c.y + 1.0570 * c.z};
}

class Blackbody {
 public:
  Blackbody() : table_(kSize) {
    for (int i = 0; i < kSize; ++i) table_[i] = planck_rgb(temperature(i));
  }

  // Linear sRGB radiance of a blackbody at temperature T, in arbitrary units.
  Vec3 radiance(double T) const {
    if (T <= kTmin) return {};
    double f = std::log(T / kTmin) / std::log(kTmax / kTmin) * (kSize - 1);
    if (f >= kSize - 1) return table_.back() * (T / kTmax);  // Rayleigh-Jeans tail
    int i = static_cast<int>(f);
    return mix(table_[i], table_[i + 1], f - i);
  }

  // Just the colour, scaled to unit luminance.
  Vec3 tint(double T) const {
    Vec3 c = radiance(T);
    return c / luminance(c);
  }

 private:
  static constexpr int kSize = 2048;
  static constexpr double kTmin = 300.0, kTmax = 1e6;

  static double temperature(int i) { return kTmin * std::pow(kTmax / kTmin, i / double(kSize - 1)); }

  static Vec3 planck_rgb(double T) {
    constexpr double c2 = 1.4387769e-2;  // h c / k_B in m K
    Vec3 xyz;
    for (double nm = 360.0; nm <= 830.0; nm += 5.0) {
      double lam = nm * 1e-9;
      xyz += cie_xyz(nm) * (1e-30 / (std::pow(lam, 5) * std::expm1(c2 / (lam * T))));
    }
    Vec3 rgb = xyz_to_linear_srgb(xyz);
    return {std::max(rgb.x, 0.0), std::max(rgb.y, 0.0), std::max(rgb.z, 0.0)};
  }

  std::vector<Vec3> table_;
};

// ---------------------------------------------------------------------------
// Hashing and value noise (disk turbulence, stars, the Milky Way)

inline uint32_t hash32(uint32_t x) {
  x ^= x >> 16;
  x *= 0x7feb352dU;
  x ^= x >> 15;
  x *= 0x846ca68bU;
  x ^= x >> 16;
  return x;
}

inline uint32_t hash3(int x, int y, int z, uint32_t seed = 0) {
  return hash32(static_cast<uint32_t>(x) * 0x8da6b343U ^
                hash32(static_cast<uint32_t>(y) * 0xd8163841U ^
                       hash32(static_cast<uint32_t>(z) * 0xcb1ab31fU ^ seed)));
}

inline double to_unit(uint32_t h) { return (h >> 8) * (1.0 / 16777216.0); }

double value_noise(const Vec3& p) {
  double fx = std::floor(p.x), fy = std::floor(p.y), fz = std::floor(p.z);
  int ix = static_cast<int>(fx), iy = static_cast<int>(fy), iz = static_cast<int>(fz);
  double tx = p.x - fx, ty = p.y - fy, tz = p.z - fz;
  tx = tx * tx * (3.0 - 2.0 * tx);
  ty = ty * ty * (3.0 - 2.0 * ty);
  tz = tz * tz * (3.0 - 2.0 * tz);
  auto v = [&](int dx, int dy, int dz) { return to_unit(hash3(ix + dx, iy + dy, iz + dz)); };
  double x00 = v(0, 0, 0) + (v(1, 0, 0) - v(0, 0, 0)) * tx;
  double x10 = v(0, 1, 0) + (v(1, 1, 0) - v(0, 1, 0)) * tx;
  double x01 = v(0, 0, 1) + (v(1, 0, 1) - v(0, 0, 1)) * tx;
  double x11 = v(0, 1, 1) + (v(1, 1, 1) - v(0, 1, 1)) * tx;
  double y0 = x00 + (x10 - x00) * ty;
  double y1 = x01 + (x11 - x01) * ty;
  return y0 + (y1 - y0) * tz;
}

// Fractal sum of octaves, normalised to [0, 1].
double fbm(Vec3 p, int octaves) {
  double sum = 0.0, amp = 0.5, norm = 0.0;
  for (int i = 0; i < octaves; ++i) {
    sum += amp * value_noise(p);
    norm += amp;
    amp *= 0.5;
    p = p * 2.03 + Vec3(17.1, 5.3, 11.7);
  }
  return sum / norm;
}

// ---------------------------------------------------------------------------
// Settings

enum class Sky { Stars, Grid, None };
enum class Shift { Full, Gravity, None };

struct Settings {
  int width = 1280, height = 720;
  double fov_deg = 48.0;        // vertical field of view
  double distance = 30.0;       // camera radius, in M
  double incl_deg = 84.0;       // angle between the line of sight and the disk axis
  double azimuth_deg = 0.0;     // camera position around the disk axis
  bool disk = true;
  double disk_in = kIsco;       // inner edge of the disk, in M
  double disk_out = 20.0;       // outer edge of the disk, in M
  double disk_temp = 4500.0;    // hottest temperature in the disk, in its own frame (K)
  double exposure = 1.0;
  int aa = 2;                   // aa x aa samples per pixel
  double step = 0.03;           // integration step as a fraction of r
  Sky sky = Sky::Stars;
  Shift shift = Shift::Full;
  bool flat = false;            // straight light rays, for comparison
  bool bloom = true;
  int frames = 1;
  double frame_dt = 1.5;        // coordinate time between frames, in M
  double orbit_deg = 0.0;       // how far the camera circles during an animation
  int threads = std::max(1u, std::thread::hardware_concurrency());
  double mass_msun = 10.0;      // only used for the physical scales printed at start
  std::string output = "blackhole.png";
  bool test = false;
};

// ---------------------------------------------------------------------------
// Renderer

class Renderer {
 public:
  explicit Renderer(const Settings& s)
      : s_(s), gravity_(!s.flat), cam_lapse_(s.flat ? 1.0 : std::sqrt(1.0 - 2.0 * kM / s.distance)) {
    // Thin-disk temperature profile with a torque-free inner edge,
    // T ~ r^-3/4 (1 - sqrt(r_in/r))^1/4, which peaks at r = (49/36) r_in.
    double peak = std::pow(36.0 / 49.0, 0.75) * std::pow(1.0 / 7.0, 0.25);
    disk_scale_ = s.disk_temp / peak;
    disk_norm_ = 0.4 / luminance(bb_.radiance(s.disk_temp));
  }

  // Renders one frame into a linear HDR buffer.
  std::vector<Vec3> render(double azimuth_deg, double time, bool progress) const {
    const int w = s_.width, h = s_.height, aa = s_.aa;
    double elev = (90.0 - s_.incl_deg) * kPi / 180.0;
    double az = azimuth_deg * kPi / 180.0;
    Vec3 eye = Vec3(std::cos(elev) * std::sin(az), std::sin(elev), std::cos(elev) * std::cos(az)) *
               s_.distance;
    Vec3 fwd = normalize(eye * -1.0);
    Vec3 right = cross(fwd, Vec3(0, 1, 0));
    if (length(right) < 1e-9) right = Vec3(std::cos(az), 0, -std::sin(az));  // looking down the axis
    right = normalize(right);
    Vec3 up = cross(right, fwd);
    double tan_half = std::tan(0.5 * s_.fov_deg * kPi / 180.0);
    double aspect = double(w) / h;
    // Direction, in the camera's own frame, of the ray through image point (px, py).
    auto ray = [&](double px, double py) {
      double u = (2.0 * px / w - 1.0) * tan_half * aspect;
      double v = (1.0 - 2.0 * py / h) * tan_half;
      return normalize(fwd + right * u + up * v);
    };

    std::vector<Vec3> img(size_t(w) * h);
    std::atomic<int> next_row{0}, rows_done{0};
    auto worker = [&] {
      for (int y; (y = next_row++) < h;) {
        for (int x = 0; x < w; ++x) {
          Vec3 sum;
          for (int i = 0; i < aa * aa; ++i) {
            uint32_t rnd = hash3(x, y, i, 0xa11a5);
            double jx = (i % aa + to_unit(rnd)) / aa;
            double jy = (i / aa + to_unit(hash32(rnd))) / aa;
            sum += trace(eye, ray(x + jx, y + jy), time).light;
          }
          Vec3 pixel = sum / double(aa * aa);
          if (s_.sky == Sky::Stars) pixel += star_light(eye, ray, x + 0.5, y + 0.5, time);
          img[size_t(y) * w + x] = pixel;
        }
        int done = ++rows_done;
        if (progress && done * 20 / h != (done - 1) * 20 / h)
          std::fprintf(stderr, "\r  rendering %3d%%", 100 * done / h);
      }
    };
    std::vector<std::thread> pool;
    for (int t = 0; t < s_.threads; ++t) pool.emplace_back(worker);
    for (auto& t : pool) t.join();
    if (progress) std::fprintf(stderr, "\n");
    return img;
  }

 private:
  struct Sample {
    Vec3 light;             // disk and diffuse sky: everything except the stars
    double transmit = 1.0;  // fraction of the sky's light that gets through the disk
    bool escaped = false;
    Vec3 sky_dir;           // where the ray ends up on the sky
  };

  // Follows one photon back from the camera and collects the light it brings.
  Sample trace(const Vec3& eye, const Vec3& dir, double time) const {
    Photon ph(eye, observer_to_ray(eye, dir, gravity_), gravity_);
    Sample out;
    auto hit_disk = [&](const Vec3& p0, const Vec3&) {
      if (!s_.disk || p0.y * ph.pos.y > 0.0 || p0.y == ph.pos.y) return true;
      Vec3 hit = p0 + (ph.pos - p0) * (p0.y / (p0.y - ph.pos.y));  // crossing of the plane y = 0
      double r = std::hypot(hit.x, hit.z);
      if (r < s_.disk_in || r > s_.disk_out) return true;
      double alpha;
      out.light += disk_light(hit, r, ph.h.y, time, &alpha) * out.transmit;
      out.transmit *= 1.0 - alpha;
      return out.transmit > 1e-3;
    };
    if (propagate(ph, s_.step, kEscapeRadius, hit_disk) == Fate::Escaped) {
      out.escaped = true;
      out.sky_dir = normalize(ph.vel);
      if (s_.sky == Sky::Stars) out.light += milky_way(out.sky_dir) * out.transmit;
      if (s_.sky == Sky::Grid) out.light += grid(out.sky_dir) * out.transmit;
    }
    return out;
  }

  // Light from the disk element at `hit` (radius r), and its opacity.
  // hy is the component of the backwards ray's h along the disk axis.
  Vec3 disk_light(const Vec3& hit, double r, double hy, double time, double* alpha) const {
    double x = s_.disk_in / r;
    double temp = disk_scale_ * std::pow(x, 0.75) * std::pow(1.0 - std::sqrt(x), 0.25);

    // Gas on circular geodesics turns at the Kepler rate (as seen from far away).
    double omega = std::sqrt(kM / (r * r * r));
    double g = 1.0;  // nu_camera / nu_emitted
    switch (s_.shift) {
      case Shift::Full:
        // The emitter has u^t = 1/sqrt(1 - 3M/r), and the photon's angular
        // momentum per unit energy about the axis is b_y = -h_y / |v_inf|
        // (minus because we trace it backwards; |v_inf| = cam_lapse_). So
        // g = sqrt(1 - 3M/r) / (cam_lapse * (1 - Omega b_y)): it combines the
        // gravitational redshift, time dilation and the Doppler shift.
        g = std::sqrt(1.0 - 3.0 * kM / r) /
            (cam_lapse_ * std::max(1e-6, 1.0 + omega * hy / cam_lapse_));
        break;
      case Shift::Gravity:  // as if the gas were at rest
        g = std::sqrt(1.0 - 2.0 * kM / r) / cam_lapse_;
        break;
      case Shift::None:
        break;
    }

    // Turbulent clumps, stretched along the orbit and carried round with the gas.
    double phi = std::atan2(-hit.z, hit.x) - omega * time;
    double clump =
        fbm(Vec3(std::cos(phi) * 2.5, std::sin(phi) * 2.5, std::log(r) * 9.0), 5);
    double density = smoothstep(0.3, 0.7, clump);
    double edge = smoothstep(s_.disk_out, s_.disk_out - 0.3 * (s_.disk_out - s_.disk_in), r);
    *alpha = edge * (0.8 + 0.2 * density);
    return bb_.radiance(g * temp) * (disk_norm_ * *alpha * (0.35 + 1.1 * density));
  }

  // Stars are points, so lensing cannot smear them: it moves them and makes
  // them brighter or fainter. Two extra rays a fraction of a pixel away give
  // the Jacobian J of the map from image to sky, so each star is drawn as a
  // sharp dot where it lands, its flux multiplied by the magnification
  // |det J_empty| / |det J|, relative to the same camera in empty space.
  template <class RayFn>
  Vec3 star_light(const Vec3& eye, const RayFn& ray, double px, double py, double time) const {
    const double d = 0.25;  // pixels
    Sample c = trace(eye, ray(px, py), time);
    if (!c.escaped || c.transmit < 1e-3) return {};
    Sample sx = trace(eye, ray(px + d, py), time);
    Sample sy = trace(eye, ray(px, py + d), time);
    if (!sx.escaped || !sy.escaped) return {};
    Vec3 d0 = ray(px, py);
    double empty_sr = length(cross(ray(px + d, py) - d0, ray(px, py + d) - d0)) / (d * d);
    return stars(c.sky_dir, (sx.sky_dir - c.sky_dir) / d, (sy.sky_dir - c.sky_dir) / d, empty_sr) *
           c.transmit;
  }

  // Stars sit in the cells of a 3-D lattice that straddle a sphere, at most
  // one per cell, which gives an even starfield with no seams or poles.
  // dx, dy: change of the sky direction per pixel; empty_sr: solid angle of
  // the pixel if there were no black hole.
  Vec3 stars(const Vec3& d, const Vec3& dx, const Vec3& dy, double empty_sr) const {
    constexpr double kCells = 180.0, kDensity = 0.15, kSigma = 0.5;  // kSigma in pixels
    Vec3 e1 = normalize(cross(d, std::abs(d.y) < 0.9 ? Vec3(0, 1, 0) : Vec3(1, 0, 0)));
    Vec3 e2 = cross(d, e1);
    double a = dot(dx, e1), b = dot(dy, e1), c = dot(dx, e2), e = dot(dy, e2);
    double det = a * e - b * c;
    if (!(std::abs(det) > 1e-30)) return {};
    double gain = std::min(empty_sr / std::abs(det), 30.0);  // capped near the caustics
    // Look far enough out on the sky to find every star within ~3 pixels.
    double reach = 3.0 * std::max(std::hypot(a, c), std::hypot(b, e)) * kCells;
    int n = std::clamp(static_cast<int>(std::ceil(reach)), 1, 3);

    Vec3 q = d * kCells;
    int cx = static_cast<int>(std::floor(q.x));
    int cy = static_cast<int>(std::floor(q.y));
    int cz = static_cast<int>(std::floor(q.z));
    Vec3 sum;
    for (int i = -n; i <= n; ++i)
      for (int j = -n; j <= n; ++j)
        for (int k = -n; k <= n; ++k) {
          uint32_t h = hash3(cx + i, cy + j, cz + k, 0x57a25);
          if (to_unit(h) > kDensity) continue;
          double ox = to_unit(h = hash32(h)), oy = to_unit(h = hash32(h)), oz = to_unit(h = hash32(h));
          Vec3 p(cx + i + ox, cy + j + oy, cz + k + oz);
          double len = length(p);
          if (std::abs(len - kCells) > 0.5) continue;
          Vec3 off = p / len - d;
          double s1 = dot(off, e1), s2 = dot(off, e2);
          double px = (e * s1 - b * s2) / det, py = (a * s2 - c * s1) / det;  // J^-1: back to pixels
          double r2 = px * px + py * py;
          if (r2 > 9.0) continue;
          double mag = to_unit(h = hash32(h));
          double flux = 0.03 + 15.0 * std::pow(mag, 24.0);  // many faint stars, a few bright ones
          double temp = 3000.0 + 22000.0 * std::pow(to_unit(hash32(h)), 3.0);
          sum += mix(Vec3(1, 1, 1), bb_.tint(temp), 0.7) *
                 (flux * gain * std::exp(-r2 / (2.0 * kSigma * kSigma)) / (2.0 * kPi * kSigma * kSigma));
        }
    return sum;
  }

  Vec3 milky_way(const Vec3& d) const {
    const Vec3 pole = normalize(Vec3(0.6, 0.78, 0.15));
    double lat = dot(d, pole);
    double cloud = fbm(d * 3.0 + Vec3(7.1, 3.3, 1.7), 5);
    double width = 0.1 + 0.15 * cloud;
    double band = std::exp(-lat * lat / (2.0 * width * width));
    double dust = smoothstep(0.45, 0.7, fbm(d * 8.0 + Vec3(1.3, 9.2, 4.4), 5));
    return mix(Vec3(0.55, 0.62, 0.85), Vec3(1.0, 0.82, 0.62), cloud) *
           (0.05 * band * (1.0 - 0.8 * dust));
  }

  // 10-degree checkerboard painted on the sky: makes the lensing easy to read.
  Vec3 grid(const Vec3& d) const {
    const double cell = kPi / 18.0;
    double theta = std::acos(std::clamp(d.y, -1.0, 1.0));
    double phi = std::atan2(d.z, d.x) + kPi;
    bool odd = (static_cast<int>(theta / cell) + static_cast<int>(phi / cell)) & 1;
    return odd ? Vec3(0.9, 0.45, 0.1) * 0.3 : Vec3(0.1, 0.35, 0.8) * 0.3;
  }

  const Settings& s_;
  Blackbody bb_;
  bool gravity_;
  double cam_lapse_;    // sqrt(1 - 2M/r_camera), or 1 in flat space
  double disk_scale_;   // turns the disk temperature profile into kelvin
  double disk_norm_;    // gives the hottest disk element, at rest, luminance 0.4
};

// ---------------------------------------------------------------------------
// Post-processing: bloom, tone mapping

void box_blur_line(Vec3* data, int n, int stride, int r, std::vector<Vec3>& line) {
  line.assign(n, Vec3());
  for (int i = 0; i < n; ++i) line[i] = data[size_t(i) * stride];
  Vec3 acc;
  for (int i = -r; i <= r; ++i) acc += line[std::clamp(i, 0, n - 1)];
  double inv = 1.0 / (2 * r + 1);
  for (int i = 0; i < n; ++i) {
    data[size_t(i) * stride] = acc * inv;
    acc += line[std::min(i + r + 1, n - 1)] - line[std::max(i - r, 0)];
  }
}

// Three box blurs in a row approximate a Gaussian.
void blur(std::vector<Vec3>& img, int w, int h, int r) {
  std::vector<Vec3> line;
  for (int pass = 0; pass < 3; ++pass) {
    for (int y = 0; y < h; ++y) box_blur_line(&img[size_t(y) * w], w, 1, r, line);
    for (int x = 0; x < w; ++x) box_blur_line(&img[x], h, w, r, line);
  }
}

// Bright light scatters a little inside any real lens or eye.
void add_bloom(std::vector<Vec3>& img, int w, int h) {
  const double threshold = 0.9;
  std::vector<Vec3> bright(img.size());
  for (size_t i = 0; i < img.size(); ++i) {
    double l = luminance(img[i]);
    if (l > threshold) bright[i] = img[i] * ((l - threshold) / l);
  }
  const std::pair<double, double> layers[] = {{0.003, 0.20}, {0.015, 0.10}, {0.05, 0.06}};
  for (auto [radius, weight] : layers) {
    std::vector<Vec3> b = bright;
    blur(b, w, h, std::max(1, static_cast<int>(radius * h)));
    for (size_t i = 0; i < img.size(); ++i) img[i] += b[i] * weight;
  }
}

double aces(double x) {  // filmic curve (Narkowicz' fit of ACES)
  return std::clamp(x * (2.51 * x + 0.03) / (x * (2.43 * x + 0.59) + 0.14), 0.0, 1.0);
}

double srgb_encode(double c) {
  return c <= 0.0031308 ? 12.92 * c : 1.055 * std::pow(c, 1.0 / 2.4) - 0.055;
}

std::vector<uint8_t> tone_map(const std::vector<Vec3>& img, double exposure) {
  std::vector<uint8_t> out;
  out.reserve(img.size() * 3);
  for (const Vec3& c : img)
    for (double v : {c.x, c.y, c.z})
      out.push_back(static_cast<uint8_t>(std::lround(255.0 * srgb_encode(aces(v * exposure)))));
  return out;
}

// ---------------------------------------------------------------------------
// Image files

bool write_ppm(const std::string& path, int w, int h, const std::vector<uint8_t>& rgb) {
  std::ofstream f(path, std::ios::binary);
  f << "P6\n" << w << ' ' << h << "\n255\n";
  f.write(reinterpret_cast<const char*>(rgb.data()), static_cast<std::streamsize>(rgb.size()));
  return bool(f);
}

uint32_t crc32(const std::vector<uint8_t>& data) {
  static const auto table = [] {
    std::array<uint32_t, 256> t{};
    for (uint32_t i = 0; i < 256; ++i) {
      uint32_t c = i;
      for (int k = 0; k < 8; ++k) c = (c & 1) ? 0xedb88320U ^ (c >> 1) : c >> 1;
      t[i] = c;
    }
    return t;
  }();
  uint32_t crc = 0xffffffffU;
  for (uint8_t b : data) crc = table[(crc ^ b) & 0xff] ^ (crc >> 8);
  return ~crc;
}

void put_u32(std::vector<uint8_t>& v, uint32_t x) {
  for (int s = 24; s >= 0; s -= 8) v.push_back(static_cast<uint8_t>(x >> s));
}

// A minimal PNG writer. The pixels go into uncompressed ("stored") deflate
// blocks: the files are bigger than they need to be, but any viewer opens
// them and there is no zlib dependency.
bool write_png(const std::string& path, int w, int h, const std::vector<uint8_t>& rgb) {
  std::vector<uint8_t> raw;
  raw.reserve(size_t(h) * (3 * size_t(w) + 1));
  for (int y = 0; y < h; ++y) {
    raw.push_back(0);  // no filter
    auto row = rgb.begin() + size_t(y) * 3 * w;
    raw.insert(raw.end(), row, row + 3 * size_t(w));
  }

  std::vector<uint8_t> zlib = {0x78, 0x01};
  size_t pos = 0;
  do {
    size_t n = std::min<size_t>(65535, raw.size() - pos);
    zlib.push_back(pos + n == raw.size() ? 1 : 0);
    zlib.push_back(static_cast<uint8_t>(n));
    zlib.push_back(static_cast<uint8_t>(n >> 8));
    zlib.push_back(static_cast<uint8_t>(~n));
    zlib.push_back(static_cast<uint8_t>(~n >> 8));
    zlib.insert(zlib.end(), raw.begin() + pos, raw.begin() + pos + n);
    pos += n;
  } while (pos < raw.size());
  uint32_t a = 1, b = 0;
  for (uint8_t c : raw) {
    a = (a + c) % 65521;
    b = (b + a) % 65521;
  }
  put_u32(zlib, (b << 16) | a);

  std::vector<uint8_t> ihdr;
  put_u32(ihdr, static_cast<uint32_t>(w));
  put_u32(ihdr, static_cast<uint32_t>(h));
  ihdr.insert(ihdr.end(), {8, 2, 0, 0, 0});  // 8-bit RGB

  std::ofstream f(path, std::ios::binary);
  const uint8_t signature[8] = {0x89, 'P', 'N', 'G', '\r', '\n', 0x1a, '\n'};
  f.write(reinterpret_cast<const char*>(signature), 8);
  auto chunk = [&](const char* type, const std::vector<uint8_t>& data) {
    std::vector<uint8_t> body(type, type + 4);
    body.insert(body.end(), data.begin(), data.end());
    std::vector<uint8_t> len, crc;
    put_u32(len, static_cast<uint32_t>(data.size()));
    put_u32(crc, crc32(body));
    f.write(reinterpret_cast<const char*>(len.data()), 4);
    f.write(reinterpret_cast<const char*>(body.data()), static_cast<std::streamsize>(body.size()));
    f.write(reinterpret_cast<const char*>(crc.data()), 4);
  };
  chunk("IHDR", ihdr);
  chunk("IDAT", zlib);
  chunk("IEND", {});
  return bool(f);
}

bool ends_with(const std::string& s, const std::string& suffix) {
  return s.size() >= suffix.size() &&
         std::equal(suffix.rbegin(), suffix.rend(), s.rbegin(),
                    [](unsigned char a, unsigned char b) { return std::tolower(a) == std::tolower(b); });
}

bool write_image(const std::string& path, int w, int h, const std::vector<uint8_t>& rgb) {
  return ends_with(path, ".ppm") ? write_ppm(path, w, h, rgb) : write_png(path, w, h, rgb);
}

// "out.png", 7 -> "out_0007.png"
std::string frame_path(const std::string& path, int frame) {
  char num[16];
  std::snprintf(num, sizeof num, "_%04d", frame);
  size_t dot_pos = path.find_last_of('.');
  size_t slash = path.find_last_of("/\\");
  if (dot_pos == std::string::npos || (slash != std::string::npos && dot_pos < slash))
    return path + num;
  return path.substr(0, dot_pos) + num + path.substr(dot_pos);
}

// ---------------------------------------------------------------------------
// Self-test: the integrator against exact results of general relativity

// Exact bending angle of a light ray whose closest approach is r0:
//   alpha = 2 * Integral_0^u0 du / sqrt(1/b^2 - u^2 + 2Mu^3) - pi,   u0 = 1/r0.
// Substituting u = u0 (1 - t^2) and dividing out the root at u0 leaves a
// smooth integrand, 2 sqrt(u0 / Q(u)) with Q(u) = u + u0 - 2M(u^2 + u u0 + u0^2).
double exact_deflection(double r0) {
  double u0 = 1.0 / r0;
  auto f = [u0](double t) {
    double u = u0 * (1.0 - t * t);
    return 2.0 * std::sqrt(u0 / (u + u0 - 2.0 * kM * (u * u + u * u0 + u0 * u0)));
  };
  const int n = 4000;  // Simpson's rule
  double sum = f(0.0) + f(1.0);
  for (int i = 1; i < n; ++i) sum += f(double(i) / n) * (i % 2 ? 4.0 : 2.0);
  return 2.0 * sum / (3.0 * n) - kPi;
}

int run_self_test(double step) {
  int failures = 0;
  auto verdict = [&](double err, double tol) {
    bool ok = err < tol;
    failures += !ok;
    return ok ? "ok" : "FAIL";
  };
  auto no_disk = [](const Vec3&, const Vec3&) { return true; };

  std::printf("Self-test, step = %.3g r\n\n", step);

  // 1. Edge of the shadow seen by a static observer at radius r (Synge 1966):
  //    sin(alpha) = 3 sqrt(3) M sqrt(1 - 2M/r) / r.
  //    This checks both the geodesics and the camera model.
  std::printf("Angular radius of the shadow (static observer):\n");
  for (double r : {4.0, 6.0, 10.0, 30.0, 100.0}) {
    auto captured = [&](double alpha) {
      Vec3 pos(0, 0, r), dir(std::sin(alpha), 0, -std::cos(alpha));
      Photon ph(pos, observer_to_ray(pos, dir, true), true);
      return propagate(ph, step, kEscapeRadius, no_disk) == Fate::Captured;
    };
    double lo = 0.0, hi = 0.5 * kPi;
    for (int i = 0; i < 50; ++i) {
      double mid = 0.5 * (lo + hi);
      (captured(mid) ? lo : hi) = mid;
    }
    double exact = std::asin(kCriticalB * std::sqrt(1.0 - 2.0 * kM / r) / r);
    double err = std::abs(lo - exact) / exact;
    std::printf("  r = %5.1f M   traced %9.5f deg   exact %9.5f deg   rel. error %.1e  %s\n", r,
                lo * 180 / kPi, exact * 180 / kPi, err, verdict(err, 1e-4));
  }

  // 2. Bending of light passing the hole, against the exact integral.
  std::printf("\nDeflection of a ray with closest approach r0:\n");
  for (double r0 : {3.2, 3.5, 4.0, 6.0, 10.0, 30.0, 100.0, 1000.0}) {
    double b = r0 / std::sqrt(1.0 - 2.0 * kM / r0);
    Photon ph(Vec3(b, 0, -1e6), Vec3(0, 0, 1), true);
    double turned = 0.0;
    auto count_turn = [&](const Vec3&, const Vec3& v0) {
      turned += std::atan2(length(cross(v0, ph.vel)), dot(v0, ph.vel));
      return true;
    };
    Fate fate = propagate(ph, step, 1e6, count_turn);
    double exact = exact_deflection(r0);
    double err = fate == Fate::Escaped ? std::abs(turned - exact) / exact : 1.0;
    std::printf("  r0 = %6.1f M  b = %8.3f M   traced %10.6f rad   exact %10.6f rad   "
                "rel. error %.1e  %s\n",
                r0, b, turned, exact, err, verdict(err, 1e-4));
  }
  {
    double b = 1000.0 / std::sqrt(1.0 - 2.0 * kM / 1000.0);
    std::printf("  (weak field: Einstein's 4M/b = %.6f rad for b = %.1f M)\n", 4.0 * kM / b, b);
  }

  // 3. Conserved quantities on a ray that whirls round the photon sphere.
  std::printf("\nConservation along a ray that orbits the photon sphere:\n");
  {
    double r = 30.0;
    double edge = std::asin(kCriticalB * std::sqrt(1.0 - 2.0 * kM / r) / r);
    Vec3 pos(0, 0, r), dir(std::sin(edge * (1 + 1e-9)), 0, -std::cos(edge * (1 + 1e-9)));
    Photon ph(pos, observer_to_ray(pos, dir, true), true);
    double h0 = length(ph.h);
    auto invariant = [&] {
      double rr = length(ph.pos);
      return dot(ph.vel, ph.vel) - 2.0 * kM * dot(ph.h, ph.h) / (rr * rr * rr);
    };
    double c0 = invariant();
    double phi = 0.0;
    auto track = [&](const Vec3& p0, const Vec3&) {
      phi += std::atan2(length(cross(p0, ph.pos)), dot(p0, ph.pos));
      return true;
    };
    propagate(ph, step, kEscapeRadius, track);
    double dh = std::abs(length(cross(ph.pos, ph.vel)) - h0) / h0;
    double dc = std::abs(invariant() - c0) / c0;
    std::printf("  swept %.1f turns; drift of |h| %.1e  %s; of |v|^2 - 2M h^2/r^3 %.1e  %s\n",
                phi / (2 * kPi), dh, verdict(dh, 1e-6), dc, verdict(dc, 1e-4));
  }

  std::printf("\n%s\n", failures ? "SOME CHECKS FAILED" : "All checks passed.");
  return failures ? 1 : 0;
}

// ---------------------------------------------------------------------------
// Command line

void print_usage() {
  std::printf(
      "Ray-traced Schwarzschild black hole with an accretion disk.\n"
      "Lengths are in units of the mass M (G = c = 1): horizon 2, photon sphere 3, ISCO 6.\n"
      "\n"
      "usage: blackhole [options]\n"
      "\n"
      "  -o FILE            output image, .png or .ppm (default blackhole.png)\n"
      "  --width N          image width in pixels (1280)\n"
      "  --height N         image height in pixels (720)\n"
      "  --fov DEG          vertical field of view (48)\n"
      "  --dist R           camera distance from the hole, in M (30)\n"
      "  --incl DEG         line of sight to disk axis: 90 = edge-on, 0 = face-on (84)\n"
      "  --azimuth DEG      camera position around the axis (0)\n"
      "  --disk-in R        inner edge of the disk, > 3 (6, the ISCO)\n"
      "  --disk-out R       outer edge of the disk (20)\n"
      "  --temp K           peak disk temperature in the gas frame (4500)\n"
      "  --no-disk          no accretion disk: just the hole and the sky\n"
      "  --sky S            stars | grid | none (stars)\n"
      "  --shift S          full | gravity | none: which frequency shifts to apply (full)\n"
      "  --flat             light travels in straight lines (no lensing), to compare\n"
      "  --exposure X       brightness multiplier before tone mapping (1)\n"
      "  --aa N             N x N samples per pixel (2)\n"
      "  --no-bloom         no glow around bright light\n"
      "  --step F           integration step as a fraction of r (0.03)\n"
      "  --frames N         render an animation of N frames (1)\n"
      "  --frame-dt T       time between frames, in M (1.5)\n"
      "  --orbit DEG        camera turns this far around the hole during the animation (0)\n"
      "  --threads N        worker threads (all cores)\n"
      "  --mass MSUN        mass in solar masses, for the physical scales printed (10)\n"
      "  --test             check the integrator against exact GR results, then exit\n"
      "  --help             this text\n");
}

// Returns 0 to go on, otherwise the process exit code + 1.
int parse_args(int argc, char** argv, Settings& s) {
  auto fail = [](const std::string& msg) {
    std::fprintf(stderr, "blackhole: %s (see --help)\n", msg.c_str());
    return 2;
  };
  for (int i = 1; i < argc; ++i) {
    std::string a = argv[i];
    auto value = [&]() -> const char* { return i + 1 < argc ? argv[++i] : nullptr; };
    auto num = [&](double& out) {
      const char* v = value();
      char* end = nullptr;
      if (v) out = std::strtod(v, &end);
      return v && end != v && *end == '\0';
    };
    auto integer = [&](int& out) {
      double d;
      if (!num(d) || d != std::floor(d) || std::abs(d) > 1e6) return false;
      out = static_cast<int>(d);
      return true;
    };
    bool ok = true;
    if (a == "--help" || a == "-h") {
      print_usage();
      return 1;
    } else if (a == "-o" || a == "--output") {
      const char* v = value();
      ok = v != nullptr;
      if (ok) s.output = v;
    } else if (a == "--width") ok = integer(s.width);
    else if (a == "--height") ok = integer(s.height);
    else if (a == "--fov") ok = num(s.fov_deg);
    else if (a == "--dist") ok = num(s.distance);
    else if (a == "--incl") ok = num(s.incl_deg);
    else if (a == "--azimuth") ok = num(s.azimuth_deg);
    else if (a == "--disk-in") ok = num(s.disk_in);
    else if (a == "--disk-out") ok = num(s.disk_out);
    else if (a == "--temp") ok = num(s.disk_temp);
    else if (a == "--no-disk") s.disk = false;
    else if (a == "--exposure") ok = num(s.exposure);
    else if (a == "--aa") ok = integer(s.aa);
    else if (a == "--no-bloom") s.bloom = false;
    else if (a == "--step") ok = num(s.step);
    else if (a == "--flat") s.flat = true;
    else if (a == "--frames") ok = integer(s.frames);
    else if (a == "--frame-dt") ok = num(s.frame_dt);
    else if (a == "--orbit") ok = num(s.orbit_deg);
    else if (a == "--threads") ok = integer(s.threads);
    else if (a == "--mass") ok = num(s.mass_msun);
    else if (a == "--test") s.test = true;
    else if (a == "--sky") {
      const char* v = value();
      std::string m = v ? v : "";
      ok = m == "stars" || m == "grid" || m == "none";
      s.sky = m == "grid" ? Sky::Grid : m == "none" ? Sky::None : Sky::Stars;
    } else if (a == "--shift") {
      const char* v = value();
      std::string m = v ? v : "";
      ok = m == "full" || m == "gravity" || m == "none";
      s.shift = m == "gravity" ? Shift::Gravity : m == "none" ? Shift::None : Shift::Full;
    } else {
      return fail("unknown option " + a);
    }
    if (!ok) return fail("bad or missing value for " + a);
  }

  if (s.width < 1 || s.height < 1 || s.width > 16384 || s.height > 16384)
    return fail("image size must be between 1 and 16384 pixels");
  if (!(s.fov_deg > 0 && s.fov_deg < 179)) return fail("--fov must be between 0 and 179");
  if (!(s.distance > 2.05 * kM)) return fail("--dist must be outside the horizon (> 2.05)");
  if (!(s.incl_deg >= 0 && s.incl_deg <= 180)) return fail("--incl must be between 0 and 180");
  if (!(s.disk_in > kPhotonSphere)) return fail("--disk-in must be > 3: no circular orbits inside");
  if (!(s.disk_out > s.disk_in)) return fail("--disk-out must be larger than --disk-in");
  if (!(s.disk_temp > 0)) return fail("--temp must be positive");
  if (!(s.exposure > 0)) return fail("--exposure must be positive");
  if (s.aa < 1 || s.aa > 16) return fail("--aa must be between 1 and 16");
  if (!(s.step > 0 && s.step <= 0.2)) return fail("--step must be in (0, 0.2]");
  if (s.frames < 1) return fail("--frames must be at least 1");
  if (s.threads < 1) return fail("--threads must be at least 1");
  if (!(s.mass_msun > 0)) return fail("--mass must be positive");
  if (s.disk && s.disk_in < kIsco)
    std::fprintf(stderr, "note: orbits inside r = 6 (the ISCO) are unstable; gas there would plunge\n");
  return 0;
}

void print_scales(const Settings& s) {
  double km = s.mass_msun * kMetresPerMsun / 1000.0;
  double sec = s.mass_msun * kSecondsPerMsun;
  auto period = [&](double r) { return 2.0 * kPi * std::sqrt(r * r * r / kM) * sec; };
  auto speed = [](double r) { return std::sqrt(kM / (r - 2.0 * kM)); };  // seen by a static observer
  double shadow = std::asin(std::min(1.0, kCriticalB * std::sqrt(1.0 - 2.0 * kM / s.distance) /
                                              s.distance));
  if (s.distance < kPhotonSphere) shadow = kPi - shadow;
  std::printf("Schwarzschild black hole of %g solar masses (M = GM/c^2 = %.4g km)\n", s.mass_msun, km);
  std::printf("  event horizon    r = 2M      = %10.4g km\n", 2 * km);
  std::printf("  photon sphere    r = 3M      = %10.4g km\n", 3 * km);
  std::printf("  ISCO             r = 6M      = %10.4g km   orbit %.3g s at %.2fc\n", 6 * km,
              period(kIsco), speed(kIsco));
  std::printf("  shadow edge      b = 3sqrt3 M = %9.4g km\n", kCriticalB * km);
  std::printf("  camera at        r = %gM = %.4g km, sees the shadow %.2f deg across\n",
              s.distance, s.distance * km, 2 * shadow * 180 / kPi);
  if (s.disk)
    std::printf("  disk from r = %gM to %gM; gas at the inner edge orbits every %.3g s\n", s.disk_in,
                s.disk_out, period(s.disk_in));
}

}  // namespace

int main(int argc, char** argv) {
  Settings s;
  if (int code = parse_args(argc, argv, s)) return code - 1;
  if (s.test) return run_self_test(s.step);

  print_scales(s);
  std::fflush(stdout);
  Renderer renderer(s);
  for (int f = 0; f < s.frames; ++f) {
    auto t0 = std::chrono::steady_clock::now();
    double azimuth = s.azimuth_deg + (s.frames > 1 ? s.orbit_deg * f / s.frames : 0.0);
    std::vector<Vec3> hdr = renderer.render(azimuth, f * s.frame_dt, true);
    if (s.bloom) add_bloom(hdr, s.width, s.height);
    std::string path = s.frames > 1 ? frame_path(s.output, f) : s.output;
    if (!write_image(path, s.width, s.height, tone_map(hdr, s.exposure))) {
      std::fprintf(stderr, "blackhole: could not write %s\n", path.c_str());
      return 1;
    }
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    std::printf("wrote %s (%dx%d, %.1f s)\n", path.c_str(), s.width, s.height, secs);
  }
  return 0;
}
