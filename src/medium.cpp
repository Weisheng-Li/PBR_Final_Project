#include <dirt/medium.h>
#include <dirt/parser.h>
#include <dirt/onb.h>
#include <dirt/ray.h>
#include <dirt/scene.h>
#include <dirt/sampler.h>

HenyeyGreenstein::HenyeyGreenstein(const json &j)
{
  g = j.value("g", g);
}

// HenyeyGreenstein Method Definitions
float HenyeyGreenstein::p(const Vec3f &wo, const Vec3f &wi) const
{
  float cosTheta = dot(normalize(wo), normalize(wi));
  float denom = 1.0 + g * g + 2.0 * g * cosTheta;
  return INV_FOURPI * (1.0 - g * g) / (denom * std::sqrt(denom));
}

float HenyeyGreenstein::sample(const Vec3f &wo, Vec3f &wi, const Vec2f &u) const
{
  float cosTheta;
  if (std::abs(g) < 1e-3)
  {
    cosTheta = 1.0 - 2.0 * u.x;
  }
  else
  {
    float sqrTerm = (1.0 - g * g) / (1.0 + g - 2.0 * g * u.x);
    cosTheta = -(1.0 + g * g - sqrTerm * sqrTerm) / (2.0 * g);
  }

  // Compute direction _wi_ for Henyey--Greenstein sample
  float sinTheta = std::sqrt(std::max(0.0f, 1.0f - cosTheta * cosTheta));
  float phi = 2 * M_PI * u.y;

  ONB<float> onb;
  onb.build_from_w(wo);
  wi = normalize(onb.toWorld(SphericalDirection(sinTheta, cosTheta, phi)));
  return p(wo, wi);
}

// !!! Be aware that DIRT sample wi given wo, but SGGX paper sample wo given wi

// Specular SGGX phase function, serves as both the evaluation and PDF
float SGGXPhase::p(const Vec3f &wo, const Vec3f &wi) const {
  Vec3f half_vec = normalize(normalize(wo) + normalize(wi));
  float result = NDF(half_vec) / (4 * projection(wo));
  return NDF(half_vec) / (4 * projection(wo));
}

float SGGXPhase::sample(const Vec3f &wo, Vec3f &wi, const Vec2f &u) const {
  assert(abs(length(wo) - 1.f) < 0.001);

  // Project S into the basis
  Mat44f S_local = Mat44f(0.0f);

  ONB<float> onb;
  onb.build_from_w(wo);
  Vec4f wi4 = {onb.w, 0.f};
  Vec4f wj4 = {onb.v, 0.f};
  Vec4f wk4 = {onb.u, 0.f};

  // (row, col)
  S_local(0, 0) = dot(wk4, S * wk4);
  S_local(0, 1) = dot(wk4, S * wj4);
  S_local(0, 2) = dot(wk4, S * wi4);
  S_local(1, 0) = dot(wk4, S * wj4);
  S_local(1, 1) = dot(wj4, S * wj4);
  S_local(1, 2) = dot(wj4, S * wi4);
  S_local(2, 0) = dot(wk4, S * wi4);
  S_local(2, 1) = dot(wj4, S * wi4);
  S_local(2, 2) = dot(wi4, S * wi4);

  Vec3f Mk = {
    sqrt(determinant3x3(S_local) / (S_local(1,1)*S_local(2,2) - S_local(1,2) * S_local(1,2))),
    0,
    0,
  };

  Vec3f Mj = {
    -(S_local(2,0)*S_local(2,1) - S_local(0,1)*S_local(2,2)) / 
    sqrt(S_local(1,1)*S_local(2,2) - S_local(1,2)*S_local(1,2)),
    sqrt(S_local(1,1)*S_local(2,2) - S_local(1,2)*S_local(1,2)),
    0,
  };
  Mj *= 1 / sqrt(S_local(2,2));

  Vec3f Mi = {
    S_local(0,2), S_local(1,2), S_local(2,2)
  };
  Mi *= 1 / sqrt(S_local(2,2));

  // Sample on sphere (maybe?)
  float uu = sqrt(u.x) * cos(2 * M_PI * u.y);
  float vv =  sqrt(u.x) * sin(2 * M_PI * u.y);
  float ww = sqrt(1 - uu*uu - vv*vv);

  // compute the normal
  Vec3f M_sum = uu * Mk + vv * Mj + ww * Mi;
  Vec3f n = normalize(onb.toWorld(M_sum));

  // now simply reflect around the normal (since micro-phase function is dirac delta)
  wi = -wo + 2 * n * dot(wo, n);

  return p(wo, wi);
}

float SGGXPhase::NDF(const Vec3f &wh) const {
  Vec3f wh_n = normalize(wh);
  Vec4f wh4 = {wh_n, 0.f};
  float wSw = dot(wh4, inverse3x3(S) * wh4);
  float denom = M_PI * sqrt(determinant3x3(S)) * wSw * wSw;
  return 1.f / denom;
}

float SGGXPhase::projection(const Vec3f &wi) const {
  Vec3f wi_n = normalize(wi);
  Vec4f wi4 = {wi_n, 0.f};
  return sqrt(dot(wi4, S * wi4));
}

HomogeneousMedium::HomogeneousMedium(const json &j)
{
  phase = parsePhase(j.at("phase"));
  sigma_a = j.value("sigma_a", sigma_a);
  sigma_s = j.value("sigma_s", sigma_s);
  sigma_t = sigma_s + sigma_a;
}

float HomogeneousMedium::Tr(const Ray3f &ray_, Sampler &sampler) const
{
  Ray3f ray = ray_.normalizeRay();
  return std::exp(-sigma_t * (ray.maxt - ray.mint));
}

float HomogeneousMedium::Sample(const Ray3f &ray_, Sampler &sampler, MediumInteraction &mi) const
{
  Ray3f ray = ray_.normalizeRay();
  float dist = -std::log(1.0f - sampler.next1D()) / sigma_t;
  float t = std::min(dist, ray.maxt);
  bool sampledMedium = t < ray.maxt;
  if (sampledMedium)
    mi = MediumInteraction(ray(t), -ray.d, this);
  return sampledMedium ? sigma_s / sigma_t : 1.0f;
}

float HomogeneousMedium::density(const Vec3f &p) const
{
  return sigma_t;
}

PerlinMedium::PerlinMedium(const json &j)
{
  phase = parsePhase(j.at("phase"));
  sigma_a = j.value("sigma_a", sigma_a);
  sigma_s = j.value("sigma_s", sigma_s);
  sigma_t = sigma_s + sigma_a;

  spatialScale = j.value("spatial_scale", spatialScale);
  densityScale = abs(j.value("density_scale", densityScale));
  densityOffset = j.value("density_offset", densityOffset);

  assert(densityScale + densityOffset > 0.0f);

  invMaxDensity = 1.0f / (sigma_t * (densityScale + densityOffset));
}

float PerlinMedium::Tr(const Ray3f &ray_, Sampler &sampler) const
{
  Ray3f ray = ray_.normalizeRay();
  float Tr = 1;
  float t = ray.mint;
  int i = 0;
  while (true) {
    t -= std::log(1.0 - sampler.next1D()) * invMaxDensity;
    if (t * invMaxDensity >= ray.maxt) break;
    Tr *= 1.0 - std::max((float)0,  density(ray(t)) * invMaxDensity);
    if (Tr < Epsilon) break; // Guard against infinite looping. Should we set Tr = 0?
    i++;
    if (i % 1000000 == 0) warning("[PerlinMedium::Tr] %d-th iteration. Tr = %f\n", i, Tr);
  }

  return Tr;
}

float PerlinMedium::Sample(const Ray3f &ray_, Sampler &sampler, MediumInteraction &mi) const
{
  Ray3f ray = ray_.normalizeRay();
  float t = ray.mint;
  int i = 0;
  while (true)
  {
    t -= std::log(1.0 - sampler.next1D()) * invMaxDensity;
    if (ray.maxt <= t) break;
    if (sampler.next1D() < density(ray(t)) * invMaxDensity)
    {
      mi = MediumInteraction(ray(t), -ray.d, this);
      return sigma_s / sigma_t;
    }
    i++;
    if (i % 1000000 == 0) warning("[PerlinMedium::Sample] %d-th iteration. t = %f\n", i, t);
  }
  return 1.0f;
}

float PerlinMedium::density(const Vec3f &p) const
{
  Vec3f pScaled(p.x * spatialScale.x, p.y * spatialScale.y, p.z * spatialScale.z);
  return sigma_t * std::max(0.0f, densityScale * perlin.noise(pScaled) + densityOffset);
}

HomogeneousAnisotropicMedium::HomogeneousAnisotropicMedium(const json &j)
{
  rho = j.value("rho", rho);
  albedo = j.value("albedo", albedo);

  // identity matrix without the last diagonal 1
  S = Mat44f(1);
  S(3,3) = 0;

  phase = make_shared<SGGXPhase>(S);
}

float HomogeneousAnisotropicMedium::Tr(const Ray3f &ray_, Sampler &sampler) const
{
  Ray3f ray = ray_.normalizeRay();
  return std::exp(-sigma_t(-ray.d) * (ray.maxt - ray.mint));
}

float HomogeneousAnisotropicMedium::Sample(const Ray3f &ray_, Sampler &sampler, MediumInteraction &mi) const
{
  Ray3f ray = ray_.normalizeRay();
  float dist = -std::log(1.0f - sampler.next1D()) / sigma_t(-ray.d);
  float t = std::min(dist, ray.maxt);
  bool sampledMedium = t < ray.maxt;
  if (sampledMedium)
    mi = MediumInteraction(ray(t), -ray.d, this);
  return sampledMedium ? sigma_s(-ray.d) / sigma_t(-ray.d) : 1.0f;
}

// useless. In homogeneous medium density doesn't depend on p.
float HomogeneousAnisotropicMedium::density(const Vec3f &p) const
{
  return rho;
}

float HomogeneousAnisotropicMedium::sigma_t(const Vec3f &wi) const {
  return rho * projection(wi);
}

float HomogeneousAnisotropicMedium::sigma_s(const Vec3f &wi) const {
  return albedo * rho * projection(wi);
}

float HomogeneousAnisotropicMedium::projection(const Vec3f &wi) const {
  Vec3f wi_n = normalize(wi);
  Vec4f wi4 = {wi_n, 0};
  return sqrt(dot(wi4, S * wi4));
}

std::shared_ptr<const Medium> MediumInterface::getMedium(const Ray3f ray, const HitInfo &hit) const
{
    if (dot(hit.sn, ray.d) < 0)
      return inside;
    else
      return outside;
}

Color3f TrL(const Scene &scene, Sampler &sampler, const Ray3f &ray_)
{
  Ray3f ray = ray_.normalizeRay();
  float Tr = 1.0;

  // Advance through the ray in a straight line, one intersection at a time.
  while (true)
  {
    /** Loop invariants at this point:
     *  * ray.maxt = infinity
     *  * ray.medium is the medium that the ray is just about to travel on.
     *      (The neighborhood of ray.o at direction ray.d.
     *      NOT the intersection between the ray and the scene.)
     *  * Tr = transmittance between the initial ray origin (ray_.o) and the current ray origin (ray.o).
     */

    /* Part 1: Advance Tr to next intersection */
  
    HitInfo hit;
    bool hitSurface = scene.intersect(ray, hit);

    if (hitSurface) ray.maxt = length(hit.p - ray.o) + 2.0 * Epsilon;

    // Travel to intersection (or to infinity): multiply by medium transmittance
    if (ray.medium) Tr *= ray.medium->Tr(ray, sampler);

    // At this point, even though the Tr value has advanced to the next
    // intersection, ray.o has not. It will do so in Part 3.

    /* Part 2: Terminate ray on certain events */

    // hit an emitter
    if (hitSurface && hit.mat != nullptr)
      return hit.mat->isEmissive() ? Tr * hit.mat->emitted(ray, hit) : Color3f(0.0f);

    // if transmittance below threshold exit
    if (Tr < Epsilon) break;

    // escaped scene (assume no infinite lights)
    if (!hitSurface)
      return Tr * scene.background(ray);

    /* Part 3: Prepare for next iteration */

    // set next medium based on whether we are entering or exiting the surface
    if (hit.mi->IsMediumTransition())
      ray.medium = hit.mi->getMedium(ray, hit);

    // update ray origin to start at the intersection point
    ray.o = ray(hit.t + Epsilon); // Epsilon for extra safety against infinite loops of intersecting the same point.
                                  // Probably unnecessary, especially if ray.mint is already an Epsilon.
    // reset ray max to infinity
    ray.maxt = std::numeric_limits<float>::infinity();
  }

  return Color3f(0.0f);
}
