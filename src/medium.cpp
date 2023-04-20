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
