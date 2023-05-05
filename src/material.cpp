/*
    This file is part of Dirt, the Dartmouth introductory ray tracer.

    Copyright (c) 2017-2019 by Wojciech Jarosz

    Dirt is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Dirt is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#define _USE_MATH_DEFINES

#include <dirt/material.h>
#include <dirt/texture.h>
#include <dirt/parser.h>
#include <dirt/scene.h>
#include <dirt/surface.h>
#include <dirt/onb.h>
#include <cmath>

shared_ptr<const Material> Material::defaultMaterial()
{
	return nullptr;
}


inline bool refract(const Vec3f &v, const Vec3f &n, float iorIOverT, Vec3f &refracted, float &cosTheta2)
{
	Vec3f uv = normalize(v);
	float dt = dot(uv,n);
	float discrim = 1.0f - iorIOverT * iorIOverT * (1.0f - dt * dt);
	if (discrim > 0)
	{
		cosTheta2 = std::sqrt(discrim);
		refracted = iorIOverT * (uv - n * dt) - n * cosTheta2;
		return true;
	}
	else
	{
		return false;
	}
}

inline Vec3f reflect(const Vec3f &v, const Vec3f &n)
{
	return v - 2 * dot(v, n) * n;
}

inline float conductorFresnel(float cosine, float ior, float ext) 
{
	float tmp = (ior*ior + ext*ext);
	float R_s = (tmp - 2*ior*cosine + cosine*cosine) / (tmp + 2*ior*cosine + cosine*cosine);
	float R_p = (tmp*cosine*cosine - 2*ior*cosine + 1)/(tmp*cosine*cosine + 2*ior*cosine + 1);
	return (R_s + R_p)/2.0f;
}

Lambertian::Lambertian(const json & j)
{
	albedo = parseTexture(j.at("albedo"));
}

bool Lambertian::scatter(const Ray3f &ray, const HitInfo &hit, const Vec2f &sample, Color3f &attenuation, Ray3f &scattered) const
{
	// TODO: Implement Lambertian reflection
	//       You should assign the albedo to ``attenuation'', and
	//       you should assign the scattered ray to ``scattered''
	//       The origin of the scattered ray should be at the hit point,
	//       and the scattered direction is the shading normal plus a random
	//       point on a sphere (please look at the text book for this)
	//       You can get the hit point using hit.p, and the shading normal using hit.sn

	//       Hint: You can use the function randomInUnitSphere() to get a random
	//       point in a sphere. IMPORTANT: You want to add a random point *on*
	//       a sphere, not *in* the sphere (the text book gets this wrong)
	//       If you normalize the point, you can force it to be on the sphere always, so
	//       add normalize(randomInUnitSphere()) to your shading normal
	scattered = Ray3f(hit.p, hit.sn + randomOnUnitSphere(sample));
	attenuation = albedo->value(hit);
	return true;
}

bool Lambertian::sample(const Vec3f & dirIn, const HitInfo &hit, const Vec2f &sample, ScatterRecord &srec) const 
{
	ONBf onb;
	onb.build_from_w(hit.sn);
    srec.scattered = onb.toWorld(randomCosineHemisphere(sample));
	srec.isSpecular = false;
    srec.attenuation = albedo->value(hit);
	return true;
}

Color3f Lambertian::eval(const Vec3f & dirIn, const Vec3f &scattered, const HitInfo & hit) const 
{
	return albedo->value(hit) * max(0.0f, dot(scattered, hit.sn) / M_PI);
}

float Lambertian::pdf(const Vec3f & dirIn, const Vec3f &scattered, const HitInfo & hit) const 
{
	return max(0.0f, dot(scattered, hit.sn) / M_PI);
}

Metal::Metal(const json & j)
{
	albedo = parseTexture(j.at("albedo"));
	roughness = parseTexture(j.at("roughness"));
}

bool Metal::scatter(const Ray3f &ray, const HitInfo &hit, const Vec2f &sample, Color3f &attenuation, Ray3f &scattered) const
{
	// TODO: Implement metal reflection
	//       This function proceeds similar to the lambertian material, except that the
	//       scattered direction is different.
	//       Instead of adding a point on a sphere to the normal as before, you should add the point
	//       to the *reflected ray direction*.
	//       You can reflect a vector by the normal using reflect(vector, hit.sn); make sure the vector is normalized.
	//       Different to before you can't just use randomInUnitSphere directly; the sphere should be scaled by roughness.
	//       (see text book). In other words, if roughness is 0, the scattered direction should just be the reflected direction.
	//       
	//       This procedure could produce directions below the surface. Handle this by returning false if the scattered direction and the shading normal
	//       point in different directions (i.e. their dot product is negative)
	Vec3f reflected = reflect(normalize(ray.d), hit.sn);
	scattered = Ray3f(hit.p, reflected + luminance(roughness->value(hit)) * randomOnUnitSphere(sample));
	attenuation = albedo->value(hit);
	return (dot(scattered.d, hit.sn) > 0);
}

bool Metal::sample(const Vec3f & dirIn, const HitInfo &hit, const Vec2f &sample, ScatterRecord &srec) const 
{
	srec.scattered = reflect(normalize(dirIn), hit.sn);
	srec.attenuation = albedo->value(hit);
	srec.isSpecular = true;
	return true;
}

Dielectric::Dielectric(const json & j)
{
	ior = j.value("ior", ior);
}

bool Dielectric::scatter(const Ray3f &ray, const HitInfo &hit, const Vec2f &sample, Color3f &attenuation, Ray3f &scattered) const
{
	Vec3f normal;
	float eta1, eta2;

	// ensure ior and normal are correctly oriented for computing reflection and refraction
	if (dot(ray.d, hit.sn) > 0)
	{
		normal = -hit.sn;	
		eta1 = ior;
		eta2 = 1.0f;
	}
	else
	{
		normal = hit.sn;
		eta1 = 1.0f;
		eta2 = ior;
	}

	attenuation = Color3f(1);

	// compute reflected + refracted ray
	float cosTheta2, cosTheta1 = dot(ray.d, -normal) / length(ray.d);
	Vec3f refracted, reflected = reflect(ray.d, hit.sn);
	if (!refract(ray.d, normal, eta1 / eta2, refracted, cosTheta2))
	{
		// no refraction, only reflection
		scattered = Ray3f(hit.p, reflected);
		return true;
	}

	// compute fresnel solutions
	float rho_parallel = ((eta2 * cosTheta1) - (eta1 * cosTheta2)) / ((eta2 * cosTheta1) + (eta1 * cosTheta2));
	float rho_perp = ((eta1 * cosTheta1) - (eta2 * cosTheta2)) / ((eta1 * cosTheta1) + (eta2 * cosTheta2));
	float Freflected = (rho_parallel * rho_parallel + rho_perp * rho_perp) / 2.0f;

	// sample scattered or reflected ray
	scattered = sample.x < Freflected 
		? Ray3f(hit.p, reflected)
		: Ray3f(hit.p, refracted);

	return true;
}

bool Dielectric::sample(const Vec3f & dirIn, const HitInfo &hit, const Vec2f &sample, ScatterRecord &srec) const 
{
	Vec3f normal;
	float eta1, eta2;

	// ensure ior and normal are correctly oriented for computing reflection and refraction
	if (dot(dirIn, hit.sn) > 0)
	{
		normal = -hit.sn;	
		eta1 = ior;
		eta2 = 1.0f;
	}
	else
	{
		normal = hit.sn;
		eta1 = 1.0f;
		eta2 = ior;
	}
	
	srec.isSpecular = true;
	srec.attenuation = Color3f(1.f);

	// compute reflected + refracted ray
	float cosTheta2, cosTheta1 = dot(dirIn, -normal) / length(dirIn);
	Vec3f refracted, reflected = reflect(dirIn, hit.sn);
	if (!refract(dirIn, normal, eta1 / eta2, refracted, cosTheta2))
	{
    // no refraction, only reflection
    srec.scattered = reflected;
    return true;
	}

  // compute fresnel solutions
  float rho_parallel = ((eta2 * cosTheta1) - (eta1 * cosTheta2)) / ((eta2 * cosTheta1) + (eta1 * cosTheta2));
  float rho_perp = ((eta1 * cosTheta1) - (eta2 * cosTheta2)) / ((eta1 * cosTheta1) + (eta2 * cosTheta2));
  float Freflected = (rho_parallel * rho_parallel + rho_perp * rho_perp) / 2.0f;

  // sample scattered or reflected ray
  srec.scattered = sample.x < Freflected ? reflected : refracted;

  return true;
}

Color3f Dielectric::eval(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const 
{
  return Color3f(0.f);
}

float Dielectric::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const 
{
  return 1.f;
}

DiffuseLight::DiffuseLight(const json & j)
{
	emit = j.value("emit", emit);
}

Color3f DiffuseLight::emitted(const Ray3f &ray, const HitInfo &hit) const
{
	// only emit from the normal-facing side
	if (dot(ray.d, hit.sn) > 0)
		return Color3f(0,0,0);
	else
		return emit;
}

BlendMaterial::BlendMaterial(const json & j)
{
	a = parseMaterial(j.at("a"));
	b = parseMaterial(j.at("b"));
	amount = parseTexture(j.at("amount"));
}

bool BlendMaterial::scatter(const Ray3f &ray, const HitInfo &hit, const Vec2f &sample, Color3f &attenuation, Ray3f &scattered) const
{
	float t = luminance(amount->value(hit));
	if (randf() < t)
		return b->scatter(ray, hit, sample, attenuation, scattered);
	else
		return a->scatter(ray, hit, sample, attenuation, scattered);
}

bool BlendMaterial::sample(const Vec3f & dirIn, const HitInfo &hit, const Vec2f &sample, ScatterRecord &srec) const 
{
	float t = luminance(amount->value(hit));
	if (randf() < t)
		return b->sample(dirIn, hit, sample, srec);
	else
		return a->sample(dirIn, hit, sample, srec);
}

Color3f BlendMaterial::eval(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
    float t = luminance(amount->value(hit));
    return lerp(a->eval(dirIn, scattered, hit),
                b->eval(dirIn, scattered, hit), t);
}

float BlendMaterial::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const 
{
	float t = luminance(amount->value(hit));
    return lerp(a->pdf(dirIn, scattered, hit),
                b->pdf(dirIn, scattered, hit), t);
}


Phong::Phong(const json & j)
{
	albedo = parseTexture(j.at("albedo"));
	exponent = j.value("exponent", exponent);
}

bool Phong::scatter(const Ray3f &ray, const HitInfo &hit, const Vec2f &sample, Color3f &attenuation, Ray3f &scattered) const
{
	ScatterRecord srec;
	if (!this->sample(ray.d, hit, sample, srec))
		return false;
	attenuation = srec.attenuation;
	scattered = Ray3f(hit.p, srec.scattered);
	return true;
}

bool Phong::sample(const Vec3f & dirIn, const HitInfo &hit, const Vec2f &sample, ScatterRecord &srec) const 
{
	ONBf onb;
	onb.build_from_w(reflect(dirIn, hit.sn));
  	srec.scattered = onb.toWorld(randomCosinePowerHemisphere(exponent, sample));
	srec.isSpecular = false;
  	srec.attenuation = albedo->value(hit);
	return dot(hit.sn, srec.scattered) > 0;
}

Color3f Phong::eval(const Vec3f & dirIn, const Vec3f &scattered, const HitInfo & hit) const 
{
	return albedo->value(hit) * pdf(dirIn, scattered, hit);
}

float Phong::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const 
{
  Vec3f mirrorDir = normalize(reflect(dirIn, hit.sn));
	float cosine = max(dot(normalize(scattered), mirrorDir), 0.0f);
	return (exponent + 1.0f)/(2.0f*M_PI) * powf(cosine, exponent);
}

BlinnPhong::BlinnPhong(const json & j)
{
	albedo = parseTexture(j.at("albedo"));
	exponent = j.value("exponent", exponent);
}

bool BlinnPhong::scatter(const Ray3f &ray, const HitInfo &hit, const Vec2f &sample, Color3f &attenuation, Ray3f &scattered) const
{
	ScatterRecord srec;
	if (!this->sample(ray.d, hit, sample, srec))
		return false;
	attenuation = srec.attenuation;
	scattered = Ray3f(hit.p, srec.scattered);
	return true;
}

bool BlinnPhong::sample(const Vec3f & dirIn, const HitInfo &hit, const Vec2f &sample, ScatterRecord &srec) const 
{
	ONBf onb;
	onb.build_from_w(hit.sn);
	Vec3f normal = onb.toWorld(randomCosinePowerHemisphere(exponent, sample));
	srec.scattered = normalize(reflect(dirIn, normal));
	srec.isSpecular = false;
	srec.attenuation = albedo->value(hit);
	return dot(hit.sn, srec.scattered) > 0;
}

Color3f BlinnPhong::eval(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const 
{
	return albedo->value(hit) * pdf(dirIn, scattered, hit);
}

float BlinnPhong::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const 
{
	float iDotN = dot(-dirIn, hit.sn);
	if (iDotN <= 0.0f || dot(scattered, hit.sn) <= 0.0f)
		return 0.0f;
	Vec3f normal = normalize(-normalize(dirIn) + normalize(scattered));
	float cosine = max(dot(normal, hit.sn), 0.0f);
	float normalPdf = (exponent + 1.0f)/(2.0f*M_PI)*powf(cosine, exponent);
  return normalPdf/(4.0f*dot(-dirIn, normal));
}


Beckmann::Beckmann(const json & j)
{
	albedo = parseTexture(j.at("albedo"));
	ab = j.value("ab", ab);
	ior = j.value("ior", ior);
	ext = j.value("ext", ext);
}

bool Beckmann::scatter(const Ray3f &ray, const HitInfo &hit, const Vec2f &sample, Color3f &attenuation, Ray3f &scattered) const
{
	ScatterRecord srec;
	if (!this->sample(ray.d, hit, sample, srec))
		return false;
	attenuation = srec.attenuation;
	scattered = Ray3f(hit.p, srec.scattered);
	return true;
}

bool Beckmann::sample(const Vec3f & dirIn, const HitInfo &hit, const Vec2f &sample, ScatterRecord &srec) const 
{
	ONBf onb;
	onb.build_from_w(hit.sn);
	Vec3f hn = randomBeckmannNormal(ab, sample);
	Vec3f normal = onb.toWorld(hn);
	srec.scattered = normalize(reflect(dirIn, normal));
	srec.isSpecular = false;
	srec.attenuation = albedo->value(hit);
	return dot(hit.sn, srec.scattered) > 0;
}


Color3f Beckmann::eval(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const 
{
	ONBf onb;
	onb.build_from_w(hit.sn);
	Vec3f normal = normalize(-normalize(dirIn) + normalize(scattered));
	Vec3f v1 = -normalize(dirIn);
	Vec3f v2 = normalize(scattered);
	float G2 = G1(-v1, normal, hit) * G1(v2, normal, hit);
	float F = conductorFresnel(abs(dot(v2, normal)), ior, ext);

	float cosTheta = abs(dot(normal, hit.sn));
	float sinTheta = sqrt(1 - cosTheta*cosTheta);
	float tanTheta = sinTheta/cosTheta;
	float nPdf = (powf(M_E, -tanTheta*tanTheta / (ab*ab))) / (M_PI*ab*ab*cosTheta*cosTheta*cosTheta*cosTheta);

	return F * albedo->value(hit) * G2 * pdf(dirIn, scattered, hit);
}

float Beckmann::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const 
{
	float iDotN = dot(-dirIn, hit.sn);
	if (iDotN <= 0.0f || dot(scattered, hit.sn) <= 0.0f)
		return 0.0f;
	Vec3f normal = normalize(-normalize(dirIn) + normalize(scattered));
	
	float cosTheta = abs(dot(normal, hit.sn));
	float sinTheta = sqrt(1 - cosTheta*cosTheta);
	float tanTheta = sinTheta/cosTheta;
	float nPdf = (powf(M_E, -tanTheta*tanTheta / (ab*ab))) / (M_PI*ab*ab*cosTheta*cosTheta*cosTheta*cosTheta);
	return nPdf/(4.0f*abs(dot(normalize(dirIn), normalize(normal))));
	return abs(dot(normalize(hit.sn), normal)) * nPdf/(4.0f*abs(dot(normalize(dirIn), normalize(normal))));
}

float Beckmann::G1(Vec3f v1, Vec3f v2, const HitInfo & hit) const
{
	float cosTheta1 = abs(dot(v1, hit.sn));
	float sinTheta1 = sqrt(1 - cosTheta1*cosTheta1);
	float tanTheta1 = sinTheta1/cosTheta1;
	float alpha = 1/(ab*tanTheta1);

	if (alpha < 1.6f) {
		return (3.535f*alpha + 2.181f*alpha*alpha)/(1.0f+2.276f*alpha+2.577f*alpha*alpha);
	} else {
		return 1.0f;
	}
}


OrenNayar::OrenNayar(const json & j)
{
   albedo = parseTexture(j.at("albedo"));
   sigma = j.value("sigma", sigma);
}
 
bool OrenNayar::scatter(const Ray3f &ray, const HitInfo &hit, const Vec2f &sample, Color3f &attenuation, Ray3f &scattered) const
{
   scattered = Ray3f(hit.p, hit.sn + normalize(randomInUnitSphere()));
   attenuation = albedo->value(hit);
   return true;
}
 
bool OrenNayar::sample(const Vec3f & dirIn, const HitInfo &hit, const Vec2f &sample, ScatterRecord &srec) const
{
   ONBf onb;
   onb.build_from_w(hit.sn);
   srec.scattered = onb.toWorld(randomCosineHemisphere(sample));
   srec.isSpecular = false;
   srec.attenuation = albedo->value(hit);
   return true;
}
 
Color3f OrenNayar::eval(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
   float A = 1.0f - (sigma*sigma)/(2.0f*(sigma*sigma + 0.33f));
   float B = (0.45f*sigma*sigma)/(sigma*sigma + 0.09f);
  
   Vec3f inPar = dot(-dirIn, hit.sn)*hit.sn;
   Vec3f inOrth = normalize(-dirIn - inPar);
   Vec3f outPar = dot(scattered, hit.sn)*hit.sn;
   Vec3f outOrth = normalize(scattered - outPar);
 
   float cosTerm = max(-10.0f, dot(inOrth, outOrth));
 
   float sinTanTerm = 0.f;
   float cosI = abs(dot(-normalize(dirIn), normalize(hit.sn)));
   float cosO = abs(dot(normalize(scattered), normalize(hit.sn)));
   if (cosI > cosO) {
       sinTanTerm = sqrt(1 - cosO*cosO) * (sqrt(1 - cosI*cosI))/cosI;
   } else {
       sinTanTerm = sqrt(1 - cosI*cosI) * (sqrt(1 - cosO*cosO))/cosO;
   }

   return albedo->value(hit) * (A + B*cosTerm * sinTanTerm) * max(0.0f, dot(scattered, hit.sn) / M_PI);
}
 
float OrenNayar::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
   return max(0.0f, dot(scattered, hit.sn) / M_PI);
}

Layered::Layered(const json & j) {
	nb_layers = j.value("nb_layers", 1);

	// Add the external IOR (air)
	m_etas.push_back(Color3f(1.000277f));
	m_kappas.push_back(Color3f(0.f));

	for (auto it = j["layers"].begin(); it != j["layers"].end(); it++) {
		Color3f eta_k = it.value().value("eta", Color3f(1.0f));
		Color3f kappa_k = it.value().value("kappa", Color3f(0.0f));
		float alpha_k = it.value().value("alpha", 0.0f);

		// assume no participating media between interfaces yet
		float depth = 0.0f;
		float sigma_s = 0.0f;
		float sigma_a = 0.0f;
		float g = 0.9f;

		// store them in the vector
		m_etas.push_back(eta_k);
        m_kappas.push_back(kappa_k);
        m_alphas.push_back(alpha_k);

        // Update the media
        m_depths.push_back(depth);
        m_sigma_s.push_back(sigma_s);
        m_sigma_a.push_back(sigma_a);
        m_gs.push_back(g);
	}

	// load the FGD
	m_FGD = FGD("data/FGD.bin");
}

/* Roughness to linear space conversions */
inline float roughnessToVariance(float a) {
   return a / (1.0f-a);
}
inline float varianceToRoughness(float v) {
   return v / (1.0f+v);
}

inline float avg(const Color3f& color) {
	return (1.f / 3.f) * (color.r + color.g + color.b);
};

// Source: https://stackoverflow.com/questions/27229371/inverse-error-function-in-c
float erfinv(float x){
   float tt1, tt2, lnx, sgn;
   sgn = (x < 0) ? -1.0f : 1.0f;

   x = (1 - x)*(1 + x);        // x = 1 - x*x;
   lnx = logf(x);

   tt1 = 2/(M_PI*0.147) + 0.5f * lnx;
   tt2 = 1/(0.147) * lnx;

   return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}

/* Beckmann Distribution for Microfacet Model (Mitsuba Version but in world coordinate) */

// sn is the surface normal, n is the sampled normal. Assume they are both normalized.
inline float BechmannEval(const Vec3f& sn, const Vec3f& n, const float alpha) {
	if (dot(sn, n) <= 0)
            return 0.0f;

	float cosTheta2 = dot(sn, n) * dot(sn, n);
	float beckmannExponent = ((1 - cosTheta2) / (alpha * alpha)) / cosTheta2;

	float result;
	/* Beckmann distribution function for Gaussian random surfaces - [Walter 2005] evaluation */
	result = exp(-beckmannExponent) /
		(M_PI * alpha * alpha * cosTheta2 * cosTheta2);

	/* Prevent potential numerical issues in other stages of the model */
	if (result * dot(sn, n) < 1e-20f)
		result = 0;

	return result;
}

/* Smith's shadowing-masking function G1 for each of the supported microfacet distributions */
// n is a sampled microfacet normal, v is an arbitrary vector
inline float smithG1(const Vec3f& sn, const Vec3f& n, const float alpha, const Vec3f& v) {
	/* Ensure consistent orientation (can't see the back
		of the microfacet from the front and vice versa) */
	if (dot(v, n) * dot(v, sn) <= 0)
		return 0.0f;

	/* Perpendicular incidence -- no shadowing/masking */
	float sinTheta = length(cross(v, sn));
	float cosTheta = dot(v, sn);
	float tanTheta = abs(sinTheta/cosTheta);
	if (tanTheta == 0.0f)
		return 1.0f;

	float a = 1.0f / (alpha * tanTheta);
	if (a >= 1.6f)
		return 1.0f;

	/* Use a fast and accurate (<0.35% rel. error) rational
		approximation to the shadowing-masking function */
	float aSqr = a*a;
	return (3.535f * a + 2.181f * aSqr)
			/ (1.0f + 2.276f * a + 2.577f * aSqr);
}

/* Draw a sample from the distribution of visible normals */
inline Vec3f sampleVisible(const Vec3f& sn, const Vec3f &_wi, const float alpha, const Vec2f &sample) {
	/* Step 1: stretch wi */
	ONBf onb;
	onb.build_from_w(sn);
    Vec3f wi = onb.toLocal(_wi);

	wi = normalize(Vec3f(
		alpha * wi.x,
		alpha * wi.y,
		wi.z
	));

	/* Get polar coordinates */
	float theta = 0, phi = 0;
	if (wi.z < (float) 0.99999) {
		theta = std::acos(wi.z);
		phi = std::atan2(wi.y, wi.x);
	}
	float sinPhi = std::sin(phi);
	float cosPhi = std::cos(phi);

	/* Step 2: simulate P22_{wi}(slope.x, slope.y, 1, 1) */
	Vec2f slope;
	{
        const float SQRT_PI_INV = 1 / std::sqrt(M_PI);

        /* Special case (normal incidence) */
		if (theta < 1e-4f) {
			float sine = sin(2 * M_PI * sample.y);
			float cosine = cos(2 * M_PI * sample.y);
			float r = std::sqrt(-log(1.0f-sample.x));
			slope = Vec2f(r * cosine, r * sine);
		} else {
			/* The original inversion routine from the paper contained
				discontinuities, which causes issues for QMC integration
				and techniques like Kelemen-style MLT. The following code
				performs a numerical inversion with better behavior */
			float tanThetaI = std::tan(theta);
			float cotThetaI = 1 / tanThetaI;

			/* Search interval -- everything is parameterized
				in the erf() domain */
			float a = -1, c = erf(cotThetaI);
			float sample_x = std::max(sample.x, (float) 1e-6f);

			/* Start with a good initial guess */
			//Float b = (1-sample_x) * a + sample_x * c;

			/* We can do better (inverse of an approximation computed in Mathematica) */
			float fit = 1 + theta*(-0.876f + theta * (0.4265f - 0.0594f*theta));
			float b = c - (1+c) * std::pow(1-sample_x, fit);

			/* Normalization factor for the CDF */
			float normalization = 1 / (1 + c + SQRT_PI_INV*
				tanThetaI*std::exp(-cotThetaI*cotThetaI));

			int it = 0;
			while (++it < 10) {
				/* Bisection criterion -- the oddly-looking
					boolean expression are intentional to check
					for NaNs at little additional cost */
				if (!(b >= a && b <= c))
					b = 0.5f * (a + c);

				/* Evaluate the CDF and its derivative
					(i.e. the density function) */
				float invErf = erfinv(b);
				float value = normalization*(1 + b + SQRT_PI_INV*
					tanThetaI*std::exp(-invErf*invErf)) - sample_x;
				float derivative = normalization * (1
					- invErf*tanThetaI);

				if (std::abs(value) < 1e-5f)
					break;

				/* Update bisection intervals */
				if (value > 0)
					c = b;
				else
					a = b;

				b -= value / derivative;
			}

			/* Now convert back into a slope value */
			slope.x = erfinv(b);

			/* Simulate Y component */
			slope.y = erfinv(2.0f*std::max(sample.y, (float) 1e-6f) - 1.0f);
		}
    }

	/* Step 3: rotate */
	slope = Vec2f(
		cosPhi * slope.x - sinPhi * slope.y,
		sinPhi * slope.x + cosPhi * slope.y);

	/* Step 4: unstretch */
	slope.x *= alpha;
	slope.y *= alpha;

	/* Step 5: compute normal */
	float normalization = (float) 1 / std::sqrt(slope.x*slope.x
			+ slope.y*slope.y + (float) 1.0);

	Vec3f result = Vec3f(
		-slope.x * normalization,
		-slope.y * normalization,
		normalization
	);

	// Convert back go global coordinate
	return normalize(onb.toWorld(result));
}

bool Layered::scatter(const Ray3f &ray, const HitInfo &hit, const Vec2f &sample, Color3f &attenuation, Ray3f &scattered) const
{
	return false;
}

bool Layered::sample(const Vec3f & dirIn, const HitInfo &hit, const Vec2f &sample, ScatterRecord &srec) const
{
	/* Constants */
	Vec3f wi = normalize(-dirIn);
	Vec3f n = normalize(hit.sn);

	if (dot(wi, n) < 0) {
		return false;
	}

	/* Evaluate the adding method to get coeffs and variances */
	Color3f* coeffs = new Color3f[nb_layers];
	float* alphas = new float[nb_layers];
	int nb_valid = 0;
	float normalProj = dot(wi, n);
	computeAddingDoubling(normalProj, coeffs, alphas, nb_valid);

	/* Convert Spectral coefficients to floats to select BRDF lobe to sample */
	float* weights = new float[nb_valid];
	float cum_w = 0.0;
	for(int i=0; i<nb_valid; ++i) {
		weights[i] = avg(coeffs[i]);
		cum_w += weights[i];
	}

	/* Select a random BRDF lobe */
	float sel_w = sample.x * cum_w - weights[0];
	int   sel_i = 0;
	for(sel_i=0; sel_w>0.0 && sel_i<nb_valid; sel_i++) {
		sel_w -= weights[sel_i+1];
	}

	/* Sample a microfacet normal */
	Vec3f m = sampleVisible(n, wi, alphas[sel_i], Vec2f(randf(), randf()));

	delete[] coeffs;
	delete[] alphas;
	delete[] weights;

	/* Perfect specular reflection based on the microfacet normal */
	srec.scattered = reflect(-wi, m);
	if(dot(srec.scattered, n) <= 0.0f) {
		return false;
	}

	srec.isSpecular = false;
	// not used by Monte Carlo path tracer，just set to 1
	srec.attenuation = Vec3f(1.f, 1.f, 1.f);

	return true;
}

Color3f Layered::eval(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	/* Constants */
	Vec3f wi = normalize(-dirIn);
	Vec3f wo = normalize(scattered);
	Vec3f n = normalize(hit.sn);
	Vec3f H = normalize(wo+wi);

	/* Result */
	Color3f f(0.0f);

	/* Evaluate the adding method to get coeffs and variances */
	Color3f* coeffs = new Color3f[nb_layers];
	float* alphas = new float[nb_layers];
	int nb_valid = 0;
	float normalProj = dot(wi, n);
	computeAddingDoubling(normalProj, coeffs, alphas, nb_valid);

	/* Sum the contribution of all the interfaces */
	for(int index=0; index<nb_valid; ++index) {
		// Skip zero contributions
		if(length(coeffs[index]) == 0) {
			continue;
		}

		// Fetch current roughness
		const float a = alphas[index];

		// Evaluate microfacet model
		const float D = BechmannEval(n, H, a);
		const float G = smithG1(n, H, a, wi) * smithG1(n, H, a, wo);

		// Add to the contribution
		f += D*G * coeffs[index] / (4.0f * dot(wi, n));
	}

	delete[] coeffs;
	delete[] alphas;

	return f;
}

float Layered::pdf(const Vec3f & dirIn, const Vec3f & scattered, const HitInfo & hit) const
{
	Vec3f wi = normalize(-dirIn);
	Vec3f wo = normalize(scattered);
	Vec3f n = normalize(hit.sn);

	if (dot(wi, n) <= 0 || dot(wo, n) <= 0)
		return 0.0f;

	/* Calculate the reflection half-vector */
	Vec3f H = normalize(wo + wi);

	/* Evaluate the adding method to get coeffs and variances */
	Color3f* coeffs = new Color3f[nb_layers];
	float* alphas = new float[nb_layers];
	int nb_valid = 0;

	float normalProj = dot(wi, n);
	computeAddingDoubling(normalProj, coeffs, alphas, nb_valid);

	/* Convert Spectral coefficients to float for pdf weighting */
	float pdf = 0.0;
	float cum_w = 0.0;
	for(int i=0; i<nb_valid; ++i) {
		// Skip zero contributions
		if(length(coeffs[i]) == 0) {
			continue;
		}

		/* Evlaluate weight */
		auto weight = avg(coeffs[i]);
		cum_w += weight;

		/* Use Beckmann as NDF */

		/* Evaluate the pdf */
		float DG = BechmannEval(n, H, alphas[i]) * smithG1(n, H, alphas[i], wo) / (4.0f * dot(wi, n));
		pdf += weight * DG;
	}

	delete[] coeffs;
	delete[] alphas;
	
	if(cum_w > 0.0f) {
		return pdf / cum_w;
	} else {
		return 0.0f;
	}
}

void Layered::evalFresnel(float ct, float alpha, const Color3f& eta, const Color3f& kappa, Color3f& Rij, Color3f& Tij) const
{
	Rij = m_FGD(ct, alpha, eta, kappa);
	Tij = (length(kappa) == 0) ? Color3f(1.0) - Rij : Color3f(0.0);
}

void Layered::computeAddingDoubling(float _cti, Color3f* coeffs, float* alphas, int& nb_valid) const
{
	// Variables
	float cti  = _cti;
	Color3f R0i(0.0f), Ri0(0.0f), T0i(1.0f), Ti0(1.0f);
	float s_r0i = 0.0f, s_ri0=0.0f, s_t0i=0.0f, s_ti0=0.0f;
	float j0i=1.0f, ji0=1.0f;

	// Iterate over the layers
	for(int i=0; i<nb_layers; ++i) {

		/* Extract layer data */
		Color3f eta_1   = m_etas[i];
		Color3f eta_2   = m_etas[i+1];
		Color3f kappa_2 = m_kappas[i+1];
		Color3f eta     = eta_2 / eta_1;
		Color3f kappa   = kappa_2 / eta_1;
		float alpha      = m_alphas[i];
		float n12        = avg(eta);
		float depth      = m_depths[i];

		Color3f R12, T12, R21, T21;
		float s_r12=0.0f, s_r21=0.0f, s_t12=0.0f, s_t21=0.0f, j12=1.0f, j21=1.0f, ctt;
		if(depth > 0.0f) {
			/* Mean doesn't change with volumes */
			ctt = cti;

			/* Evaluate transmittance */
			const float sigma_t = m_sigma_a[i] + m_sigma_s[i];
			T12 = (Color3f(1.0f) + m_sigma_s[i]*depth/ctt) * exp(- (depth/ctt) * sigma_t);
			T21 = T12;
			R12 = Color3f(0.0f);
			R21 = Color3f(0.0f);

			/* Fetch precomputed variance for HG phase function */
			s_t12 = alpha;
			s_t21 = alpha;

		} else {
			/* Evaluate off-specular transmission */
			float sti = sqrt(1.0f - cti*cti);
			float stt = sti / n12;
			if(stt <= 1.0f) {
				ctt = sqrt(1.0f - stt*stt);
			} else {
				ctt = -1.0f;
			}

			/* Ray is not block by conducting interface or total reflection */
			const bool has_transmissive = ctt > 0.0f && length(kappa) == 0;

			/* Evaluate interface variance term */
			s_r12 = roughnessToVariance(alpha);
			s_r21 = s_r12;

			/* For dielectric interfaces, evaluate the transmissive roughnesses */
			if(has_transmissive) {
				const float _ctt = 1.0f; // The scaling factor overblurs the BSDF at grazing
				const float _cti = 1.0f; // angles (we cannot account for the deformation of
											// the lobe for those configurations.

				s_t12 = roughnessToVariance(alpha * 0.5f * fabs(_ctt*n12 - _cti)/(_ctt*n12));
				s_t21 = roughnessToVariance(alpha * 0.5f * fabs(_cti/n12 - _ctt)/(_cti/n12));
					j12 = (ctt/cti) * n12; // Scale due to the interface
					j21 = (cti/ctt) / n12;
			}

			/* Evaluate FGD using a modified roughness accounting for top layers */
			auto temp_alpha = varianceToRoughness(s_t0i + s_r12);

			/* Evaluate r12, r21, t12, t21 */
			evalFresnel(cti, temp_alpha, eta, kappa, R12, T12);
			if(has_transmissive) {
				R21 = R12;
				T21 = T12 ; // We don't need the IOR scaling since we are
				T12 = T12 ; // computing reflectance only here.
			} else {
				R21 = Color3f(0.0f);
				T21 = Color3f(0.0f);
				T12 = Color3f(0.0f);
			}
		}

		/* Multiple scattering forms */
		const Color3f denom = (Color3f(1.0f) - Ri0*R12);
		const Color3f m_R0i = (avg(denom) <= 0.0f)? Color3f(0.0f) : (T0i*R12*Ti0) / denom;
		const Color3f m_Ri0 = (avg(denom) <= 0.0f)? Color3f(0.0f) : (T21*Ri0*T12) / denom;
		const Color3f m_Rr  = (avg(denom) <= 0.0f)? Color3f(0.0f) : (Ri0*R12) / denom;
		
		/* Evaluate the adding operator on the energy */
		const Color3f e_R0i = R0i + m_R0i;
		const Color3f e_T0i = (T0i*T12) / denom;
		const Color3f e_Ri0 = R21 + m_Ri0;
		const Color3f e_Ti0 = (T21*Ti0) / denom;

		/* Scalar forms for the spectral quantities */
		const float r0i   = avg(R0i);
		const float e_r0i = avg(e_R0i);
		const float e_ri0 = avg(e_Ri0);
		const float m_r0i = avg(m_R0i);
		const float m_ri0 = avg(m_Ri0);
		const float m_rr  = avg(m_Rr);
		const float r21   = avg(R21);

		/* Evaluate the adding operator on the normalized variance */
		float _s_r0i = (r0i*s_r0i + m_r0i*(s_ti0 + j0i*(s_t0i + s_r12 + m_rr*(s_r12+s_ri0)))) ;// e_r0i;
		float _s_t0i = j12*s_t0i + s_t12 + j12*(s_r12 + s_ri0)*m_rr;
		float _s_ri0 = (r21*s_r21 + m_ri0*(s_t12 + j12*(s_t21 + s_ri0 + m_rr*(s_r12+s_ri0)))) ;// e_ri0;
		float _s_ti0 = ji0*s_t21 + s_ti0 + ji0*(s_r12 + s_ri0)*m_rr;
		_s_r0i = (e_r0i > 0.0f) ? _s_r0i/e_r0i : 0.0f;
		_s_ri0 = (e_ri0 > 0.0f) ? _s_ri0/e_ri0 : 0.0f;

		/* Store the coefficient and variance */
		if(m_r0i > 0.0f) {
			coeffs[i] = m_R0i;
			alphas[i] = varianceToRoughness(s_ti0 + j0i*(s_t0i + s_r12 + m_rr*(s_r12+s_ri0)));
		} else {
			coeffs[i] = Color3f(0.0f);
			alphas[i] = 0.0f;
		}

		/* Update energy */
		R0i = e_R0i;
		T0i = e_T0i;
		Ri0 = e_Ri0;
		Ti0 = e_Ti0;

		/* Update mean */
		cti = ctt;

		/* Update variance */
		s_r0i = _s_r0i;
		s_t0i = _s_t0i;
		s_ri0 = _s_ri0;
		s_ti0 = _s_ti0;

		/* Update jacobian */
		j0i *= j12;
		ji0 *= j21;

		/* Escape if a conductor is present */
		if(avg(kappa) > 0.0f) {
			nb_valid = i+1;
			return;
		}
	}

	nb_valid = nb_layers;
}