
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// media/homogeneous.cpp*
#include "media/homogeneous.h"
#include "sampler.h"
#include "interaction.h"
#include "paramset.h"
#include "stats.h"
#include "texture.h"

namespace pbrt {

STAT_COUNTER("Integrator/homogeneous lookups", densityLookups);

// HomogeneousMedium Method Definitions
HomogeneousMedium::HomogeneousMedium(const Spectrum &sigma_a, 
                                     const Spectrum &sigma_s, 
                                     const Spectrum &sigma_n,
                                     Float g, bool enableDualDensity,
                                     bool enableMarble,
                                     Point3f noiseScale, 
                                     Float noiseScale2, 
                                     Float a,
                                     Float w,
                                     int num_octaves,
                                     const Spectrum &sigma_a_2, 
                                     const Spectrum &sigma_s_2,
                                     const Spectrum &sigma_n_2,
                                     const MarbleScanlineOpts marbleScanlineOpts)
	: sigma_a(sigma_a)
	, sigma_s(sigma_s)
	, sigma_t(sigma_a + sigma_s + sigma_n)
	, g(g)
  , m_enableDualDensity(enableDualDensity)
  , m_enabeMarble(enableMarble)
  , m_noiseScale(noiseScale)
  , m_noiseScale2(noiseScale2)
  , m_a(a)
  , m_w(w)
  , m_num_octaves(num_octaves)
  , sigma_a_2(sigma_a_2)
  , sigma_s_2(sigma_s_2)
  , sigma_t_2(sigma_a_2 + sigma_s_2 + sigma_n_2)
  , marbleScanline(MarbleScanline(marbleScanlineOpts))
  {
    for (int i = 0; i < Spectrum::nSamples; i++) {
      max_sigma_t[i] = std::max(sigma_t[i], sigma_t_2[i]);
    }
  }

Spectrum HomogeneousMedium::Tr(const RayDifferential &ray, Sampler &sampler, uint32_t flags, TransportMode mode) const {
    ProfilePhase _(Prof::MediumTr);
    return Exp(-sigma_t * std::min(ray.tMax * ray.d.Length(), MaxFloat));
}

Spectrum HomogeneousMedium::Tr(const Ray& ray, Sampler& sampler, uint32_t flags, TransportMode mode) const
{
    return Tr(RayDifferential(ray), sampler, flags, mode);
}

Spectrum HomogeneousMedium::Sample(const RayDifferential &ray, Sampler &sampler, MemoryArena &arena, MediumInteraction *mi, uint32_t flags, TransportMode mode) const {
    ProfilePhase _(Prof::MediumSample);
    int channel = std::min((int)(sampler.Get1D() * Spectrum::nSamples),
            Spectrum::nSamples - 1);
    Float dist = -std::log(1 - sampler.Get1D()) / sigma_t[channel];
    Float t = std::min(dist / ray.d.Length(), ray.tMax);
    bool sampledMedium = t < ray.tMax;
    if (sampledMedium)
            *mi = MediumInteraction(ray(t), -ray.d, ray.time, this,
                    ARENA_ALLOC(arena, HenyeyGreenstein)(g));
      
    // Compute the transmittance and sampling density
    Spectrum Tr = Exp(-sigma_t * std::min(t, MaxFloat) * ray.d.Length());

    // Return weighting factor for scattering from homogeneous medium
    Spectrum density = sampledMedium ? (sigma_t * Tr) : Tr;
    Float pdf = 0;
    for (int i = 0; i < Spectrum::nSamples; ++i) pdf += density[i];
    pdf *= 1 / (Float)Spectrum::nSamples;
    if (pdf == 0) {
            CHECK(Tr.IsBlack());
            pdf = 1;
    }
    return sampledMedium ? (Tr * sigma_s / pdf) : (Tr / pdf);
}

Spectrum HomogeneousMedium::SampleChannel(const Ray &ray, Sampler &sampler, MemoryArena &arena, MediumInteraction *mi, int channel, uint32_t flags, TransportMode mode) const 
{
    ProfilePhase _(Prof::MediumSample);
		Ray normRay = Ray(ray.o, Normalize(ray.d), ray.d.Length() * ray.tMax);
    Float dist = -std::log(1 - sampler.Get1D()) / GetMajorant()[channel];
    Float t = std::min(dist, normRay.tMax);
    bool sampledMedium = t < normRay.tMax;
    if (sampledMedium)
            *mi = MediumInteraction(normRay(t), -normRay.d, normRay.time, this,
                    ARENA_ALLOC(arena, HenyeyGreenstein)(g));

    return Exp(-GetMajorant() * t);
}

Spectrum HomogeneousMedium::Sample(const Ray& ray, Sampler& sampler, MemoryArena& arena, MediumInteraction* mi, uint32_t flags, TransportMode mode) const
{
	return Sample(RayDifferential(ray), sampler, arena, mi, flags, mode);
}

Spectrum HomogeneousMedium::GetAbsorption(const Point3f &p) const {
  densityLookups++;
  if (m_enableDualDensity) {
    Float noise = Noise(p * m_noiseScale);
    Float m = (noise + 1.0f) / 2.0f;
    return m * sigma_a + (1.0f - m) * sigma_a_2;
  } else if (m_enabeMarble) {
    Float marble = GetMarble(p);
    return marble * sigma_a + (1.0f - marble) * sigma_a_2;
  }
  return sigma_a;
};

Spectrum HomogeneousMedium::GetScattering(const Point3f &p) const {
  densityLookups++;
  if (m_enableDualDensity) {
    Float noise = Noise(p * m_noiseScale);
    Float m = (noise + 1.0f) / 2.0f;
    return m * sigma_s + (1.0f - m) * sigma_s_2;
  } else if (m_enabeMarble) {
    Float marble = GetMarble(p);
    return marble * sigma_s + (1.0f - marble) * sigma_s_2;
  }
  return sigma_s;
};

void HomogeneousMedium::GetCoefficients(const Point3f &p, Spectrum &s, Spectrum &a, Spectrum &n) const {
  densityLookups++;
  if (m_enableDualDensity) {
    Float noise = Noise(p * m_noiseScale);
    Float m = (noise + 1.0f) / 2.0f;
    s = m * sigma_s + (1.0f - m) * sigma_s_2;
    a = m * sigma_a + (1.0f - m) * sigma_a_2;
  } else if (m_enabeMarble) {
    Float marble = GetMarble(p);
    s = marble * sigma_s + (1.0f - marble) * sigma_s_2;
    a = marble * sigma_a + (1.0f - marble) * sigma_a_2;
  } else {
    s = Spectrum(sigma_s);
    a = Spectrum(sigma_a);
  }
  n = GetMajorant() - (s + a);
};

Spectrum HomogeneousMedium::GetMajorant() const {
  return !m_enableDualDensity || !m_enabeMarble ? sigma_t : max_sigma_t;
};

Float HomogeneousMedium::GetMarble(const Point3f &p) const {
  Float noise = 0.0f;
  Float w = 1.0f, omega = 1.0f;
  for  (int i = 0; i < m_num_octaves; i++) {
    noise +=  m_a * w * Noise(p* m_noiseScale / omega);
    w *= m_w;
    omega *= 1.99f;
  }
  return marbleScanline.get((p.x + noise) * m_noiseScale2);
};

}  // namespace pbrt
