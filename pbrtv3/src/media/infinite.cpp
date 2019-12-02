
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


// media/infinite.cpp*
#include "media/infinite.h"
#include "sampler.h"
#include "interaction.h"
#include "paramset.h"
#include "stats.h"
#include "texture.h"

namespace pbrt {

// HomogeneousMedium Method Definitions
InfiniteMedium::InfiniteMedium(const Spectrum &sigma_a, 
															 const Spectrum &sigma_s, 
															 const Spectrum &sigma_n,
															 Float g)
	: sigma_a(sigma_a), 
		sigma_s(sigma_s), 
		sigma_t(sigma_a + sigma_s + sigma_n), 
		g(g) { }

Spectrum InfiniteMedium::Tr(const RayDifferential &ray, Sampler &sampler, uint32_t flags, TransportMode mode) const {
    ProfilePhase _(Prof::MediumTr);
    return Exp(-sigma_t * std::min(ray.tMax * ray.d.Length(), MaxFloat));
}

Spectrum InfiniteMedium::Tr(const Ray& ray, Sampler& sampler, uint32_t flags, TransportMode mode) const
{
    return Tr(RayDifferential(ray), sampler, flags, mode);
}

Spectrum InfiniteMedium::Sample(const RayDifferential &ray, Sampler &sampler, MemoryArena &arena, MediumInteraction *mi, uint32_t flags, TransportMode mode) const {
    ProfilePhase _(Prof::MediumSample);
    int channel = std::min((int)(sampler.Get1D() * Spectrum::nSamples),
            Spectrum::nSamples - 1);
    Float dist = -std::log(1 - sampler.Get1D()) / sigma_t[channel];
    Float t = std::min(dist / ray.d.Length(), ray.tMax);
		*mi = MediumInteraction(ray(t), -ray.d, ray.time, this,
														ARENA_ALLOC(arena, HenyeyGreenstein)(g));
      
    // Compute the transmittance and sampling density
    Spectrum Tr = Exp(-sigma_t * std::min(t, MaxFloat) * ray.d.Length());

    // Return weighting factor for scattering from homogeneous medium
    Spectrum density = (sigma_t * Tr);
    Float pdf = 0;
    for (int i = 0; i < Spectrum::nSamples; ++i) pdf += density[i];
    pdf *= 1 / (Float)Spectrum::nSamples;
    if (pdf == 0) {
			CHECK(Tr.IsBlack());
			pdf = 1;
    }
    return (Tr * sigma_s / pdf);
}

Spectrum InfiniteMedium::SampleChannel(const Ray &ray, Sampler &sampler, MemoryArena &arena, MediumInteraction *mi, int channel, uint32_t flags, TransportMode mode) const 
{
    ProfilePhase _(Prof::MediumSample);
    Float dist = -std::log(1 - sampler.Get1D()) / GetMajorant()[channel];
    Float t = dist / ray.d.Length();
		*mi = MediumInteraction(ray(t), -ray.d, ray.time, this,
														ARENA_ALLOC(arena, HenyeyGreenstein)(g));
    return Exp(-GetMajorant() * t * ray.d.Length());
}

Spectrum InfiniteMedium::Sample(const Ray& ray, Sampler& sampler, MemoryArena& arena, MediumInteraction* mi, uint32_t flags, TransportMode mode) const
{
	return Sample(RayDifferential(ray), sampler, arena, mi, flags, mode);
}

Spectrum InfiniteMedium::GetAbsorption(const Point3f &p) const {
  return sigma_a;
};

Spectrum InfiniteMedium::GetScattering(const Point3f &p) const {
  return sigma_s;
};

Spectrum InfiniteMedium::GetMajorant() const {
  return sigma_t;
};

}  // namespace pbrt
