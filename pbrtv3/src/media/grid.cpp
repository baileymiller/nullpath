
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


// media/grid.cpp*
#include "media/grid.h"
#include "paramset.h"
#include "sampler.h"
#include "stats.h"
#include "interaction.h"

namespace pbrt {

      STAT_RATIO("Media/Grid steps per Tr() call", nTrSteps, nTrCalls);
      STAT_COUNTER("Integrator/grid lookups", densityLookups);

      Spectrum GridDensityMedium::SampleChannel(const Ray &rWorld, 
        Sampler &sampler, MemoryArena &arena, MediumInteraction *mi, 
        int channel, uint32_t flags, TransportMode mode) const {
        ProfilePhase _(Prof::MediumSample);
        Ray normRay = Ray(rWorld.o, Normalize(rWorld.d), rWorld.d.Length() * rWorld.tMax);
        Ray ray = WorldToMedium(normRay);
	      const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
        Float tMin, tMax;
	      if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(1.f);
        Float t = tMin - std::log(1 - sampler.Get1D()) /GetMajorant()[channel];
        bool sampledMedium = t < tMax;
        t = std::min(t, tMax);
        if (sampledMedium)
          *mi = MediumInteraction(normRay(t), -normRay.d, normRay.time, 
                              this, ARENA_ALLOC(arena, HenyeyGreenstein)(g));
        return Exp(-GetMajorant() * (t-tMin));
      }
       
      Spectrum GridDensityMedium::GetAbsorption(const Point3f &p) const {
        Spectrum sig_a;
        if (!m_enableDualDensity) sig_a = sigma_a;
        else sig_a = Interpolate(sigma_a, sigma_a_2, p);
        return Density(WorldToMedium(p)) * sig_a;
      }

      Spectrum GridDensityMedium::GetScattering(const Point3f &p) const {
        Spectrum sig_s;
        if (!m_enableDualDensity) sig_s = sigma_s;
        else sig_s = Interpolate(sigma_s, sigma_s_2, p);
        return Density(WorldToMedium(p)) * sig_s;
      }

      void GridDensityMedium::GetCoefficients(const Point3f &p, Spectrum &s, Spectrum &a, Spectrum &n) const {
        densityLookups++;
        Float density = Density(WorldToMedium(p));
        if (!m_enableDualDensity) {
          s = sigma_s * density;
          a = sigma_a * density;
        } else {
          s = Interpolate(sigma_s, sigma_s_2, p) * density;
          a = Interpolate(sigma_a, sigma_a_2, p) * density;
        }
        n = GetMajorant() - (s + a);
      }

      Spectrum GridDensityMedium::GetMajorant() const {
        return sigma_t / invMaxDensity;
      }

      Spectrum GridDensityMedium::Interpolate(Spectrum sigma_x, Spectrum sigma_y, Point3f p) const {
        Point3f msP = WorldToMedium(p);
        if (m_interpType == InterpType::HEIGHT_Y) {
          return msP.y * sigma_x + (1.0f - msP.y) * sigma_y;
        } else {
          return sigma_a;
        }
      }

	// GridDensityMedium Method Definitions
	Float GridDensityMedium::Density(const Point3f &p) const {
		// Compute voxel coordinates and offsets for _p_
		Point3f pSamples(p.x * nx - .5f, p.y * ny - .5f, p.z * nz - .5f);
		Point3i pi = (Point3i)Floor(pSamples);
		Vector3f d = pSamples - (Point3f)pi;

		// Trilinearly interpolate density values to compute local density
		Float d00 = Lerp(d.x, D(pi), D(pi + Vector3i(1, 0, 0)));
		Float d10 = Lerp(d.x, D(pi + Vector3i(0, 1, 0)), D(pi + Vector3i(1, 1, 0)));
		Float d01 = Lerp(d.x, D(pi + Vector3i(0, 0, 1)), D(pi + Vector3i(1, 0, 1)));
		Float d11 = Lerp(d.x, D(pi + Vector3i(0, 1, 1)), D(pi + Vector3i(1, 1, 1)));
		Float d0 = Lerp(d.y, d00, d10);
		Float d1 = Lerp(d.y, d01, d11);
		return Lerp(d.z, d0, d1);
	}

	Spectrum GridDensityMedium::Sample(const Ray &rWorld, Sampler &sampler,
		MemoryArena &arena,
		MediumInteraction *mi, uint32_t flags, TransportMode mode) const {
		ProfilePhase _(Prof::MediumSample);
		Ray ray = WorldToMedium(
			Ray(rWorld.o, Normalize(rWorld.d), rWorld.tMax * rWorld.d.Length()));
		// Compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
		const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
		Float tMin, tMax;
		if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(1.f);

		// Run delta-tracking iterations to sample a medium interaction
		Float t = tMin;
		while (true) {
			t -= std::log(1 - sampler.Get1D()) * invMaxDensity / sigma_t[0];
			if (t >= tMax) break;
			if (Density(ray(t)) * invMaxDensity > sampler.Get1D()) {
				// Populate _mi_ with medium interaction information and return
				PhaseFunction *phase = ARENA_ALLOC(arena, HenyeyGreenstein)(g);
				*mi = MediumInteraction(rWorld(t), -rWorld.d, rWorld.time, this,
					phase);
				return sigma_s / sigma_t[0];
			}
		}
		return Spectrum(1.f);
	}

	Spectrum GridDensityMedium::Tr(const Ray &rWorld, Sampler &sampler, uint32_t flags, TransportMode mode) const {
		ProfilePhase _(Prof::MediumTr);
		++nTrCalls;

		Ray ray = WorldToMedium(
			Ray(rWorld.o, Normalize(rWorld.d), rWorld.tMax * rWorld.d.Length()));
		// Compute $[\tmin, \tmax]$ interval of _ray_'s overlap with medium bounds
		const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
		Float tMin, tMax;
		if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(1.f);

		// Perform ratio tracking to estimate the transmittance value
		Float Tr = 1, t = tMin;
		while (true) {
			++nTrSteps;
			t -= std::log(1 - sampler.Get1D()) * invMaxDensity / sigma_t[0];
			if (t >= tMax) break;
			Float density = Density(ray(t));
			Tr *= 1 - std::max((Float)0, density * invMaxDensity);
			// Added after book publication: when transmittance gets low,
			// start applying Russian roulette to terminate sampling.
			const Float rrThreshold = .1;
			if (Tr < rrThreshold) {
				Float q = std::max((Float).05, 1 - Tr);
				if (sampler.Get1D() < q) return 0;
				Tr /= 1 - q;
			}
		}
		return Spectrum(Tr);
	}

}  // namespace pbrt
