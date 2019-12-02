
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MEDIA_GRID_H
#define PBRT_MEDIA_GRID_H

// media/grid.h*
#include "medium.h"
#include "transform.h"
#include "stats.h"

namespace pbrt {

	STAT_MEMORY_COUNTER("Memory/Mod Volume density grid", densityBytes);
        enum InterpType {
          HEIGHT_Y,
          UNKNOWN
        };

	// GridDensityMedium Declarations
	class GridDensityMedium : public Medium {
	public:
                
		// GridDensityMedium Public Methods
		GridDensityMedium(const Spectrum &sigma_a, const Spectrum &sigma_s, Float g,
			int nx, int ny, int nz, const Transform &mediumToWorld,
			const Float *d,
                        bool enableDualDensity = false,
                        const Spectrum &sigma_n = Spectrum(0.0f),
                        const Spectrum &sigma_a_2 = Spectrum(0.0f),
                        const Spectrum &sigma_s_2 = Spectrum(0.0f),
                        const Spectrum &sigma_n_2 = Spectrum(0.0f),
                        std::string interp = "")
			: sigma_a(sigma_a),
			sigma_s(sigma_s),
                        sigma_n(sigma_n),
                        sigma_a_2(sigma_a_2),
                        sigma_s_2(sigma_s_2),
                        sigma_n_2(sigma_n_2),
                        m_interpType(StringToInterpType(interp)),
                        m_enableDualDensity(enableDualDensity),
			g(g),
			nx(nx),
			ny(ny),
			nz(nz),
      MediumToWorld(mediumToWorld),
			WorldToMedium(Inverse(mediumToWorld)),
			density(new Float[nx * ny * nz]) {
			densityBytes += nx * ny * nz * sizeof(Float);
			memcpy((Float *)density.get(), d, sizeof(Float) * nx * ny * nz);
			// Precompute values for Monte Carlo sampling of _GridDensityMedium_
			sigma_t = sigma_a + sigma_s + sigma_n;
      Spectrum sigma_t_2 = sigma_a_2 + sigma_s_2 + sigma_n_2;
      if (enableDualDensity) {
        for (int i = 0; i < Spectrum::nSamples; i++) {
          sigma_t[i] = std::max(sigma_t[i], sigma_t_2[i]);
        }
      }
			Float maxDensity = 0;
			for (int i = 0; i < nx * ny * nz; ++i)
				maxDensity = std::max(maxDensity, density[i]);
			invMaxDensity = 1 / maxDensity;
		}

		Float Density(const Point3f &p) const;
		Float D(const Point3i &p) const {
			Bounds3i sampleBounds(Point3i(0, 0, 0), Point3i(nx, ny, nz));
			if (!InsideExclusive(p, sampleBounds)) return 0;
			return density[(p.z * ny + p.y) * nx + p.x];
		}
		Spectrum Sample(const Ray &ray, Sampler &sampler, MemoryArena &arena, MediumInteraction *mi) const;
                Spectrum SampleChannel(const Ray &rWorld, Sampler &sampler, MemoryArena &arena,    
                                      MediumInteraction *mi, int channel) const;

                Spectrum GetAbsorption(const Point3f &p) const;
                Spectrum GetScattering(const Point3f &p) const;
                void GetCoefficients(const Point3f &p, Spectrum &s, Spectrum &a, Spectrum &n) const;

                Spectrum GetMajorant() const;
                Spectrum Interpolate(Spectrum x, Spectrum y, Point3f p) const;

                InterpType StringToInterpType(std::string type) {
                  if (type == "height_y") {
                    return InterpType::HEIGHT_Y;
                  } else {
                    return InterpType::UNKNOWN;
                  }
                };


		Spectrum Tr(const Ray &ray, Sampler &sampler) const;

    Float GetG() const { return g; };

	private:
		// GridDensityMedium Private Data
		const Spectrum sigma_a, sigma_s, sigma_n;
                const Spectrum sigma_a_2, sigma_s_2, sigma_n_2;
                bool m_enableDualDensity;
                const InterpType m_interpType;

                Spectrum sigma_t; 

		const Float g;
		const int nx, ny, nz;

                const Transform MediumToWorld;
		const Transform WorldToMedium;

		std::unique_ptr<Float[]> density;
		Float invMaxDensity;
	};

}  // namespace pbrt

#endif  // PBRT_MEDIA_GRID_H
