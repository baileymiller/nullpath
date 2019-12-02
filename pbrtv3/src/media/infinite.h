
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

#ifndef PBRT_MEDIA_INFINITE_H
#define PBRT_MEDIA_INFINITE_H

// media/infinite.h*
#include "medium.h"

namespace pbrt {

// HomogeneousMedium Declarations
class InfiniteMedium : public Medium {
  public:
    // InfiniteMedium Public Methods
    InfiniteMedium(const Spectrum &sigma_a, const Spectrum &sigma_s, 
                      const Spectrum &sigma_n,
                      Float g);
  
  Spectrum Tr(const Ray &ray, Sampler &sampler, uint32_t flags = 0xFFFFFFFF, 
              TransportMode mode = TransportMode::Radiance) const;

  Spectrum Tr(const RayDifferential& ray, Sampler& sampler, 
              uint32_t flags = 0xFFFFFFFF, 
              TransportMode mode = TransportMode::Radiance) const;
	
  Spectrum Sample(const Ray &ray, Sampler &sampler, 
                  MemoryArena &arena, MediumInteraction *mi, 
                  uint32_t flags = 0xFFFFFFFF, 
                  TransportMode mode = TransportMode::Radiance) const;
  	
  Spectrum SampleChannel(const Ray &ray, Sampler &sampler, MemoryArena &arena, 
                         MediumInteraction *mi, 
                         int channel, uint32_t flags = 0xFFFFFFFF, 
                         TransportMode mode = TransportMode::Radiance) const;
	
  Spectrum Sample(const RayDifferential &ray, Sampler &sampler, 
                  MemoryArena &arena, MediumInteraction *mi, 
                  uint32_t flags = 0xFFFFFFFF, 
                  TransportMode mode = TransportMode::Radiance) const;

  Spectrum GetAbsorption(const Point3f &p) const;
  Spectrum GetScattering(const Point3f &p) const;
  Spectrum GetMajorant() const;

  private:
    // InfiniteMedium Private Data
    const Spectrum sigma_a, sigma_s, sigma_t;
    const Float g;
};

}  // namespace pbrt

#endif  // PBRT_MEDIA_INFINITE_H
