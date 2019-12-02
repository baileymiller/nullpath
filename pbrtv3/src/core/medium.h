
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

#ifndef PBRT_CORE_MEDIUM_H
#define PBRT_CORE_MEDIUM_H

// core/medium.h*
#include "pbrt.h"
#include "geometry.h"
#include "material.h"
#include "spectrum.h"
#include <memory>

namespace pbrt {

// different non-exponential modes
// these modes correspond to the following contributions
// surface->surface, medium->surface, surface->medium, medium->medium
enum TransmittanceMode
{
  NONEXP_SSC_MSC_SMC_MMC, // the original formulation
  NONEXP_SST_MST_SMC_MMC, // non-reciprocal energy conserving
  NONEXP_SST_MSC_SMT_MMC, // the counterpart of non-reciprocal to make it energy conserving - we dont generally use it as such
  NONEXP_SST_MSC_SMC_MMC, // Reciprocal collision
  NONEXP_SST_MST_SMT_MMC, // Reciprocal tracklength
  NONEXP_SST_MST_SMT_MMT, // Reciprocal exclusive
  MODE_UNKNOWN      // used to indicate an unknown flag
};

enum TransmittanceFunctionType
{
  TRANSMITTANCE_EXPONENTIAL,
  TRANSMITTANCE_DAVIS,
  TRANSMITTANCE_DAVIS_2PARAM,
  TRANSMITTANCE_QUADRATIC_RAMP,
  TRANSMITTANCE_LINEAR_RAMP,
  TRANSMITTANCE_ERLANG2,
  TRANSMITTANCE_UNKNOWN
};

// Media Declarations
class PhaseFunction {
  public:
  // PhaseFunction Interface
  virtual ~PhaseFunction();
  virtual Float p(const Vector3f &wo, const Vector3f &wi) const = 0;
  virtual Float Sample_p(const Vector3f &wo, Vector3f *wi,
         const Point2f &u) const = 0;
  virtual std::string ToString() const = 0;
};

inline std::ostream &operator<<(std::ostream &os, const PhaseFunction &p) {
  os << p.ToString();
  return os;
}

bool GetMediumScatteringProperties(const std::string &name, Spectrum *sigma_a,
           Spectrum *sigma_s);

// Media Inline Functions
inline Float PhaseHG(Float cosTheta, Float g) {
  Float denom = 1 + g * g + 2 * g * cosTheta;
  return Inv4Pi * (1 - g * g) / (denom * std::sqrt(denom));
}

// Medium Declarations
class Medium {
public:
  // Medium Interface
  virtual ~Medium() {}
  virtual Spectrum Tr(const Ray &ray, Sampler &sampler, uint32_t flags, TransportMode mode = TransportMode::Radiance) const = 0;
  virtual Spectrum Tr(const RayDifferential& ray, Sampler& sampler, uint32_t flags, TransportMode mode = TransportMode::Radiance) const
  {
  std::cout << "Error: Tr() with RayDifferentials is not implemented by default." << std::endl;
  return Spectrum(1.0f);
  }
  virtual Spectrum Sample(const Ray &ray, Sampler &sampler, MemoryArena &arena, MediumInteraction *mi, uint32_t flags = 0xFFFFFFFF, TransportMode mode = TransportMode::Radiance) const = 0;
  
  virtual Spectrum SampleChannel(const Ray &ray, Sampler &sampler, MemoryArena &arena, MediumInteraction *mi, int channel, uint32_t flags = 0xFFFFFFFF, TransportMode mode = TransportMode::Radiance) const {
    std::cout << "ERROR: SampleChannel is not implemented by default." << std::endl;
    return Spectrum(1.0f);
  };

  virtual Spectrum Sample(const RayDifferential& ray, Sampler &sampler, MemoryArena &arena, MediumInteraction *mi, uint32_t flags = 0xFFFFFFFF, TransportMode mode = TransportMode::Radiance) const
  {
  std::cout << "ERROR: Sample with RayDifferential is not implemented by default." << std::endl;
  return Spectrum(1.0f);
  } 

  virtual Spectrum GetAbsorption(const Point3f &p) const {
    std::cout << "ERROR: GetAbsorption is not implemented by default." << std::endl;
    return Spectrum(1.0f);
  };

  virtual Spectrum GetScattering(const Point3f &p) const {
    std::cout << "ERROR: GetScattering is not implemented by default."<< std::endl;
    return Spectrum(1.0f);
  };

  virtual void GetCoefficients(const Point3f &p, Spectrum &s, Spectrum &a, Spectrum &n) const {
    std::cout << "ERROR: GetCoefficients is not implemented by default."<< std::endl;
  };

  virtual Spectrum GetMajorant() const {
    std::cout << "ERROR: GetMajorant is not implemented by default."<< std::endl;
    return Spectrum(1.0f);
  };

  virtual Float GetG() const { 
    std::cout << "ERROR: GetG is not implemented by default." << std::endl;
    return 0.0f;
  };

  virtual bool SetMinAndMaxPosition(Ray &ray) const {
    std::cout << "ERROR: SetMinAndMaxPosition is not implemented by default" << std::endl;
    return true;
  };

};

// HenyeyGreenstein Declarations
class HenyeyGreenstein : public PhaseFunction {
  public:
  // HenyeyGreenstein Public Methods
  HenyeyGreenstein(Float g) : g(g) {}
  Float p(const Vector3f &wo, const Vector3f &wi) const;
  Float Sample_p(const Vector3f &wo, Vector3f *wi,
       const Point2f &sample) const;
  std::string ToString() const {
  return StringPrintf("[ HenyeyGreenstein g: %f ]", g);
  }

  private:
  const Float g;
};

// MediumInterface Declarations
struct MediumInterface {
  MediumInterface() : inside(nullptr), outside(nullptr) {}
  // MediumInterface Public Methods
  MediumInterface(const Medium *medium) : inside(medium), outside(medium) {}
  MediumInterface(const Medium *inside, const Medium *outside)
  : inside(inside), outside(outside) {}
  bool IsMediumTransition() const { return inside != outside; }
  const Medium *inside, *outside;
};

}  // namespace pbrt

#endif  // PBRT_CORE_MEDIUM_H
