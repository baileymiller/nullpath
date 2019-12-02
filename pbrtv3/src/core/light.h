
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

#ifndef PBRT_CORE_LIGHT_H
#define PBRT_CORE_LIGHT_H

// core/light.h*
#include "pbrt.h"
#include "memory.h"
#include "interaction.h"

namespace pbrt {

// LightFlags Declarations
enum class LightFlags : int {
    DeltaPosition = 1,
    DeltaDirection = 2,
    Area = 4,
    Infinite = 8
};

inline bool IsDeltaLight(int flags) {
    return flags & (int)LightFlags::DeltaPosition ||
           flags & (int)LightFlags::DeltaDirection;
}

// Light Declarations
class Light {
  public:
    // Light Interface
    virtual ~Light();
    Light(int flags, const Transform &LightToWorld,
          const MediumInterface &mediumInterface, int nSamples = 1);
    virtual Spectrum Sample_Li(const Interaction &ref, const Point2f &u,
                               Vector3f *wi, Float *pdf,
                               VisibilityTester *vis) const = 0;
    virtual Spectrum Power() const = 0;
    virtual void Preprocess(const Scene &scene) {}
    virtual Spectrum Le(const RayDifferential &r) const;
    virtual Float Pdf_Li(const Interaction &ref, const Vector3f &wi) const = 0;
    virtual Spectrum Sample_Le(const Point2f &u1, const Point2f &u2, Float time,
                               Ray *ray, Normal3f *nLight, Float *pdfPos,
                               Float *pdfDir) const = 0;
    virtual void Pdf_Le(const Ray &ray, const Normal3f &nLight, Float *pdfPos,
                        Float *pdfDir) const = 0;

    virtual bool Intersects(const Interaction &ref, const Vector3f &wi) const {
      return true;
    }

    // Light Public Data
    const int flags;
    const int nSamples;
    const MediumInterface mediumInterface;

  protected:
    // Light Protected Data
    const Transform LightToWorld, WorldToLight;
};

class VisibilityTester {
  public:
    VisibilityTester()
        : m_enableShadows(true)
        , m_transportMode(TransportMode::Radiance)
    {}

    // VisibilityTester Public Methods
    VisibilityTester(const Interaction &p0, const Interaction &p1, TransportMode transportMode = TransportMode::Radiance, bool enableShadows = true)
        : p0(p0), p1(p1), m_transportMode(transportMode), m_enableShadows(enableShadows) {}
    const Interaction &P0() const { return p0; }
    const Interaction &P1() const { return p1; }
	void setTransportMode(TransportMode mode) { m_transportMode = mode; }
	void setEnableShadows(bool enableShadows = true) { m_enableShadows = enableShadows; }
    bool Unoccluded(const Scene &scene) const;
	Spectrum Tr(const Scene &scene, Sampler &sampler, uint32_t flags = 0) const;

  private:
    Interaction p0, p1;
    TransportMode m_transportMode;
	bool m_enableShadows;
};

class AreaLight : public Light {
  public:
    // AreaLight Interface
    AreaLight(const Transform &LightToWorld, const MediumInterface &medium,
              int nSamples);
    virtual Spectrum L(const Interaction &intr, const Vector3f &w) const = 0;

    virtual Point3f Sample(const Point2f &u, Float *pdf) const = 0;

    virtual std::shared_ptr<Shape> GetShape() const = 0;
};

}  // namespace pbrt

#endif  // PBRT_CORE_LIGHT_H
