
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

#ifndef PBRT_CORE_INTEGRATOR_H
#define PBRT_CORE_INTEGRATOR_H

// core/integrator.h*
#include "pbrt.h"
#include "primitive.h"
#include "spectrum.h"
#include "light.h"
#include "reflection.h"
#include "sampler.h"
#include "material.h"

namespace pbrt {

// Integrator Declarations
class Integrator {
  public:
    // Integrator Interface
    virtual ~Integrator();
    virtual void Render(const Scene &scene) = 0;
};

Spectrum UniformSampleAllLights(const Interaction &it, const Scene &scene,
                                MemoryArena &arena, Sampler &sampler,
                                const std::vector<int> &nLightSamples,
                                bool handleMedia = false,
	                            bool enableShadows = true);
Spectrum UniformSampleOneLight(const Interaction &it, const Scene &scene,
                               MemoryArena &arena, Sampler &sampler,
                               bool handleMedia = false,
                               const Distribution1D *lightDistrib = nullptr,
							   bool enableShadows = true);
Spectrum EstimateDirect(const Interaction &it, const Point2f &uShading,
                        const Light &light, const Point2f &uLight,
                        const Scene &scene, Sampler &sampler,
                        MemoryArena &arena, bool handleMedia = false, bool enableShadows = true,
                        bool specular = false);

std::unique_ptr<Distribution1D> ComputeLightPowerDistribution(
    const Scene &scene);

// Common functionality//

inline Float SQR(Float x) { return x * x; };

MediumInteraction SampleDistEquiangular(const Scene &scene, 
                                  MemoryArena &arena,
                                  Point3f lightPos,
                                  Ray ray,
                                  Float aRaySurfaceDist,
                                  Sampler &sampler,
                                  float &pdf);

Float DistEquiangularPDF(const Scene &scene, 
    Point3f lightPos,
                          Ray ray,
                          float aRaySurfaceDist,
                          float distance);

Float LightPdf(const Scene &scene, 
               SurfaceInteraction isect,
               Vector3f *wi);

Float LightPdf(const Scene &scene, 
               MediumInteraction mi,
               Vector3f *wi);

Float EALightPdf(const Scene &scene, 
                 std::shared_ptr<AreaLight> light,
                 SurfaceInteraction lightIsect,
                 Ray toScatter, 
                 Interaction scatterPoint);

double SampleLightDirection(
    const Scene &scene, MediumInteraction mi,
    Sampler &sampler, Vector3f *wi);

double MISPhaseLightDirection(
    const Scene &scene, MediumInteraction mi,
    Sampler &sampler, Vector3f wo, Vector3f *wi,
    bool &sampledLight);

Spectrum SampleSpectralNextEvent(
    const Scene &scene, SurfaceInteraction isect, 
    Sampler &sampler, Ray ray, MemoryArena &arena, 
    Spectrum &W, Spectrum &W_NEE, int channel = -1);

Spectrum SampleSpectralNextEvent(
    const Scene &scene, MediumInteraction mi, 
    Sampler &sampler, Ray ray, MemoryArena &arena, 
    Spectrum &W, Spectrum &W_NEE, int channel = -1);

Spectrum SampleNextEvent(
    const Scene &scene, SurfaceInteraction isect,
    Sampler &sampler, Ray ray, int hero, MemoryArena &arena, 
    std::vector<double> &P_OVER_F, std::vector<double> &P_OVER_F_NEE,
    bool equalizeMajorant = false,
    bool enableDeltaTrack = false);

Spectrum SampleNextEvent(
    const Scene &scene, MediumInteraction mi,
    Sampler &sampler, Ray ray, int hero, MemoryArena &arena, 
    std::vector<double> &P_OVER_F, std::vector<double> &P_OVER_F_NEE,
    bool equalizeMajorant = false,
    bool enableDeltaTrack = false);

bool CalculateNullPath(const Scene &scene, 
                       Sampler &sampler, 
                       MemoryArena &arena,
                       int hero,
                       Interaction &start,
                       Interaction &stop,
                       std::vector<double> &P_OVER_F,
                       std::vector<double> &P_OVER_F_FORCED_NULL,
                       std::string &path,
                       bool equalizeMajorant = false,
                       bool printPath = false);

Spectrum SampleNEE(const Scene &scene,
    Sampler &sampler, MediumInteraction mi, Ray ray,
    Ray lastRealRay,
    int hero,
    MemoryArena &arena, 
    std::vector<double> &P_OVER_F_NEE, std::vector<double> &P_OVER_F_EA,
    Spectrum scatterAtHitPoint,
    bool equalizeMajorant = false, bool printPath = false);

Spectrum SampleEA(
		const Scene &scene, Sampler &sampler, 
    Ray ray, 
    int hero, MemoryArena &arena,
		std::vector<double> &P_OVER_F, std::vector<double> &P_OVER_F_EA,
    bool equalizeMajorant = false,
    bool printPath = false);

Spectrum HomogenizeMajorants(Spectrum &majorants, int &maxChannel);

Spectrum SampleHomogenizedMajorants(const Ray &rWorld, Sampler &sampler, 
                                    MemoryArena &arena, MediumInteraction *mi, 
                                    int channel); 

void UpdateWeights(std::vector<double> &P_OVER_F, Spectrum F, Spectrum P);

std::vector<double> SetWeights(std::vector<double> &vec, double val);

std::vector<double> SetWeights(std::vector<double> &vec, std::vector<double> val);

std::vector<double> MultiplyVec(std::vector<double> &vec, double val);

std::vector<double> MultiplyVec(std::vector<double> &vec, Spectrum s);

std::vector<double> MultiplyVec(std::vector<double> &vec, 
                                std::vector<double> val);

std::vector<double> DivideVec(std::vector<double> &vec, double val);

std::vector<double> DivideVecVal(const std::vector<double> &vec, double val);

double AverageVec(const std::vector<double> &vec);

double MaxVal(const std::vector<double> &vec);

bool AllZerosVec(const std::vector<double> &vec);

std::vector<double> SpecToVec(Spectrum spec);

std::vector<double> SpecToVec(Spectrum spec, std::vector<double> &vec);

std::vector<double> FloatToVec(Float f, std::vector<double> &size);

Spectrum VecToSpec(const std::vector<double> &vec);

void PrintVec(const std::vector<double> &vec);


Spectrum FlagFilter(Spectrum s, 
                    bool ignoreNoScatter,
										int viewNumScatterEvents,
                    int numScatterEvents);

bool ExitEarly(bool singleScatter, int viewNumScatterEvents, int numScatterEvents);

// End Common functionality//

// SamplerIntegrator Declarations
class SamplerIntegrator : public Integrator {
  public:
    // SamplerIntegrator Public Methods
    SamplerIntegrator(std::shared_ptr<const Camera> camera,
                      std::shared_ptr<Sampler> sampler,
                      const Bounds2i &pixelBounds)
        : camera(camera), sampler(sampler), pixelBounds(pixelBounds) {}
    virtual void Preprocess(const Scene &scene, Sampler &sampler) {}
    void Render(const Scene &scene);
    virtual Spectrum Li(const RayDifferential &ray, const Scene &scene,
                        Sampler &sampler, MemoryArena &arena,
                        int depth = 0) const = 0;
    Spectrum SpecularReflect(const RayDifferential &ray,
                             const SurfaceInteraction &isect,
                             const Scene &scene, Sampler &sampler,
                             MemoryArena &arena, int depth) const;
    Spectrum SpecularTransmit(const RayDifferential &ray,
                              const SurfaceInteraction &isect,
                              const Scene &scene, Sampler &sampler,
                              MemoryArena &arena, int depth) const;

  protected:
    // SamplerIntegrator Protected Data
    std::shared_ptr<const Camera> camera;

  private:
    // SamplerIntegrator Private Data
    std::shared_ptr<Sampler> sampler;
    const Bounds2i pixelBounds;
};

}  // namespace pbrt

#endif  // PBRT_CORE_INTEGRATOR_H
