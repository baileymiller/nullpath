#ifndef PBRT_INTEGRATORS_DELTA_H
#define PBRT_INTEGRATORS_DELTA_H

// integrators/nullpath.h*
#include "pbrt.h"
#include "integrator.h"
#include "lightdistrib.h"

namespace pbrt {

// SpectralIntegrator Declarations
class DeltaIntegrator : public SamplerIntegrator {
  public:
    // SpectralIntegrator Public Methods
    DeltaIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
                      std::shared_ptr<Sampler> sampler,
                      const Bounds2i &pixelBounds, Float rrThreshold = 1,
                      const std::string &lightSampleStrategy = "spatial", 
                      bool enableShadows = true,
                      bool enableNEE = false,
                      bool singleScatter = false,
                      bool ignoreNoScatter = false)
        : SamplerIntegrator(camera, sampler, pixelBounds),
          maxDepth(maxDepth),
          rrThreshold(rrThreshold),
          lightSampleStrategy(lightSampleStrategy),
		  m_enableShadows(enableShadows),
          m_enableNEE(enableNEE),
          m_singleScatter(singleScatter),
          m_ignoreNoScatter(ignoreNoScatter) { }
    void Preprocess(const Scene &scene, Sampler &sampler);
    Spectrum Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const;

  private:
    // SpectralIntegrator Private Data
    const int maxDepth;
    const Float rrThreshold;
    const std::string lightSampleStrategy;
    bool m_enableShadows;
    bool m_enableNEE;
    bool m_singleScatter;
    bool m_ignoreNoScatter;
    std::unique_ptr<LightDistribution> lightDistribution;

    int PickSingleWavelength(Sampler &sampler) const;
};

DeltaIntegrator *CreateDeltaIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_DELTA_H
