#ifndef PBRT_INTEGRATORS_NULLPATH_SIMPLE_H
#define PBRT_INTEGRATORS_NULLPATH_SIMPLE_H

// integrators/nullpath.h*
#include "pbrt.h"
#include "integrator.h"
#include "lightdistrib.h"
#include <sstream>

namespace pbrt {

// NullPathSimpleIntegrator Declarations
class NullPathSimpleIntegrator : public SamplerIntegrator {
  public:
    // NullPathSimpleIntegrator Public Methods
    NullPathSimpleIntegrator(int maxDepth, 
                      std::shared_ptr<const Camera> camera,
                      std::shared_ptr<Sampler> sampler,
                      const Bounds2i &pixelBounds,
                      Float rrThreshold = 1,
                      const std::string &lightSampleStrategy = "spatial")
        : SamplerIntegrator(camera, sampler, pixelBounds),
          maxDepth(maxDepth),
          rrThreshold(rrThreshold),
          lightSampleStrategy(lightSampleStrategy) {}

    void Preprocess(const Scene &scene, Sampler &sampler);
    Spectrum Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const;

  private:
    // NullPathSimpleIntegrator Private Data
    const int maxDepth;
    const Float rrThreshold;
    const std::string lightSampleStrategy;
    std::unique_ptr<LightDistribution> lightDistribution;
};


NullPathSimpleIntegrator *CreateNullPathSimpleIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_NULLPATH_SIMPLE_H
