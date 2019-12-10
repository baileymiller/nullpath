#ifndef PBRT_INTEGRATORS_SPECTRAL_MIS_H
#define PBRT_INTEGRATORS_SPECTRAL_MIS_H

// integrators/spectralmis.h
#include "pbrt.h"
#include "integrator.h"
#include "lightdistrib.h"
#include <sstream>

namespace pbrt {

// SpectralMISIntegrator Declarations
class SpectralMISIntegrator : public SamplerIntegrator {
  public:
    // SpectralMISIntegrator Public Methods
    SpectralMISIntegrator(
      int maxDepth, 
      std::shared_ptr<const Camera> camera,
      std::shared_ptr<Sampler> sampler,
      const Bounds2i &pixelBounds,
      Float rrThreshold = 1,
      const std::string &lightSampleStrategy = "spatial",
      bool enableSingleScatter = true,
      bool enableRatioTracking = true
    )
    : SamplerIntegrator(camera, sampler, pixelBounds)
    , maxDepth(maxDepth)
    , rrThreshold(rrThreshold)
    , lightSampleStrategy(lightSampleStrategy)
    , enableSingleScatter(enableSingleScatter)
    , enableRatioTracking(enableRatioTracking) {}

    void Preprocess(const Scene &scene, Sampler &sampler);
    Spectrum Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const;

  private:
    // SpectralMISIntegrator Private Data
    const int maxDepth;
    const Float rrThreshold;
    const std::string lightSampleStrategy;
    std::unique_ptr<LightDistribution> lightDistribution;
    const bool enableSingleScatter;
    const bool enableRatioTracking;
};


SpectralMISIntegrator *CreateSpectralMISIntegrator(
  const ParamSet &params,
  std::shared_ptr<Sampler> sampler,
  std::shared_ptr<const Camera> camera
);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_SPECTRAL_MIS_H
