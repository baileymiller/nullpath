#ifndef PBRT_INTEGRATORS_SPECTRAL_H
#define PBRT_INTEGRATORS_SPECTRAL_H

// integrators/nullpath.h*
#include "pbrt.h"
#include "integrator.h"
#include "lightdistrib.h"

namespace pbrt {

// SpectralIntegrator Declarations
class SpectralIntegrator : public SamplerIntegrator {
  public:
    // SpectralIntegrator Public Methods
    SpectralIntegrator(int maxDepth, std::shared_ptr<const Camera> camera,
                      std::shared_ptr<Sampler> sampler,
                      const Bounds2i &pixelBounds, Float rrThreshold = 1,
                      const std::string &lightSampleStrategy = "spatial",
                      bool enableShadows = true,
                      bool enableSurfaceNEE = false,
                      bool enableNaiveDelta = false,
                      bool enableNEE = false,
											bool enableSplitNEE = false,
                      bool enableOnlyNEE = false,
                      bool singleScatter = false,
                      int viewNumScatterEvents = -1,
                      bool ignoreNoScatter = false)
        : SamplerIntegrator(camera, sampler, pixelBounds),
          maxDepth(maxDepth),
          rrThreshold(rrThreshold),
          lightSampleStrategy(lightSampleStrategy),
		      m_enableShadows(enableShadows),
          m_enableSurfaceNEE(enableSurfaceNEE),
          m_enableNaiveDelta(enableNaiveDelta),
          m_enableNEE(enableNEE),
          m_enableOnlyNEE(enableOnlyNEE),
					m_enableSplitNEE(enableSplitNEE),
          m_singleScatter(singleScatter),
          m_viewNumScatterEvents(viewNumScatterEvents),
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
    bool m_enableSurfaceNEE;
    bool m_enableNaiveDelta;
    bool m_enableNEE;
		bool m_enableSplitNEE;
    bool m_enableOnlyNEE;
    bool m_singleScatter;
    int m_viewNumScatterEvents;
    bool m_ignoreNoScatter;
    std::unique_ptr<LightDistribution> lightDistribution;

		Spectrum CalculateMISWeight(Spectrum L, Spectrum W, Spectrum W_NEE) const;
};

SpectralIntegrator *CreateSpectralIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_SPECTRAL_H
