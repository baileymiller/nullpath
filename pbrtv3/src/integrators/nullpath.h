#ifndef PBRT_INTEGRATORS_NULLPATH_H
#define PBRT_INTEGRATORS_NULLPATH_H

// integrators/nullpath.h*
#include "pbrt.h"
#include "integrator.h"
#include "lightdistrib.h"
#include <sstream>

namespace pbrt {

// NullPathIntegrator Declarations
class NullPathIntegrator : public SamplerIntegrator {
  public:
    // VolPathIntegrator Public Methods
    NullPathIntegrator(int maxDepth, 
                      std::shared_ptr<const Camera> camera,
                      std::shared_ptr<Sampler> sampler,
                      const Bounds2i &pixelBounds,
                      Spectrum p_channel,
                      Float rrThreshold = 1,
                      const std::string &lightSampleStrategy = "spatial",
                      bool enableShadows = true,
											bool enableSurfaceNEE = false,
                      bool enableEA = false,
                      bool enableOnlyEA = false,
                      bool enableNEE = false, 
                      bool enableSplitNEE = false,
                      bool enableOnlyNEE = false,
                      bool enableLightSample = false,
                      bool enableEqualMajorant = false,
                      bool enableDeltaTrack = false,
                      bool singleScatter = false,
                      int viewNumScatterEvents = -1,
                      float constantEmission = 0.0f,
                      bool ignoreNoScatter = false,
                      bool enableBalanceHeuristic = true,
                      bool printPath = false,
                      bool printHits = false,
                      bool printEAPath = false,
                      bool showTwoStrategyMISWeight = false,
                      bool splitTwoStrategyMISValue = false,
                      bool showDirPdf = false)
        : SamplerIntegrator(camera, sampler, pixelBounds),
          maxDepth(maxDepth),
          m_p_channel(p_channel),
          rrThreshold(rrThreshold),
          lightSampleStrategy(lightSampleStrategy),
          m_enableShadows(enableShadows),
					m_enableSurfaceNEE(enableSurfaceNEE),
          m_enableEA(enableEA),
	  m_enableOnlyEA(enableOnlyEA),
          m_enableNEE(enableNEE),
          m_enableSplitNEE(enableSplitNEE),
          m_enableOnlyNEE(enableOnlyNEE),
          m_enableLightSample(enableLightSample),
          m_enableEqualMajorant(enableEqualMajorant),
          m_enableDeltaTrack(enableDeltaTrack),
          m_singleScatter(singleScatter),
	  m_viewNumScatterEvents(viewNumScatterEvents),
	  m_constantEmission(constantEmission),
          m_ignoreNoScatter(ignoreNoScatter),
          m_enableBalanceHeuristic(enableBalanceHeuristic),
          m_printPath(printPath),
          m_printHits(printHits),
          m_printEAPath(printEAPath),
          m_showTwoStrategyMISWeight(showTwoStrategyMISWeight),
          m_splitTwoStrategyMISValue(splitTwoStrategyMISValue),
          m_showDirPdf(showDirPdf) { }

    void Preprocess(const Scene &scene, Sampler &sampler);
    Spectrum Li(const RayDifferential &ray, const Scene &scene,
                Sampler &sampler, MemoryArena &arena, int depth) const;

  private:
    // NullPathIntegrator Private Data
    const int maxDepth;
    const Spectrum m_p_channel;
    const Float rrThreshold;
    const std::string lightSampleStrategy;
    bool m_enableShadows;
		bool m_enableSurfaceNEE;
    bool m_enableEA;
    bool m_enableOnlyEA;

    bool m_enableNEE;
    bool m_enableSplitNEE;
    bool m_enableOnlyNEE;
    bool m_enableLightSample;
    bool m_enableEqualMajorant;
    bool m_enableDeltaTrack;

    bool m_singleScatter;
    int m_viewNumScatterEvents;
    float m_constantEmission;
    bool m_ignoreNoScatter;
    bool m_enableBalanceHeuristic;
    bool m_printPath;
    bool m_printHits;
    bool m_printEAPath;
    bool m_showTwoStrategyMISWeight;
    bool m_splitTwoStrategyMISValue;
    bool m_showDirPdf;

    std::unique_ptr<LightDistribution> lightDistribution;

    Spectrum CalculateMISWeight(std::vector<double> P_OVER_F, 
                                Spectrum L = Spectrum(1.0f)) const;

    Spectrum CalculateTwoStrategyMISWeight(std::vector<double> P_OVER_F_1, 
					   std::vector<double> P_OVER_F_2,
					   Spectrum L = Spectrum(1.0f),
					   bool isStrat1 = false) const;

    int PickHeroWavelength(Sampler &sampler) const;

    void PrintPath(std::stringstream *path) const;
};


NullPathIntegrator *CreateNullPathIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera);

}  // namespace pbrt

#endif  // PBRT_INTEGRATORS_NULLPATH_H
