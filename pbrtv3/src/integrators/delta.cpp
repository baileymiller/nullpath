// integrators/spectral.cpp*
#include "integrators/delta.h"
#include "bssrdf.h"
#include "camera.h"
#include "film.h"
#include "interaction.h"
#include "paramset.h"
#include "scene.h"
#include "stats.h"
#include "sampling.h"

namespace pbrt {

STAT_FLOAT_DISTRIBUTION("Integrator/Path length", pathLength);
STAT_COUNTER("Integrator/Volume interactions", volumeInteractions);
STAT_COUNTER("Integrator/Surface interactions", surfaceInteractions);

// SpectralIntegrator Method Definitions
void DeltaIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
    lightDistribution =
        CreateLightSampleDistribution(lightSampleStrategy, scene);
}

int DeltaIntegrator::PickSingleWavelength(Sampler &sampler) const {
  Float rand = sampler.Get1D();
  // Assume there are three channels.
  return (int)(rand * 3);
}

Spectrum DeltaIntegrator::Li(const RayDifferential &r, const Scene &scene,
                               Sampler &sampler, MemoryArena &arena,
                               int depth) const {
    ProfilePhase p(Prof::SamplerIntegratorLi);
    Spectrum L(0.0f);
    Spectrum W(0.0f);

    RayDifferential ray(r);
    bool specularBounce = false;
    bool scattered = false;
    int numScatterEvents = 0;
    int bounces = 0; 
    
    int sampleChannel = PickSingleWavelength(sampler);
    W[sampleChannel] = 3.0f; // boost weight to correspond to prob. of picking channel.

    Float etaScale = 1;
    for(bounces = 0;; bounces++) {
      // Intersect _ray_ with scene and store intersection in _isect_
      SurfaceInteraction isect;
      bool foundIntersection = scene.Intersect(ray, &isect);

      // Sample the participating medium, if present
      MediumInteraction mi;
      if (ray.medium) ray.medium->SampleChannel(ray, sampler, arena, &mi, sampleChannel);
      
      // Handle an interaction with a medium or a surface
      if (mi.IsValid()) {
        ++volumeInteractions;

        Spectrum majorants = ray.medium->GetMajorant();
        Spectrum scatter, absorption, null;
        ray.medium->GetCoefficients(mi.p, scatter, absorption, null);

        Float pScatter = scatter[sampleChannel] / majorants[sampleChannel]; 
        Float pAbsorb = absorption[sampleChannel] / majorants[sampleChannel];
      
        Vector3f wi;
        float xi = sampler.Get1D();
        if (xi < pAbsorb) return Spectrum(0.0f);
        else if (xi < 1.0f - pScatter) wi = ray.d;
        else {
          if (scattered && m_singleScatter) break;
        
          // Generate a new sampling direction.
          Vector3f wo = -ray.d;
          if (m_enableNEE) {
            bool sampledLight;
            Float weight = MISPhaseLightDirection(scene, mi, sampler, wo, &wi, sampledLight);
            W *= mi.phase->p(wo, wi) / weight; 
          } else { 
            mi.phase->Sample_p(wo, &wi, sampler.Get2D());
          }

          numScatterEvents++;
          scattered = true;
        }
        ray = mi.SpawnRay(wi);
      } else {
        ++surfaceInteractions;
        
        // Terminate path if ray escaped or _maxDepth_ was reached
        if (foundIntersection)  {
          L += FlagFilter(W * isect.Le(-ray.d), 
                          m_ignoreNoScatter,
													-1,
                          numScatterEvents);
        } else {
          for(const auto &light : scene.infiniteLights) {
            L += FlagFilter(W * light->Le(ray),
                            m_ignoreNoScatter,
														-1,
                            numScatterEvents);
          }
        }
        
        if (!foundIntersection || !L.IsBlack()) break;

       // Compute scattering functions and skip over medium boundaries
        isect.ComputeScatteringFunctions(ray, arena, true);
        if (!isect.bsdf) {
          ray = isect.SpawnRay(ray.d);
          bounces--;
          continue;
        }

        // Sample BSDF to get new path direction
        Vector3f wo = -ray.d, wi;
        Float pdf;
        BxDFType flags;
        Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf,
                                          BSDF_ALL, &flags);
        if (f.IsBlack() || pdf == 0.f) return Spectrum(0.0f);
        W *= f * AbsDot(wi, isect.shading.n) / pdf;
        specularBounce = (flags & BSDF_SPECULAR) != 0;
        if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)) {
            Float eta = isect.bsdf->eta;
            // Update the term that tracks radiance scaling for refraction
            // depending on whether the ray is entering or leaving the
            // medium.
            etaScale *=
                (Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
          }
          ray = isect.SpawnRay(wi);
      }
      if (W.IsBlack()) break;  
    }
    return L;
}

DeltaIntegrator *CreateDeltaIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    bool enableShadows = params.FindOneBool("enable-shadows", true);
	if (!enableShadows)
		std::cout << "\tINFO: Shadows were disabled for the renderer.!!" << std::endl;
    int np;
    const int *pb = params.FindInt("pixelbounds", &np);
    Bounds2i pixelBounds = camera->film->GetSampleBounds();
    if (pb) {
        if (np != 4)
            Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                  np);
        else {
            pixelBounds = Intersect(pixelBounds,
                                    Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
            if (pixelBounds.Area() == 0)
                Error("Degenerate \"pixelbounds\" specified.");
        }
    }
    Float rrThreshold = params.FindOneFloat("rrthreshold", 1.);
    std::string lightStrategy =
        params.FindOneString("lightsamplestrategy", "spatial");
    
    bool enableNEE = params.FindOneBool("enable-NEE", false);
    bool singleScatter = params.FindOneBool("singleScatter", false);
    bool ignoreNoScatter = params.FindOneBool("ignoreNoScatter", false);

    return new DeltaIntegrator(maxDepth, camera, sampler, pixelBounds,
                                 rrThreshold, lightStrategy, enableShadows,
                                 enableNEE, singleScatter, ignoreNoScatter);
}

}  // namespace pbrt
