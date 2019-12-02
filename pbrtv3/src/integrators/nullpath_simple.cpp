// integrators/nullpath_simple.cpp*
#include "integrators/nullpath_simple.h"
#include "bssrdf.h"
#include "camera.h"
#include "film.h"
#include "interaction.h"
#include "paramset.h"
#include "scene.h"
#include "stats.h"

namespace pbrt {

STAT_FLOAT_DISTRIBUTION("Integrator/Path length", pathLength);
STAT_COUNTER("Integrator/Volume interactions", volumeInteractions);
STAT_COUNTER("Integrator/Surface interactions", surfaceInteractions);

// NullPathSimpleIntegrator Method Definitions
void NullPathSimpleIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
    lightDistribution =
        CreateLightSampleDistribution(lightSampleStrategy, scene);
}

void
NullPathSimpleIntegrator::PrintPath(std::stringstream *path) const {
  if (m_printPath) std::cout << "\nPATH: " << path->str() << std::endl;
}

Spectrum NullPathSimpleIntegrator::Li(const RayDifferential &r, const Scene &scene,
                               Sampler &sampler, MemoryArena &arena,
                               int depth) const {
  ProfilePhase p(Prof::SamplerIntegratorLi);
  Spectrum L(0.0f);

  // P_r / F_r   P_g / F_r   P_b / F_r   
  // P_r / F_g   P_g / F_g   P_b / F_g   
  // P_r / F_b   P_g / F_b   P_b / F_b   
  std::vector<double> P_OVER_F = std::vector<double>(9, 1.0f);

  int hero = 0.0f;

  RayDifferential ray(r);
  bool specularBounce = false;
  bool lastScatterWasSurface = false;
  int bounces = 0;
  int numScatterEvents = 0;
  int realPathLength = 0;
  bool lastEventWasReal = true;
  Ray lastRealRay;

  // Added after book publication: etaScale tracks the accumulated effect
  // of radiance scaling due to rays passing through refractive
  // boundaries (see the derivation on p. 527 of the third edition). We
  // track this value in order to remove it from beta when we apply
  // Russian roulette; this is worthwhile, since it lets us sometimes
  // avoid terminating refracted rays that are about to be refracted back
  // out of a medium and thus have their beta value increased.
  Float etaScale = 1;
  bool scattered = false;
  for (bounces = 0;; bounces++) {

      // Intersect _ray_ with scene and store intersection in _isect_
      SurfaceInteraction isect;
      bool foundIntersection = scene.Intersect(ray, &isect);

      // Keeps track of the last ray traced from a real collision or boundary.
      if (lastEventWasReal) lastRealRay = ray;

      // Sample the participating medium, if present
      MediumInteraction mi;
      Spectrum Tr(1.0f), majorants;
      if (ray.medium) {
        // Sample a distance.
        Tr = ray.medium->SampleChannel(ray, sampler, arena, &mi, hero);
       
        // update weights with the transmittance and pdf.
        majorants = ray.medium->GetMajorant();
        Spectrum freeFlightPdf = mi.IsValid() ? majorants * Tr : Tr; 
        UpdateWeights(P_OVER_F, Tr, freeFlightPdf);
      }
    
      // Handle an interaction with a medium or a surface
      if (mi.IsValid()) {
        ++volumeInteractions;

        // Get the scattering and absorption coefficients. 
        Spectrum scatter, absorption, null;
        ray.medium->GetCoefficients(mi.p, scatter, absorption, null);
          
        Float pScatter = scatter[hero] / majorants[hero];
        Float pAbsorb = absorption[hero] / majorants[hero];
        
        // Sample a vertex type. 
        Vector3f wi;
        float xi = sampler.Get1D();
        if (xi < pAbsorb) { 
          return L;
        } else if (xi < 1.0f - pScatter) { 
          lastEventWasReal = false;

          // Update path to account for null event. 
          UpdateWeights(P_OVER_F, null, null / majorants);

          wi = ray.d;
        } else {
          lastEventWasReal = true;
          lastScatterWasSurface = false;
          numScatterEvents++;
          realPathLength++;

          bool reachedMax = realPathLength >= maxDepth && maxDepth > 0;
          if (reachedMax) return L;
         
          // Update path to account for scattering. 
          UpdateWeights(P_OVER_F, scatter, scatter / majorants);
          
          // Sample a direction.
          Vector3f wo = -ray.d;
          mi.phase->Sample_p(wo, &wi, sampler.Get2D());
          
          // UpdateWeights(P_OVER_F, scatterFunc, scatterPdf); phaseFunc == phasePdf
          scattered = true;
        }
      }
      ray = mi.SpawnRay(wi); 
    } else {
      ++surfaceInteractions; 
      ++realPathLength;
        // If intersection goes out of scene, look up the environment lights.
      if (!foundIntersection) {
        for(auto light : scene.infiniteLights) {
            L += CalculateMISWeight(P_OVER_F, light->Le(ray));
        }
        return L;
      }
      
      if (!isect.Le(-ray.d).IsBlack()) {
        L += CalculateMISWeight(P_OVER_F, isect.Le(-ray.d)); 
        return L;
      }
      
      // Compute scattering functions.
      isect.ComputeScatteringFunctions(ray, arena, true);

      // If no bsdf at intersection, medium boundary so continue. 
      if (!isect.bsdf) {
        --realPathLength;
        ray = isect.SpawnRay(ray.d);
        continue;
     } 

      lastScatterWasSurface = true;

      bool reachedMax = realPathLength >= maxDepth && maxDepth > 0;
      if (reachedMax) return L;

      // Sample BSDF to get new path direction
      Vector3f wo = -ray.d, wi;
      Float bsdfPdf;
      BxDFType flags;
      Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &bsdfPdf,
                                        BSDF_ALL, &flags);
      // If the path contribution function is 0 or pdf is 0 return Spectrum(0.0f).
      if (f.IsBlack() || bsdfPdf == 0.0f) {
        if (m_printPath) PrintPath(&symbPath);
        return L;
      }
      specularBounce = (flags & BSDF_SPECULAR) != 0;

      if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)) {
            Float eta = isect.bsdf->eta;
            // Update the term that tracks radiance scaling for refraction
            // depending on whether the ray is entering or leaving the
            // medium.
            etaScale *=
                (Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
      }

      Spectrum bsdfFunc = f * AbsDot(wi, isect.shading.n);
      UpdateWeights(P_OVER_F, bsdfFunc, bsdfPdf);

      ray = isect.SpawnRay(wi);
    }
  }
}

NullPathIntegrator *CreateNullPathIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera) {
    int maxDepth = params.FindOneInt("maxdepth", 0);
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

    return new NullPathIntegrator(maxDepth, 
                                camera, 
                                sampler, 
                                pixelBounds,
                                p_channel, 
                                rrThreshold, 
                                lightStrategy);
}

}  // namespace pbrt
