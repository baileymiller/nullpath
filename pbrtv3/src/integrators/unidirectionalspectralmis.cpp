// integrators/spectral_mis.cpp*
#include "integrators/unidirectionalspectralmis.h"
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

// UnidirectionalSpectralMISIntegrator Method Definitions
void UnidirectionalSpectralMISIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
    lightDistribution =
        CreateLightSampleDistribution(lightSampleStrategy, scene);
}

Spectrum UnidirectionalSpectralMISIntegrator::Li(const RayDifferential &r, const Scene &scene,
                                   Sampler &sampler, MemoryArena &arena,
                                   int depth) const {
  ProfilePhase p(Prof::SamplerIntegratorLi);
  RayDifferential ray(r);
  Spectrum L(0.0f);
  bool lastScatterWasSurface;
  bool specularBounce = false;
  int numScatterEvents = 0;

  // Choose a random color channel randomly. This color channel
  // will be used to make all sampling decisions for this path.
  int hero;
  float h = sampler.Get1D() * 3.0f;
  if (h < 1.0f) {
    hero = 0;
  } else if (h < 2.0f) {
    hero = 1;
  } else {
    hero = 2;
  }

  // Keep track of the real path length since we don't count
  // null collsions as bounces given a max depth (i.e. for single
  // scattering).
  int realPathLength = 0;

  // Store the ratios of the probability to function value.
  // We store 9 values since for each color channel we need to
  // know the probability given every possible hero channel.
  // These values are used to compute an MIS weight after 
  // sampling a path.
  //
  // P_r / F_r   P_g / F_r   P_b / F_r   
  // P_r / F_g   P_g / F_g   P_b / F_g   
  // P_r / F_b   P_g / F_b   P_b / F_b   
  std::vector<double> P_OVER_F = std::vector<double>(9, 1.0f);

  // Added after book publication: etaScale tracks the accumulated effect
  // of radiance scaling due to rays passing through refractive
  // boundaries (see the derivation on p. 527 of the third edition). We
  // track this value in order to remove it from beta when we apply
  // Russian roulette; this is worthwhile, since it lets us sometimes
  // avoid terminating refracted rays that are about to be refracted back
  // out of a medium and thus have their beta value increased.
  Float etaScale = 1;
  for (int bounces = 0;; bounces++) {
      SurfaceInteraction isect;
      bool foundIntersection = scene.Intersect(ray, &isect);

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
      
      // ========== SAMPLE VERTEX EVENT=============//
      Vector3f wi;
      float xi = sampler.Get1D();
      if (xi < pAbsorb) {
        // Absoprtion Event
        return L;
      } else if (xi < 1.0f - pScatter) { 
        // Null Event 
        UpdateWeights(P_OVER_F, null, null / majorants);
        
        wi = ray.d;
      } else {
        // Scatter event 
        realPathLength++;
        numScatterEvents++;

        if (realPathLength >= maxDepth && maxDepth > 0) return L;
        if (enableSingleScatter && numScatterEvents >= 2) return L;
       
        UpdateWeights(P_OVER_F, scatter, scatter / majorants);

        // ========== SAMPLE SCATTER DIRECTION =============//
        Vector3f wo = -ray.d;
        Spectrum pdf_phase = mi.phase->Sample_p(wo, &wi, sampler.Get2D());
        Float pdf_ems = LightPdf(scene, mi, &wi);
        
        Spectrum f_phase = mi.phase->p(wo,wi);
        UpdateWeights(P_OVER_F, f_phase, pdf_phase);
      }
      ray = mi.SpawnRay(wi); 
    } else {
      ++surfaceInteractions; 

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
        ray = isect.SpawnRay(ray.d);
        continue;
      } 

      ++realPathLength;
      lastScatterWasSurface = true;
      if (realPathLength >= maxDepth && maxDepth > 0) return L;

      // Sample BSDF to get new path direction
      Vector3f wo = -ray.d, wi;
      Float bsdfPdf;
      BxDFType flags;
      Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &bsdfPdf,
                                        BSDF_ALL, &flags);

      // If the path contribution function is 0 or pdf is 0 return Spectrum(0.0f).
      if (f.IsBlack() || bsdfPdf == 0.0f) return L;
      specularBounce = (flags & BSDF_SPECULAR) != 0;
      if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)) {
            Float eta = isect.bsdf->eta;
            // Update the term that tracks radiance scaling for refraction
            // depending on whether the ray is entering or leaving the
            // medium.
            etaScale *=
                (Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
      }
      // If the path contribution function is 0 or pdf is 0 return Spectrum(0.0f).
      if (f.IsBlack() || bsdfPdf == 0.0f) {
        return L;
      }

      Spectrum bsdfFunc = f * AbsDot(wi, isect.shading.n);
      UpdateWeights(P_OVER_F, bsdfFunc, bsdfPdf);
      ray = isect.SpawnRay(wi);
    }
  }
}

UnidirectionalSpectralMISIntegrator *CreateUnidirectionalSpectralMISIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera) {
    int maxDepth = params.FindOneInt("maxdepth", 0);
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
    std::string lightStrategy = params.FindOneString("lightsamplestrategy", "spatial");
    
    bool enableSingleScatter = params.FindOneBool("singleScatter", false);

    return new UnidirectionalSpectralMISIntegrator(
      maxDepth, 
      camera, 
      sampler, 
      pixelBounds,
      rrThreshold, 
      lightStrategy,
      enableSingleScatter
    );
}

}  // namespace pbrt
