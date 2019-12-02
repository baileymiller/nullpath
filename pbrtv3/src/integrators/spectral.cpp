// integrators/spectral.cpp*
#include "integrators/spectral.h"
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
void SpectralIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
    lightDistribution =
        CreateLightSampleDistribution(lightSampleStrategy, scene);
}

Spectrum SpectralIntegrator::CalculateMISWeight(Spectrum L, 
																							  Spectrum W, Spectrum W_NEE) const {
	Spectrum res;
	for (int i = 0; i < 3; i++) {
    Float P_OVER_F = W[i] = 0.0f ? 0.0f : 1.0f / W[i];
    Float P_OVER_F_NEE = W_NEE[i] = 0.0f ? 0.0f : 1.0f / W_NEE[i];
		res[i] = L[i] / (P_OVER_F + P_OVER_F_NEE);
	}
	return res;
}

Spectrum SpectralIntegrator::Li(const RayDifferential &r, const Scene &scene,
                               Sampler &sampler, MemoryArena &arena,
                               int depth) const {
    ProfilePhase p(Prof::SamplerIntegratorLi);
    Spectrum W(1.0f), W_NEE(1.0f), L(0.0f);
    RayDifferential ray(r);
    bool specularBounce = false;
		bool lastScatterWasSurface = false;
    int numScatterEvents = 0;
    int bounces = 0; 
    Float etaScale = 1;

    int channel = -1;
    if (m_enableNaiveDelta) {
      // Choose a channel randomly
      channel = static_cast<int>(sampler.Get1D() * 3);
      for (int i = 0; i < 3; i++) {
        if (channel ==i) {
          W[i] = 3.0f;
          W_NEE[i] = 3.0f;
        } else {
          W[i] = 0.0f;
          W_NEE[i] = 0.0f;
        }
      }
    }

    for(bounces = 0;; bounces++) {
      // Intersect _ray_ with scene and store intersection in _isect_
      SurfaceInteraction isect;
      bool foundIntersection = scene.Intersect(ray, &isect);

      // Sample the participating medium, if present
      MediumInteraction mi;
      Float M = 0.0f;
      if (ray.medium) {
        Spectrum majorants = ray.medium->GetMajorant();
        int sampleChannel;
        for (int i=0; i<Spectrum::nSamples; i++){
          if (M < majorants[i]) {
            M = majorants[i];
            sampleChannel = i;
          }
        }
        
        if (channel >= 0) {
          M = majorants[channel];
          sampleChannel = channel;
        }

        ray.medium->SampleChannel(ray, sampler, arena, &mi, sampleChannel);
      }
      
      // Handle an interaction with a medium or a surface
      if (mi.IsValid()) {
        ++volumeInteractions;

        Spectrum scatter, absorption, _;
        ray.medium->GetCoefficients(mi.p, scatter, absorption, _);

        Float pScatter = 0.0f, pAbsorb = 0.0f, pNull = 0.0f;
        for (int i = 0; i < Spectrum::nSamples; i ++) {
          pScatter += (W[i] * scatter[i]);
          pAbsorb += (W[i] * absorption[i]); 
          pNull += (W[i] * (M - (scatter[i] + absorption[i])));
        }

        pScatter /= Spectrum::nSamples;
        pAbsorb /= Spectrum::nSamples;
        pNull /=  Spectrum::nSamples;
        Float pNorm = pScatter + pAbsorb + pNull;
        pScatter /= pNorm;
        pAbsorb /= pNorm;
        pNull /= pNorm;
       
        Vector3f wi;
        float xi = sampler.Get1D();
        if (xi < pAbsorb) {
          // Handle absorption.
          return L; 
        } else if (xi < 1.0f - pScatter) {
          // Handle null event.
          W *= (Spectrum(M) - scatter - absorption)  / (M * pNull);
					W_NEE *= (Spectrum(M) - scatter - absorption)  / (M * pNull);
          wi = ray.d;
        } else {
          numScatterEvents++;
					lastScatterWasSurface = false;

          bool exit = ExitEarly(m_singleScatter,
                                m_viewNumScatterEvents,
                                numScatterEvents);
          if (exit) return L;

          // Update path to account for scattering
          // ===========================================//
          W *= scatter / (M * pScatter);

          // Recursive rays: scatter ---> light.
          // ===========================================//
          if (m_enableOnlyNEE || m_enableSplitNEE) { 	
						Spectrum nee_W = Spectrum(W);
						Spectrum nee_W_NEE = Spectrum(W);
						
            Spectrum L_NEE = SampleSpectralNextEvent(scene, mi, 
                                                     sampler, ray, arena, 
																				             nee_W, nee_W_NEE, 
                                                     channel);

            if (m_enableOnlyNEE) {
              L += FlagFilter(L_NEE * nee_W_NEE, 
                              m_ignoreNoScatter, 
                              m_viewNumScatterEvents, 
                              numScatterEvents); 
              // If only NEE and won't count next scatter, exit early.
              bool exitEarly = ExitEarly(m_singleScatter, m_viewNumScatterEvents, numScatterEvents + 1);
              if (exitEarly) return L;
            } else {
              W_NEE = W;
						  L += FlagFilter(CalculateMISWeight(L_NEE, nee_W, nee_W_NEE),
														m_ignoreNoScatter,
														m_viewNumScatterEvents,
														numScatterEvents);
            }
          }
         
          // Sample a direction.
          //===================================================//
          Vector3f wo = -ray.d;
          bool sampledLight;
          Spectrum scatterFunc, scatterPdf;
          if (m_enableOnlyNEE || m_enableSplitNEE) {
            scatterPdf = mi.phase->Sample_p(wo, &wi, sampler.Get2D());

            // Update W_NEE path to account for sampled dir.
						Float pdf_ems = LightPdf(scene, mi, &wi);
            if (pdf_ems == 0.0f) W_NEE = Spectrum(0.0f);
            else W_NEE *= mi.phase->p(wo,wi) / pdf_ems;
          } else if (m_enableNEE) {
            scatterPdf = MISPhaseLightDirection(scene, mi, sampler, wo, &wi, 
                                                sampledLight);
          } else {
            scatterPdf = mi.phase->Sample_p(wo, &wi, sampler.Get2D());
          } 
          scatterFunc = mi.phase->p(wo,wi);
          W *= scatterFunc / scatterPdf;
 
        }
        ray = mi.SpawnRay(wi);
      } else {
        ++surfaceInteractions;
        
				bool nextEventMethod = m_enableOnlyNEE;

				bool rayFromCamera = numScatterEvents == 0 && !lastScatterWasSurface;
				bool lastScatterWasMedium = numScatterEvents > 0 && !lastScatterWasSurface;
				bool lastSurfaceUsedNEE = m_enableSurfaceNEE && !specularBounce;

				bool countMIS = false;
			  bool countDirect = false;

				if (lastScatterWasMedium) {
					// Hit a particle in the medium before reaching this point.
					if (m_enableSplitNEE) countMIS = true;	// MIS direct hits + nee
					else if (!nextEventMethod) countDirect = true; // Only doing direct hits
				} else if (rayFromCamera || !lastSurfaceUsedNEE) {
					// Reached this point from camera  or surface with no scattering events.
					countDirect = true; // Always count these direct hits.
				}

        if (!foundIntersection) {
          for(auto &light : scene.infiniteLights) {
            if (countMIS) {
              L += FlagFilter(CalculateMISWeight(light->Le(ray), W, W_NEE),
                              m_ignoreNoScatter,
                              m_viewNumScatterEvents,
                              numScatterEvents);
            } else if (countDirect) {
              L += FlagFilter(light->Le(ray) * W,
                              m_ignoreNoScatter,
                              m_viewNumScatterEvents,
                              numScatterEvents);
            }
          }
          return L;
        }

        if (!isect.Le(-ray.d).IsBlack())  {
          if (countMIS) {
            L += FlagFilter(CalculateMISWeight(isect.Le(-ray.d), W, W_NEE),
                            m_ignoreNoScatter,
                            m_viewNumScatterEvents,
                            numScatterEvents);
          } else if (countDirect){
            L += FlagFilter(isect.Le(-ray.d) * W,
                              m_ignoreNoScatter,
                              m_viewNumScatterEvents,
                              numScatterEvents);
          }
          return L;
        } 
          
        // Compute scattering functions and skip over medium boundaries
        isect.ComputeScatteringFunctions(ray, arena, true);
        if (!isect.bsdf) {
          ray = isect.SpawnRay(ray.d);
          bounces--;
          continue;
        }
				
				lastScatterWasSurface = true;

        // Sample BSDF to get new path direction
        Vector3f wo = -ray.d, wi;
        Float bsdfPdf;
        BxDFType flags;
        Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &bsdfPdf,
                                          BSDF_ALL, &flags);

        if (f.IsBlack() || bsdfPdf == 0.f) return L;
        specularBounce = (flags & BSDF_SPECULAR) != 0;

        if ((flags & BSDF_SPECULAR) && (flags & BSDF_TRANSMISSION)) {
					Float eta = isect.bsdf->eta;
					// Update the term that tracks radiance scaling for refraction
					// depending on whether the ray is entering or leaving the
					// medium.
					etaScale *=
							(Dot(wo, isect.n) > 0) ? (eta * eta) : 1 / (eta * eta);
        }

				if (!specularBounce && m_enableSurfaceNEE) {
				  Spectrum nee_W = Spectrum(W);
					Spectrum nee_W_NEE = Spectrum(W);
          Spectrum L_NEE = SampleSpectralNextEvent(scene, isect, 
                                                  sampler, ray, arena, 
                                                  nee_W, nee_W_NEE, 
                                                  channel);
          L += FlagFilter(L_NEE * nee_W_NEE, 
													m_ignoreNoScatter, 
													m_viewNumScatterEvents, 
													numScatterEvents); 
        }

        ray = isect.SpawnRay(wi);
        W *= f * AbsDot(wi, isect.shading.n) / bsdfPdf;
      }
    }
}

SpectralIntegrator *CreateSpectralIntegrator(
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
    bool enableSurfaceNEE = params.FindOneBool("enable-surface-NEE", false);
    bool enableNaiveDelta = params.FindOneBool("enable-naive-delta", false);
    bool enableNEE = params.FindOneBool("enable-NEE", false);
		bool enableSplitNEE = params.FindOneBool("enable-split-NEE", false);
    bool enableOnlyNEE = params.FindOneBool("enable-only-NEE", false);
    bool singleScatter = params.FindOneBool("singleScatter", false);
    int viewNumScatterEvents = params.FindOneInt("viewNumScatterEvents", -1);
    bool ignoreNoScatter = params.FindOneBool("ignoreNoScatter", false);

    return new SpectralIntegrator(maxDepth, camera, sampler, pixelBounds,
                                 rrThreshold, lightStrategy, enableShadows,
                                 enableSurfaceNEE, enableNaiveDelta, enableNEE, 
																 enableSplitNEE, enableOnlyNEE, 
                                 singleScatter, viewNumScatterEvents, ignoreNoScatter);
}

}  // namespace pbrt
