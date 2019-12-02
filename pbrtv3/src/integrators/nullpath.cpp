// integrators/nullpath.cpp*
#include "integrators/nullpath.h"
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

// NullPathIntegrator Method Definitions
void NullPathIntegrator::Preprocess(const Scene &scene, Sampler &sampler) {
    lightDistribution =
        CreateLightSampleDistribution(lightSampleStrategy, scene);
}

Spectrum
NullPathIntegrator::CalculateTwoStrategyMISWeight(
    std::vector<double> P_OVER_F_1, 
    std::vector<double> P_OVER_F_2,
    Spectrum L,
    bool isStrat1) const {
 
  if (L.IsBlack()) return Spectrum(0.0f);
  
  Spectrum RGB(0.0f);
  if (m_showTwoStrategyMISWeight) {
    // assume that medium is grey for now.
    // store p_1 / (p_1 + p_2) in red channel.
    double strat1_mis_weight = P_OVER_F_1[0] / (P_OVER_F_1[0] + P_OVER_F_2[0]);
    double strat2_mis_weight = P_OVER_F_2[0] / (P_OVER_F_1[0] + P_OVER_F_2[0]);
    RGB[0] = strat1_mis_weight;
    RGB[2] = strat2_mis_weight; 
  } else if (m_splitTwoStrategyMISValue) {
    for (int i = 0; i < 3; i++) {
      double sum = 0.0f;
      for (int j = 0; j < 3; j++) {
        sum += P_OVER_F_1[i * 3 + j] * m_p_channel[j]; 
        sum += P_OVER_F_2[i * 3 + j] * m_p_channel[j];
      }
      if (!isStrat1) {
        RGB[0] = sum == 0.0f ? 0.0f : L[i] / sum;
      } else {
        RGB[2] = sum == 0.0f ? 0.0f : L[i] / sum;
      }
    }
  } else {
    if (isStrat1 && AllZerosVec(P_OVER_F_1)) return Spectrum(0.0f);
    if (!isStrat1 && AllZerosVec(P_OVER_F_2)) return Spectrum(0.0f);

    for (int i = 0; i < 3; i++) {
      double sum = 0.0f;
      for (int j = 0; j < 3; j++) { 
        sum += P_OVER_F_1[i * 3 + j] * m_p_channel[j]; 
        sum += P_OVER_F_2[i * 3 + j] * m_p_channel[j];
      }
      RGB[i] = sum == 0.0f ? 0.0f : L[i] / sum;
    }
  }
  return RGB;
};

Spectrum
NullPathIntegrator::CalculateMISWeight(std::vector<double> P_OVER_F, 
                                       Spectrum L) const {
  if (L.IsBlack()) return Spectrum(0.0f);
  Spectrum RGB(0.0f);
  for (int i = 0; i < 3; i++) {
    double sum = 0.0f;
    for (int j = 0; j < 3; j++) {
      sum += P_OVER_F[i * 3 + j] * m_p_channel[j]; 
    }
    RGB[i] = sum == 0.0f ? 0.0f : L[i] / sum;
  }
  return RGB;
};

int 
NullPathIntegrator::PickHeroWavelength(Sampler &sampler) const {
  int hero = 0;
  Float cdf = m_p_channel[0];
  Float rand = sampler.Get1D();
  while (rand > cdf) {
    hero += 1;
    cdf += m_p_channel[hero];
  }
  return hero;
}

void
NullPathIntegrator::PrintPath(std::stringstream *path) const {
  if (m_printPath) std::cout << "\nPATH: " << path->str() << std::endl;
}

Spectrum NullPathIntegrator::Li(const RayDifferential &r, const Scene &scene,
                               Sampler &sampler, MemoryArena &arena,
                               int depth) const {
    ProfilePhase p(Prof::SamplerIntegratorLi);
    Spectrum L(0.0f);

   // Symbolic Path:
    std::stringstream symbPath;
    symbPath << "\\frac{f(\\overline{x})}{p(\\overline{x})}=\n";

    // P_r / F_r   P_g / F_r   P_b / F_r   
    // P_r / F_g   P_g / F_g   P_b / F_g   
    // P_r / F_b   P_g / F_b   P_b / F_b   
    
    // This path has not forced null events.
    std::vector<double> P_OVER_F = std::vector<double>(9, 1.0f);

    // This path forces null events on the last segment and will
    // sample using light sampling on the last real vertex.
    std::vector<double> P_OVER_F_NEE = std::vector<double>(9, 1.0f);

    // This path forces null events on the last segment. Also, the last
    // vertex will be a scattering event.
    std::vector<double> P_OVER_F_NULL = std::vector<double>(9, 1.0f);

    RayDifferential ray(r);
    bool specularBounce = false;
		bool lastScatterWasSurface = false;
    int bounces = 0;
    int hero = PickHeroWavelength(sampler);
    int numScatterEvents = 0;
    
    int realPathLength = 0;
		

    // Keep track of the surface distance at the last real scatter.
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
          majorants = ray.medium->GetMajorant();
          if (m_enableEqualMajorant) {
            int sampleChannel;
            majorants = HomogenizeMajorants(majorants, sampleChannel);
            Tr = SampleHomogenizedMajorants(ray, sampler, 
                                           arena, &mi, 
                                           sampleChannel);
          } else {         
            Tr = ray.medium->SampleChannel(ray, sampler, arena, &mi, hero);
          }
          Spectrum freeFlightPdf = mi.IsValid() ? majorants * Tr : Tr; 
          UpdateWeights(P_OVER_F, Tr, freeFlightPdf);
          UpdateWeights(P_OVER_F_NEE, Tr, freeFlightPdf);

          if (m_printPath) {
            symbPath << " [sampled distance]";
            symbPath << "\\frac{" <<  Tr[0] << "}";
            if (mi.IsValid()) {
              symbPath << "{" << majorants[0] << "\\cdot" << Tr[0] << "}\n";
            } else {
              symbPath << "{" << freeFlightPdf[0] << "}\n";
            }
          }
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
            if (m_constantEmission > 0.0f) {
							// Emissive media with these methods not currently supported.	
              bool unsupported = m_enableSplitNEE || m_enableOnlyNEE
                                 || m_enableOnlyEA || m_enableEA;
							assert(!unsupported);

              UpdateWeights(P_OVER_F, absorption, absorption / majorants);
              UpdateWeights(P_OVER_F_NEE, absorption, absorption / majorants);

							L += FlagFilter(CalculateMISWeight(P_OVER_F, m_constantEmission), 
															m_ignoreNoScatter,
															m_viewNumScatterEvents,
															realPathLength); 
						}
            
            if (m_printPath) {
              symbPath << " [absorb] ";
              symbPath << "\\frac{" << absorption[0] << "}"; 
              symbPath << "{\\frac{" << absorption[0] << "}{" << majorants[0] << "}}\n";
              symbPath << " [light -- emissive medium] ";
              symbPath << L[0] << "\n";
              PrintPath(&symbPath);
            }

						return L;
          } else if (xi < 1.0f - pScatter) { 
            lastEventWasReal = false;

            // Handle null event.
            UpdateWeights(P_OVER_F, null, null / majorants);

            if (m_enableDeltaTrack) {
              UpdateWeights(P_OVER_F_NEE, null, null / majorants);
            } else {
              UpdateWeights(P_OVER_F_NEE, null, Spectrum(1.0f));
            }

            // Special path weight for EA MIS.
            UpdateWeights(P_OVER_F_NULL, null, majorants);

            wi = ray.d;

            if (m_printPath) {
              symbPath << " [null] ";
              symbPath << "\\frac{" << null[0] << "}"; 
              symbPath << "{\\frac{" << null[0] << "}{" << majorants[0] << "}}\n";
            }
          } else {
            lastEventWasReal = true;
						lastScatterWasSurface = false;
            numScatterEvents++;
            realPathLength++;

            bool exit = ExitEarly(m_singleScatter, 
                                  m_viewNumScatterEvents,
                                  realPathLength);
            bool reachedMax = realPathLength >= maxDepth && maxDepth > 0;
            if (exit || reachedMax) return L;
           
            // Update path to account for scattering. 
            //===================================================//
            UpdateWeights(P_OVER_F, scatter, scatter / majorants);
            
            // Recursive rays": scatter ---> light.
            //===================================================//
            if (m_enableOnlyNEE || m_enableSplitNEE) {
              std::vector<double>  P_OVER_F_END, P_OVER_F_NEE_END;
              Spectrum L_DIRECT = SampleNextEvent(scene, mi, sampler, 
                                                  ray, hero, arena,
                                                  P_OVER_F_END, 
                                                  P_OVER_F_NEE_END,
                                                  m_enableEqualMajorant,
                                                  m_enableDeltaTrack);
              if (m_enableOnlyNEE) {
                L += FlagFilter(CalculateMISWeight(MultiplyVec(P_OVER_F_NEE_END, P_OVER_F), 
                                                    L_DIRECT), 
                                m_ignoreNoScatter,
																m_viewNumScatterEvents,
                                realPathLength);

                // If only NEE and won't count next scatter, exit early.
                bool exitEarly = ExitEarly(m_singleScatter, 
                                    m_viewNumScatterEvents,
                                    realPathLength + 1);
                if (exitEarly) return L;

              } else {
                L += FlagFilter(CalculateTwoStrategyMISWeight(
                                        MultiplyVec(P_OVER_F_END, P_OVER_F), 
                                        MultiplyVec(P_OVER_F_NEE_END, P_OVER_F),
                                        L_DIRECT,
                                        false),
                                 m_ignoreNoScatter,
																 m_viewNumScatterEvents,
                                 realPathLength);
              }
            } else if (m_enableEA) {
              std::vector<double>  P_OVER_F_NEE_END, P_OVER_F_EA_END;
              Spectrum L_DIRECT = SampleNEE(scene, sampler, mi, 
                                          ray, 
                                          lastRealRay,
                                          hero, arena, 
                                          P_OVER_F_NEE_END, P_OVER_F_EA_END,
                                          scatter,
                                          m_enableEqualMajorant,
                                          m_printEAPath); 
              L += FlagFilter(CalculateTwoStrategyMISWeight(
                                        MultiplyVec(P_OVER_F_NEE_END, P_OVER_F),
                                        MultiplyVec(P_OVER_F_EA_END, P_OVER_F_NULL),
                                        L_DIRECT,
                                        true),
                                 m_ignoreNoScatter,
																 m_viewNumScatterEvents,
                                 realPathLength);
            }
            // Sample a direction.
            //===================================================//
            Spectrum scatterFunc, scatterPdf;
            Vector3f wo = -ray.d;
            bool sampledLight;
            if (m_enableOnlyNEE || m_enableSplitNEE) {
              scatterPdf = mi.phase->Sample_p(wo, &wi, sampler.Get2D());
             
              // Update the P_OVER_F_NEE path to account for sampled dir.
              P_OVER_F_NEE = SetWeights(P_OVER_F_NEE, P_OVER_F);
              Float pdf_ems = LightPdf(scene, mi, &wi);
              UpdateWeights(P_OVER_F_NEE, mi.phase->p(wo,wi), pdf_ems);
            } else if (m_enableOnlyEA || m_enableEA) {  
              scatterPdf = mi.phase->Sample_p(wo, &wi, sampler.Get2D());
              P_OVER_F_NULL = SetWeights(P_OVER_F_NULL, P_OVER_F);
            } else if (m_enableNEE) {
              scatterPdf = MISPhaseLightDirection(scene, mi, sampler, wo, &wi, sampledLight);
            } else if (m_enableLightSample) {
              scatterPdf = SampleLightDirection(scene, mi, sampler, &wi);
            } else {
              scatterPdf = mi.phase->Sample_p(wo, &wi, sampler.Get2D());
            }

            scatterFunc = mi.phase->p(wo,wi);
            UpdateWeights(P_OVER_F, scatterFunc, scatterPdf);
             
            scattered = true;

            // Shadow connection. (scatter) ----> scatter ----> light
            //===================================================//
            if (m_enableOnlyEA || m_enableEA) {
              bool skipEA = ExitEarly(m_singleScatter, 
                                    m_viewNumScatterEvents,
                                    realPathLength + 1);

              // If skipping EA will skip on all future scatters.
              // enableOnlyEA will only use EA
              // enableEA also uses onlyNEE, so can only scatter then exit.
              if (skipEA) return L;

              // Update the P_OVER_F_NULL path to account.
              P_OVER_F_NULL = SetWeights(P_OVER_F_NULL, P_OVER_F);

              // Calculate an EA connection to a light.
              std::vector<double> P_OVER_F_NEE_END, P_OVER_F_EA_END;
              Spectrum L_FROM_PATH = SampleEA(scene, sampler, mi.SpawnRay(wi), hero, arena,
                                              P_OVER_F_NEE_END,
                                              P_OVER_F_EA_END,
                                              m_enableEqualMajorant,
                                              m_printEAPath);
                    
              // Connect path contributions to old path.
              if (m_enableOnlyEA) {
                L += FlagFilter(CalculateMISWeight(MultiplyVec(P_OVER_F_EA_END, 
                                                               P_OVER_F), 
                                                    L_FROM_PATH), 
                                m_ignoreNoScatter,
																m_viewNumScatterEvents,
                                realPathLength + 1);

              } else if (m_enableEA) {
                L += FlagFilter(CalculateTwoStrategyMISWeight(
                                        MultiplyVec(P_OVER_F_EA_END, P_OVER_F), 
                                        MultiplyVec(P_OVER_F_NEE_END, P_OVER_F),
                                        L_FROM_PATH,
                                        true),
                                 m_ignoreNoScatter,
																 m_viewNumScatterEvents, 
                                 realPathLength + 1);
              }
          }
          if (m_printPath) { 
            symbPath << " [scatter] ";
            symbPath << "\\frac{" << scatter[0] << "}"; 
            symbPath << "{\\frac{" << scatter[0] << "}{" << majorants[0] << "}}";
            symbPath << "\\frac{" << scatterFunc[0] << "}{" << scatterPdf[0] << "}\n";
            if (m_enableNEE) {
              if (sampledLight) symbPath << "sampled light-->";
              else symbPath << "sampled phase-->";
              symbPath << "phase_pdf: " << scatterFunc[0] << "  ems_pdf: " << scatterPdf[0] - (0.5 * scatterFunc[0]) << "\n";
            }
          }
          if (m_showDirPdf) {
              Spectrum s;
              s[0] = scatterFunc[0] / scatterPdf[0]; // ratio
              s[1] = scatterFunc[0]; // scatter pdf
              s[2] = scatterPdf[0] - (0.5 * scatterFunc[0]); // light pdf
              return s;
          }
        }
        ray = mi.SpawnRay(wi); 
      } else {
        ++surfaceInteractions; 
        ++realPathLength;
				
        // using nee-unidirectional MIS
				//	- scattered in medium before hit
				
				// We count direct hits (i.e. without MIS weighting) when:
				//	using nee-unidirectional MIS
				//		- ray is from camera or specular surface
				//		- last bounce was surface and no surface NEE
				//	using nee, only ea, or mis ea
				//		- ray is from camera or specular surface
				//		- last bounce was surface and no surface NEE
				//	using naive method
				//		- ray is from camera or specular surface
				//		- last bounce was surface and no surface NEE
				//    - scattered in medium before hit
				//
				//
				bool nextEventMethod = m_enableOnlyNEE || m_enableOnlyEA || m_enableEA;

				bool rayFromCamera = realPathLength == 0 && !lastScatterWasSurface;
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

        // If intersection goes out of scene, look up the environment lights.
        if (!foundIntersection) {
          for(auto light : scene.infiniteLights) {
            if (countMIS) {
              L += FlagFilter(CalculateTwoStrategyMISWeight(
                                        P_OVER_F, 
                                        P_OVER_F_NEE, 
                                        light->Le(ray),
                                        true), 
                              m_ignoreNoScatter,
                              m_viewNumScatterEvents,
                              realPathLength);
            } else if (countDirect) {
              L += FlagFilter(CalculateMISWeight(P_OVER_F, light->Le(ray)), 
                              m_ignoreNoScatter,
                              m_viewNumScatterEvents,
                              realPathLength);
            }
          }
          if (m_printPath) { 
            symbPath << " [light -- out of scene] ";
            symbPath << L[0];
            if (!m_printHits) {
              PrintPath(&symbPath);
            }
          }
          return L;
        }
        
        if (!isect.Le(-ray.d).IsBlack()) {
          if (countMIS) {
            L += FlagFilter(CalculateTwoStrategyMISWeight(
                                      P_OVER_F, 
                                      P_OVER_F_NEE, 
                                      isect.Le(-ray.d),
                                      true),
                            m_ignoreNoScatter,
                            m_viewNumScatterEvents,
                            realPathLength);
          } else if (countDirect) {
            L += FlagFilter(CalculateMISWeight(P_OVER_F, isect.Le(-ray.d)), 
                            m_ignoreNoScatter,
                            m_viewNumScatterEvents,
                            realPathLength);
          }

          if (m_printPath) { 
            symbPath << " [light -- hit surface(";
            symbPath << "p/f=" << P_OVER_F[0] << ", ";
            symbPath << "light=" << isect.Le(-ray.d) << ")]=";
            symbPath << L[0];
            PrintPath(&symbPath);
          }

          return L;
        }
        
        // Compute scattering functions.
        isect.ComputeScatteringFunctions(ray, arena, true);

        // If no bsdf at intersection, medium boundary so continue. 
        if (!isect.bsdf) {
          --realPathLength;

          if (m_printPath)  symbPath << " [Medium Boundary] \n"; 
 
          // Update the P_OVER_F_NULL path to account.
          P_OVER_F_NULL = SetWeights(P_OVER_F_NULL, P_OVER_F);
        
          ray = isect.SpawnRay(ray.d);

          // Add one scattering event since this will be included in 
          // the path.
          bool skipEA = ExitEarly(m_singleScatter, 
                                  m_viewNumScatterEvents,
                                  realPathLength+ 1);

          if (skipEA || !(m_enableOnlyEA || m_enableEA)) continue;

          // Calculate an EA connection to a light.
          std::vector<double> P_OVER_F_NEE_END, P_OVER_F_EA_END;
          // Will be used for EA shadow connection.
          Spectrum L_FROM_PATH = SampleEA(scene, sampler,
                                          ray, hero, arena,
                                          P_OVER_F_NEE_END,
                                          P_OVER_F_EA_END,
                                          m_enableEqualMajorant,
                                          m_printEAPath);
          
          // Connect path contributions to old path.
          if (m_enableOnlyEA) {
            L += FlagFilter(CalculateMISWeight(MultiplyVec(P_OVER_F_EA_END, 
                                                           P_OVER_F), 
                                                L_FROM_PATH), 
                            m_ignoreNoScatter,
                            m_viewNumScatterEvents,
                            realPathLength + 1);
          } else if (m_enableEA) {
            L += FlagFilter(CalculateTwoStrategyMISWeight(
                                    MultiplyVec(P_OVER_F_EA_END, P_OVER_F), 
                                    MultiplyVec(P_OVER_F_NEE_END, P_OVER_F),
                                    L_FROM_PATH,
                                    true),
                             m_ignoreNoScatter,
                             m_viewNumScatterEvents,
                             realPathLength + 1);
          }
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

        if (!specularBounce && m_enableSurfaceNEE) {
          std::vector<double>  P_OVER_F_END, P_OVER_F_NEE_END;
          Spectrum L_DIRECT = SampleNextEvent(scene, isect, sampler, 
                          ray, hero, arena, 
                          P_OVER_F_END, 
                          P_OVER_F_NEE_END, 
                          m_enableEqualMajorant, 
                          m_enableDeltaTrack);
					L += FlagFilter(
								CalculateMISWeight(
									MultiplyVec(
										P_OVER_F_NEE_END,
										P_OVER_F
									),
									L_DIRECT),
								m_ignoreNoScatter,
								m_viewNumScatterEvents,
							  realPathLength	
							);
				}
        
        Spectrum bsdfFunc = f * AbsDot(wi, isect.shading.n);
        UpdateWeights(P_OVER_F, bsdfFunc, bsdfPdf);

			  ray = isect.SpawnRay(wi);

        if (m_printPath) {
          symbPath << " [surface] ";
          symbPath << "\\frac{" << bsdfFunc[0] << "}"; 
          symbPath << "{" << bsdfPdf << "}\n";
        } 
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

    Spectrum p_channel = params.FindOneSpectrum(
                                  "p_channel", 
                                  Spectrum(1.0f / Spectrum::nSamples));

    // Stop tracing after first scattering event.
    bool singleScatter = params.FindOneBool("singleScatter", false);
   
    // Does not show paths that don't include at least one scattering event.
    bool ignoreNoScatter = params.FindOneBool("ignoreNoScatter", false);
	
    // Print out the path for each sample.
    bool printPath = params.FindOneBool("printPath", false);
    bool printHits = params.FindOneBool("printHits", false);

    // Print out the path for each sample.
    bool printEAPath = params.FindOneBool("printEAPath", false);

    // Only view paths with the specified number of scattering events.
    int viewNumScatterEvents = params.FindOneInt("viewNumScatterEvents", -1);

    // Enable constant emission in the medium..
    Float constantEmission = params.FindOneFloat("constantEmission", 0);

    // Create a visual of the intensity of the MIS weights.
    bool showTwoStrategyMISWeight = params.FindOneBool(
				"showTwoStrategyMISWeight", 
				false
		);
   
    // Split the estimator values for each MIS strategy.
    bool splitTwoStrategyMISValue = params.FindOneBool(
				"splitTwoStrategyMISValue", 
				false
		);

    bool showDirPdf = params.FindOneBool("showDirPdf", 0);

		// Use surface NEE
		bool enableSurfaceNEE = params.FindOneBool("enable-surface-NEE", false);

    // Use only EA sampling for distances.
    bool enableOnlyEA = params.FindOneBool("enable-only-EA", false);

    // MIS between EA and non EA sampling
    bool enableEA = params.FindOneBool("enable-EA", false);
		
    // Use directional NEE.
    bool enableNEE = params.FindOneBool("enable-NEE", false);

    // Use proper NEE.
    bool enableSplitNEE = params.FindOneBool("enable-split-NEE", false);

    // Use only shadow connections.
    bool enableOnlyNEE = params.FindOneBool("enable-only-NEE", false);

    // turn on delta tracking for NEE
    bool enableDeltaTrack = params.FindOneBool("enable-delta-track", false);

    // Choose direction based on light sampling (biased if not single scatter).
    bool enableLightSample = params.FindOneBool("enable-LightSample", false);

    // Homogenize the medium.
    bool enableEqualMajorant = params.FindOneBool("enable-equal-majorant", false);

    // Use the balance heuristic.
    bool enableBalanceHeuristic = params.FindOneBool("enable-balance-heuristic", true);

    return new NullPathIntegrator(maxDepth, 
                                camera, 
                                sampler, 
                                pixelBounds,
                                p_channel, 
                                rrThreshold, 
                                lightStrategy, 
                                enableShadows,
																enableSurfaceNEE,
                                enableEA,
                                enableOnlyEA,
                                enableNEE, 
                                enableSplitNEE, 
                                enableOnlyNEE,
                                enableLightSample, 
                                enableEqualMajorant,
                                enableDeltaTrack,
                                singleScatter,
                                viewNumScatterEvents,
                                constantEmission,
                                ignoreNoScatter,
                                enableBalanceHeuristic,
                                printPath,
                                printHits,
                                printEAPath,
                                showTwoStrategyMISWeight,
                                splitTwoStrategyMISValue,
                                showDirPdf);
}

}  // namespace pbrt
