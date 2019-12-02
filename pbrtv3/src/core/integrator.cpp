
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// core/integrator.cpp*
#include "integrator.h"
#include "scene.h"
#include "interaction.h"
#include "sampling.h"
#include "parallel.h"
#include "film.h"
#include "sampler.h"
#include "integrator.h"
#include "progressreporter.h"
#include "camera.h"
#include "stats.h"

// Includes for EA sampling //
#include "lights/diffuse.h"
// END Includes for EA sampling //

//#define DEBUG_PIXEL

#ifdef DEBUG_PIXEL
const int pixelX = 117;
const int pixelY = 135;
#endif

namespace pbrt {

STAT_COUNTER("Integrator/Camera rays traced", nCameraRays);

// Integrator Method Definitions
Integrator::~Integrator() {}

// Integrator Utility Functions
Spectrum UniformSampleAllLights(const Interaction &it, const Scene &scene,
                                MemoryArena &arena, Sampler &sampler,
                                const std::vector<int> &nLightSamples,
                                bool handleMedia, bool enableShadows) {
    ProfilePhase p(Prof::DirectLighting);
    Spectrum L(0.f);
    for (size_t j = 0; j < scene.lights.size(); ++j) {
        // Accumulate contribution of _j_th light to _L_
        const std::shared_ptr<Light> &light = scene.lights[j];
        int nSamples = nLightSamples[j];
        const Point2f *uLightArray = sampler.Get2DArray(nSamples);
        const Point2f *uScatteringArray = sampler.Get2DArray(nSamples);
        if (!uLightArray || !uScatteringArray) {
            // Use a single sample for illumination from _light_
            Point2f uLight = sampler.Get2D();
            Point2f uScattering = sampler.Get2D();
            L += EstimateDirect(it, uScattering, *light, uLight, scene, sampler,
                                arena, handleMedia);
        } else {
            // Estimate direct lighting using sample arrays
            Spectrum Ld(0.f);
            for (int k = 0; k < nSamples; ++k)
                Ld += EstimateDirect(it, uScatteringArray[k], *light,
                                     uLightArray[k], scene, sampler, arena,
                                     handleMedia, enableShadows);
            L += Ld / nSamples;
        }
    }
    return L;
}

Spectrum UniformSampleOneLight(const Interaction &it, const Scene &scene,
                               MemoryArena &arena, Sampler &sampler,
                               bool handleMedia, const Distribution1D *lightDistrib, bool enableShadows) {
    ProfilePhase p(Prof::DirectLighting);
    // Randomly choose a single light to sample, _light_
    int nLights = int(scene.lights.size());
    if (nLights == 0) return Spectrum(0.f);
    int lightNum;
    Float lightPdf;
    if (lightDistrib) {
        lightNum = lightDistrib->SampleDiscrete(sampler.Get1D(), &lightPdf);
        if (lightPdf == 0) return Spectrum(0.f);
    } else {
        lightNum = std::min((int)(sampler.Get1D() * nLights), nLights - 1);
        lightPdf = Float(1) / nLights;
    }
    const std::shared_ptr<Light> &light = scene.lights[lightNum];
    Point2f uLight = sampler.Get2D();
    Point2f uScattering = sampler.Get2D();
    return EstimateDirect(it, uScattering, *light, uLight,
                          scene, sampler, arena, handleMedia, enableShadows) / lightPdf;
}

Spectrum EstimateDirect(const Interaction &it, const Point2f &uScattering,
                        const Light &light, const Point2f &uLight,
                        const Scene &scene, Sampler &sampler,
                        MemoryArena &arena, bool handleMedia, bool enableShadows, bool specular) {
    BxDFType bsdfFlags =
        specular ? BSDF_ALL : BxDFType(BSDF_ALL & ~BSDF_SPECULAR);
    Spectrum Ld(0.f);
    // Sample light source with multiple importance sampling
    Vector3f wi;
    Float lightPdf = 0, scatteringPdf = 0;
    VisibilityTester visibility;
    Spectrum Li = light.Sample_Li(it, uLight, &wi, &lightPdf, &visibility);
	uint32_t flags = 0;
    resetVertex(flags);
    VLOG(2) << "EstimateDirect uLight:" << uLight << " -> Li: " << Li << ", wi: "
            << wi << ", pdf: " << lightPdf;
    if (lightPdf > 0 && !Li.IsBlack()) {
        // Compute BSDF or phase function's value for light sample
        Spectrum f;
        if (it.IsSurfaceInteraction()) {
            // Evaluate BSDF for light sampling strategy
            const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
            f = isect.bsdf->f(isect.wo, wi, bsdfFlags) *
                AbsDot(wi, isect.shading.n);
            scatteringPdf = isect.bsdf->Pdf(isect.wo, wi, bsdfFlags);
            VLOG(2) << "  surf f*dot :" << f << ", scatteringPdf: " << scatteringPdf;
			
			// set the first vertex to be surface and lights are always surfaces
			setVertex0Type(CustomVertexType::VERTEX_SURFACE, flags);
			setVertex1Type(CustomVertexType::VERTEX_SURFACE, flags);
        } else {
            // Evaluate phase function for light sampling strategy
            const MediumInteraction &mi = (const MediumInteraction &)it;
            Float p = mi.phase->p(mi.wo, wi);
            f = Spectrum(p);
            scatteringPdf = p;
            VLOG(2) << "  medium p: " << p;

			// set the first vertex to medium and lights are always surfaces
			setVertex0Type(CustomVertexType::VERTEX_MEDIUM, flags);
			setVertex1Type(CustomVertexType::VERTEX_SURFACE, flags);
        }
        if (!f.IsBlack()) {
            // Compute effect of visibility for light source sample
			if (handleMedia) {
                Li *= visibility.Tr(scene, sampler);
                VLOG(2) << "  after Tr, Li: " << Li;
            } else {
              if (!visibility.Unoccluded(scene)) {
                VLOG(2) << "  shadow ray blocked";
                Li = Spectrum(0.f);
              } else
                VLOG(2) << "  shadow ray unoccluded";
            }

            // Add light's contribution to reflected radiance
            if (!Li.IsBlack()) {
                if (IsDeltaLight(light.flags))
                    Ld += f * Li / lightPdf;
                else {
                    Float weight =
                        PowerHeuristic(1, lightPdf, 1, scatteringPdf);
                    Ld += f * Li * weight / lightPdf;
                }
            }
        }
    }

    // Sample BSDF with multiple importance sampling
    if (!IsDeltaLight(light.flags)) {
        resetVertex(flags);
        Spectrum f;
        bool sampledSpecular = false;
        if (it.IsSurfaceInteraction()) {
            // Sample scattered direction for surface interactions
            BxDFType sampledType;
            const SurfaceInteraction &isect = (const SurfaceInteraction &)it;
            f = isect.bsdf->Sample_f(isect.wo, &wi, uScattering, &scatteringPdf,
                                     bsdfFlags, &sampledType);
            f *= AbsDot(wi, isect.shading.n);
            sampledSpecular = (sampledType & BSDF_SPECULAR) != 0;
			setVertex0Type(CustomVertexType::VERTEX_SURFACE, flags);
			setVertex1Type(CustomVertexType::VERTEX_SURFACE, flags);
        } else {
            // Sample scattered direction for medium interactions
            const MediumInteraction &mi = (const MediumInteraction &)it;
            Float p = mi.phase->Sample_p(mi.wo, &wi, uScattering);
            f = Spectrum(p);
            scatteringPdf = p;
			setVertex0Type(CustomVertexType::VERTEX_MEDIUM, flags);
			setVertex1Type(CustomVertexType::VERTEX_SURFACE, flags);
        }
        VLOG(2) << "  BSDF / phase sampling f: " << f << ", scatteringPdf: " <<
            scatteringPdf;
        if (!f.IsBlack() && scatteringPdf > 0) {
            // Account for light contributions along sampled direction _wi_
            Float weight = 1;
            if (!sampledSpecular) {
                lightPdf = light.Pdf_Li(it, wi);
                if (lightPdf == 0) return Ld;
                weight = PowerHeuristic(1, scatteringPdf, 1, lightPdf);
            }

            // Find intersection and compute transmittance
            SurfaceInteraction lightIsect;
            Ray ray = it.SpawnRay(wi);
            Spectrum Tr(1.f);
            bool foundSurfaceInteraction =
                handleMedia ? scene.IntersectTr(ray, sampler, &lightIsect, &Tr)
                            : scene.Intersect(ray, &lightIsect);
			if (!enableShadows)
				Tr = Spectrum(1.0f);	// revert any shadowing or transmittance computation

            // Add light contribution from material sampling
            Spectrum Li(0.f);
            if (foundSurfaceInteraction) {
                if (lightIsect.primitive->GetAreaLight() == &light)
                    Li = lightIsect.Le(-wi);
            } else
                Li = light.Le(ray);
            if (!Li.IsBlack()) Ld += f * Li * Tr * weight / scatteringPdf;
        }
    }
    return Ld;
}

std::unique_ptr<Distribution1D> ComputeLightPowerDistribution(
    const Scene &scene) {
    if (scene.lights.empty()) return nullptr;
    std::vector<Float> lightPower;
    for (const auto &light : scene.lights)
        lightPower.push_back(light->Power().y());
    return std::unique_ptr<Distribution1D>(
        new Distribution1D(&lightPower[0], lightPower.size()));
}

// SamplerIntegrator Method Definitions
void SamplerIntegrator::Render(const Scene &scene) {
    Preprocess(scene, *sampler);
    // Render image tiles in parallel

    // Compute number of tiles, _nTiles_, to use for parallel rendering
    Bounds2i sampleBounds = camera->film->GetSampleBounds();
    Vector2i sampleExtent = sampleBounds.Diagonal();
    const int tileSize = 16;
    Point2i nTiles((sampleExtent.x + tileSize - 1) / tileSize,
                   (sampleExtent.y + tileSize - 1) / tileSize);
    ProgressReporter reporter(nTiles.x * nTiles.y, "Rendering");
    {
        ParallelFor2D([&](Point2i tile) {
            // Render section of image corresponding to _tile_

            // Allocate _MemoryArena_ for tile
            MemoryArena arena;

            // Get sampler instance for tile
            int tileIdx = tile.y * nTiles.x + tile.x;
            int seed;
            
            // The first case is for backwards compatibility
            if (PbrtOptions.seed == 0) {
                seed = tileIdx;
            } else {
                // We abuse a RNG to hackily convert a 64-bit seed into a 32-bit seed.
                RNG rng{PbrtOptions.seed | ((uint64_t)tileIdx) << 32};
                seed = (int)rng.UniformUInt32();
            }

            std::unique_ptr<Sampler> tileSampler = sampler->Clone(seed);

            // Compute sample bounds for tile
            int x0 = sampleBounds.pMin.x + tile.x * tileSize;
            int x1 = std::min(x0 + tileSize, sampleBounds.pMax.x);
            int y0 = sampleBounds.pMin.y + tile.y * tileSize;
            int y1 = std::min(y0 + tileSize, sampleBounds.pMax.y);
            Bounds2i tileBounds(Point2i(x0, y0), Point2i(x1, y1));
            LOG(INFO) << "Starting image tile " << tileBounds;

            // Get _FilmTile_ for tile
            std::unique_ptr<FilmTile> filmTile =
                camera->film->GetFilmTile(tileBounds);

            // Loop over pixels in tile to render them
            for (Point2i pixel : tileBounds) {
                {
                    ProfilePhase pp(Prof::StartPixel);
                    tileSampler->StartPixel(pixel);
                }

                // Do this check after the StartPixel() call; this keeps
                // the usage of RNG values from (most) Samplers that use
                // RNGs consistent, which improves reproducability /
                // debugging.
                if (!InsideExclusive(pixel, pixelBounds))
                    continue;

                do {
                    // Initialize _CameraSample_ for current sample
                    CameraSample cameraSample =
                        tileSampler->GetCameraSample(pixel);

                    // Generate camera ray for current sample
                    RayDifferential ray;
                    Float rayWeight =
                        camera->GenerateRayDifferential(cameraSample, &ray);
                    ray.ScaleDifferentials(
                        1 / std::sqrt((Float)tileSampler->samplesPerPixel));
                    ++nCameraRays;

                    // Evaluate radiance along camera ray
                    Spectrum L(0.f);
#ifdef DEBUG_PIXEL
					if (pixel.x == pixelX && pixel.y == pixelY)
					{
						// DEBUG PIXEL CODE
#endif
					if (rayWeight > 0) L = Li(ray, scene, *tileSampler, arena);
#ifdef DEBUG_PIXEL
						// DEBUG PIXEL CODE
					}
#endif

                    // Issue warning if unexpected radiance value returned
                    if (L.HasNaNs()) {
                        LOG(ERROR) << StringPrintf(
                            "Not-a-number radiance value returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Spectrum(0.f);
                    } else if (L.y() < -1e-5) {
                        LOG(ERROR) << StringPrintf(
                            "Negative luminance value, %f, returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            L.y(), pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Spectrum(0.f);
                    } else if (std::isinf(L.y())) {
                          LOG(ERROR) << StringPrintf(
                            "Infinite luminance value returned "
                            "for pixel (%d, %d), sample %d. Setting to black.",
                            pixel.x, pixel.y,
                            (int)tileSampler->CurrentSampleNumber());
                        L = Spectrum(0.f);
                    }
                    VLOG(1) << "Camera sample: " << cameraSample << " -> ray: " <<
                        ray << " -> L = " << L;

                    // Add camera ray's contribution to image
                    filmTile->AddSample(cameraSample.pFilm, L, rayWeight);

                    // Free _MemoryArena_ memory from computing image sample
                    // value
                    arena.Reset();
                } while (tileSampler->StartNextSample());
            }
            LOG(INFO) << "Finished image tile " << tileBounds;

            // Merge image tile into _Film_
            camera->film->MergeFilmTile(std::move(filmTile));
            reporter.Update();
        }, nTiles);
        reporter.Done();
    }
    LOG(INFO) << "Rendering finished";

    // Save final image after rendering
    camera->film->WriteImage();
}

Spectrum SamplerIntegrator::SpecularReflect(
    const RayDifferential &ray, const SurfaceInteraction &isect,
    const Scene &scene, Sampler &sampler, MemoryArena &arena, int depth) const {
    // Compute specular reflection direction _wi_ and BSDF value
    Vector3f wo = isect.wo, wi;
    Float pdf;
    BxDFType type = BxDFType(BSDF_REFLECTION | BSDF_SPECULAR);
    Spectrum f = isect.bsdf->Sample_f(wo, &wi, sampler.Get2D(), &pdf, type);

    // Return contribution of specular reflection
    const Normal3f &ns = isect.shading.n;
    if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
        // Compute ray differential _rd_ for specular reflection
        RayDifferential rd = isect.SpawnRay(wi);
        if (ray.hasDifferentials) {
            rd.hasDifferentials = true;
            rd.rxOrigin = isect.p + isect.dpdx;
            rd.ryOrigin = isect.p + isect.dpdy;
            // Compute differential reflected directions
            Normal3f dndx = isect.shading.dndu * isect.dudx +
                            isect.shading.dndv * isect.dvdx;
            Normal3f dndy = isect.shading.dndu * isect.dudy +
                            isect.shading.dndv * isect.dvdy;
            Vector3f dwodx = -ray.rxDirection - wo,
                     dwody = -ray.ryDirection - wo;
            Float dDNdx = Dot(dwodx, ns) + Dot(wo, dndx);
            Float dDNdy = Dot(dwody, ns) + Dot(wo, dndy);
            rd.rxDirection =
                wi - dwodx + 2.f * Vector3f(Dot(wo, ns) * dndx + dDNdx * ns);
            rd.ryDirection =
                wi - dwody + 2.f * Vector3f(Dot(wo, ns) * dndy + dDNdy * ns);
        }
        return f * Li(rd, scene, sampler, arena, depth + 1) * AbsDot(wi, ns) /
               pdf;
    } else
        return Spectrum(0.f);
}

Spectrum SamplerIntegrator::SpecularTransmit(
    const RayDifferential &ray, const SurfaceInteraction &isect,
    const Scene &scene, Sampler &sampler, MemoryArena &arena, int depth) const {
    Vector3f wo = isect.wo, wi;
    Float pdf;
    const Point3f &p = isect.p;
    const Normal3f &ns = isect.shading.n;
    const BSDF &bsdf = *isect.bsdf;
    Spectrum f = bsdf.Sample_f(wo, &wi, sampler.Get2D(), &pdf,
                               BxDFType(BSDF_TRANSMISSION | BSDF_SPECULAR));
    Spectrum L = Spectrum(0.f);
    if (pdf > 0.f && !f.IsBlack() && AbsDot(wi, ns) != 0.f) {
        // Compute ray differential _rd_ for specular transmission
        RayDifferential rd = isect.SpawnRay(wi);
        if (ray.hasDifferentials) {
            rd.hasDifferentials = true;
            rd.rxOrigin = p + isect.dpdx;
            rd.ryOrigin = p + isect.dpdy;

            Float eta = bsdf.eta;
            Vector3f w = -wo;
            if (Dot(wo, ns) < 0) eta = 1.f / eta;

            Normal3f dndx = isect.shading.dndu * isect.dudx +
                            isect.shading.dndv * isect.dvdx;
            Normal3f dndy = isect.shading.dndu * isect.dudy +
                            isect.shading.dndv * isect.dvdy;

            Vector3f dwodx = -ray.rxDirection - wo,
                     dwody = -ray.ryDirection - wo;
            Float dDNdx = Dot(dwodx, ns) + Dot(wo, dndx);
            Float dDNdy = Dot(dwody, ns) + Dot(wo, dndy);

            Float mu = eta * Dot(w, ns) - Dot(wi, ns);
            Float dmudx =
                (eta - (eta * eta * Dot(w, ns)) / Dot(wi, ns)) * dDNdx;
            Float dmudy =
                (eta - (eta * eta * Dot(w, ns)) / Dot(wi, ns)) * dDNdy;

            rd.rxDirection =
                wi + eta * dwodx - Vector3f(mu * dndx + dmudx * ns);
            rd.ryDirection =
                wi + eta * dwody - Vector3f(mu * dndy + dmudy * ns);
        }
        L = f * Li(rd, scene, sampler, arena, depth + 1) * AbsDot(wi, ns) / pdf;
    }
    return L;
}

Spectrum
NullPathIntegrator::CalculateMISWeight(std::vector<double> P_OVER_F, Spectrum L) {
  if (L.IsBlack()) return Spectrum(0.0f);
  Spectrum RGB(0.0f);
  for (int i = 0; i < 3; i++) {
    double sum = 0.0f;
    for (int j = 0; j < 3; j++) {
      sum += P_OVER_F[i * 3 + j] * (1.0f / Spectrum::nSamples); 
    }
    RGB[i] = sum == 0.0f ? 0.0f : L[i] / sum;
  }
  return RGB;
};

MediumInteraction SampleDistEquiangular(const Scene &scene, 
                            MemoryArena &arena,
														Point3f lightPos,
														Ray ray,
                            Float aRaySurfaceDist,
														Sampler &sampler,
														Float &pdf) {
	// Determine light position and ray direction.
	Point3f aLightPosition = lightPos;
	Vector3f aRayDir = Normalize(ray.d);
	Point3f aRayOrigin = ray.o;

  //Float distance = Lerp(sampler.Get1D(), 0.0f, aRaySurfaceDist);
  //pdf = 1.0f / aRaySurfaceDist;

  // Calculate distance from light plane.
	Vector3f orgToLight = aLightPosition - aRayOrigin;
	float distToProjectionPoint = Dot(aRayDir, orgToLight);
	float projectionDist = (orgToLight - (distToProjectionPoint * aRayDir)).Length();

  assert(projectionDist >= 0.0f);

  float thetaSurfacePoint = std::atan2(aRaySurfaceDist - distToProjectionPoint, 
                                       projectionDist); 
   
  // Determine range of integration.
	double thetaMin, thetaMax;
  if (distToProjectionPoint > 0) {
    thetaMin = -std::atan2(distToProjectionPoint, projectionDist);
	  thetaMax = thetaSurfacePoint;
  } else {
    thetaMin = -thetaSurfacePoint;
    thetaMax = -std::atan2(-distToProjectionPoint, projectionDist);
  }
	double thetaRange = thetaMax - thetaMin;

  if (distToProjectionPoint < 0.0f) assert(thetaMin <= 0.0f && thetaMax <= 0.0f);

  //Determine distance t.
  double sampledTheta = Lerp(sampler.Get1D(), thetaMin, thetaMax);
  double tanSampledTheta = std::tan(sampledTheta);
  double t = distToProjectionPoint < 0.0f
    ? projectionDist * -tanSampledTheta
    : projectionDist * tanSampledTheta;

  double sampledDistance = distToProjectionPoint + t;

  // Calculate the PDF.
	double distance = std::fmin(
      std::fmax(ShadowEpsilon, sampledDistance), 
      ray.tMax * ray.d.Length()
  );
  double dist2 = SQR(distToProjectionPoint - distance) + SQR(projectionDist);
  pdf = projectionDist / (thetaRange * dist2);

  if (distance <= 0.0f) {
    std::cout << "thetaMin="<< thetaMin << std::endl;
    std::cout << "thetaMax="<< thetaMax << std::endl;
    std::cout << "sampledTheta="<< sampledTheta << std::endl;
    std::cout << "sampledDistance="<< sampledDistance << std::endl;
    std::cout << "t=" << t << std::endl;
    std::cout << "projectionDist=" << projectionDist << std::endl;
    std::cout << "distToProjectionPoint=" << distToProjectionPoint << std::endl;
  }
  assert(distance > 0.0f);

  if (distance > ray.tMax * ray.d.Length()) {
    std::cout << "\nthetaMin="<< thetaMin << std::endl;
    std::cout << "thetaMax="<< thetaMax << std::endl;
    std::cout << "sampledDistance="<< sampledDistance << std::endl;
    std::cout << "sampledTheta="<< sampledTheta << std::endl;
    std::cout << "tanSampledTheta="<<tanSampledTheta << std::endl;
    std::cout << "t=" << t << std::endl;
    std::cout << "projectionDist=" << projectionDist << std::endl;
    std::cout << "distToProjectionPoint=" << distToProjectionPoint << std::endl;
    std::cout << "rayTMax=" << ray.tMax << std::endl;
    std::cout << "rayDLength="<< ray.d.Length() << std::endl;
  }
  assert(distance <= ray.tMax * ray.d.Length());
  assert(pdf > 0.0f);
  
  // Determine new position in medium.
  Float g = ray.medium->GetG();
  Point3f scatterPos = ray.o + Normalize(ray.d) * distance;
  MediumInteraction mi(scatterPos, -ray.d, ray.time, ray.medium, 
                       ARENA_ALLOC(arena, HenyeyGreenstein)(g));


  return mi;
}

Float DistEquiangularPDF(const Scene &scene, 
                          Point3f lightPos,
                          Ray ray,
                          float aRaySurfaceDist,
                          float distance) {

	// Determine light position and ray direction.
	Point3f aLightPosition = lightPos;
	Vector3f aRayDir = Normalize(ray.d);
	Point3f aRayOrigin = ray.o;

  // Calculate distance from light plane.
	Vector3f orgToLight = aLightPosition - aRayOrigin;
	float distToProjectionPoint = Dot(aRayDir, orgToLight);
	float projectionDist = (orgToLight - (distToProjectionPoint* aRayDir)).Length();
  float thetaSurfacePoint = std::atan2(aRaySurfaceDist - distToProjectionPoint, 
                                       projectionDist); 
	double thetaMin, thetaMax;
  if (distToProjectionPoint > 0) {
    thetaMin = -std::atan2(distToProjectionPoint, projectionDist);
	  thetaMax = thetaSurfacePoint;
  } else {
    thetaMin = -thetaSurfacePoint;
    thetaMax = -std::atan2(-distToProjectionPoint, projectionDist);
  }
	double thetaRange = thetaMax - thetaMin;

  // Calculate the PDF.
  double dist2 = SQR(distToProjectionPoint - distance) + SQR(projectionDist);
  Float pdf = projectionDist / (thetaRange * dist2);
  assert(pdf > 0.0f);
  return pdf;
}

Float LightPdf(const Scene &scene, SurfaceInteraction isect,
							 Vector3f *wi) {
	Float pdf_ems = 0.0f;
  VisibilityTester vis;
	int nLights = int(scene.lights.size());
	for(const auto &light : scene.lights) {
    bool intersectsLight = light->Intersects(isect, *wi);
    if (!intersectsLight) continue;
		pdf_ems += light->Pdf_Li(isect, *wi) / Float(nLights);
	}
	return pdf_ems;
}
Float LightPdf(const Scene &scene, MediumInteraction mi,
							 Vector3f *wi) {
	Float pdf_ems = 0.0f;
  VisibilityTester vis;
	int nLights = int(scene.lights.size());
	for(const auto &light : scene.lights) {
    bool intersectsLight = light->Intersects(mi, *wi);
    if (!intersectsLight) continue;
		pdf_ems += light->Pdf_Li(mi, *wi) / Float(nLights);
	}
	return pdf_ems;
}

// Used to determine the probability a light and distance were sampled.
// Does not include the null paths.
Float EALightPdf(const Scene &scene, 
                 std::shared_ptr<AreaLight> light,
                 SurfaceInteraction lightIsect,
                 Ray toScatter, 
                 Interaction scatterPoint) {

    // Probability of selecting light and light position.
    Float lightPdf;
    int nLights = int(scene.lights.size());
    std::shared_ptr<Shape> shape = light->GetShape();
    lightPdf = shape->IsSphere()
      ? 1.0f / Float(nLights)
      : 1.0f / (shape->Area() * Float(nLights));
  
    // Probability of selecting distance given light position.
    Point3f lightPos = shape->IsSphere()
      ? (*shape->ObjectToWorld)(Point3f(0, 0, 0)) // center of sphere
      : lightIsect.p;
    Float distPdf = DistEquiangularPDF(scene,
                                       lightIsect.p,
                                       toScatter,
                                       (toScatter.o - scatterPoint.p).Length(),
                                       scatterPoint.IsSurfaceInteraction());
    return lightPdf * distPdf;
}

double SampleLightDirection(const Scene &scene, MediumInteraction mi,
                           Sampler &sampler, Vector3f *wi) {
  Float pdf_ems;
  VisibilityTester vis;
  Spectrum L;
  // uniform random selection of light in scene.
  int nLights = int(scene.lights.size());
  const std::shared_ptr<Light>& light = scene.lights[(int)(sampler.Get1D() * nLights)];
  L = light->Sample_Li(mi, sampler.Get2D(), wi, &pdf_ems, &vis);
  return LightPdf(scene, mi, wi);
}

double MISPhaseLightDirection(
    const Scene &scene, MediumInteraction mi,
    Sampler &sampler,  Vector3f wo, Vector3f *wi, bool &sampledLight) {
  double pdf_phase = 0.0f;
  double pdf_ems = 0.0f;
  if (sampler.Get1D() < 0.5f) {
    sampledLight = false;
    pdf_phase = mi.phase->Sample_p(wo, wi, sampler.Get2D());
    pdf_ems = LightPdf(scene, mi, wi);
  } else {
    sampledLight = true;
    pdf_ems = SampleLightDirection(scene, mi, sampler, wi);
    pdf_phase = mi.phase->p(wo,*wi);
  }
  return 0.5 * (pdf_ems + pdf_phase);
}

Spectrum SampleSpectralNextEvent(
    const Scene &scene, MediumInteraction mi, 
    Sampler &sampler, Ray ray, MemoryArena &arena, 
    Spectrum &W, Spectrum &W_NEE, int channel) {
    
    // Sample a light source  
    Vector3f wi;
    VisibilityTester vis;
    Float pdf_ems;
    int nLights = int(scene.lights.size());
    const std::shared_ptr<Light>& light = scene.lights[
      (int)(sampler.Get1D() * nLights)
    ];
    Spectrum L = light->Sample_Li(mi, sampler.Get2D(), &wi, &pdf_ems, &vis);
    pdf_ems /= nLights;

    W_NEE *= mi.phase->p(-ray.d, wi) / pdf_ems;

    // Ray pointing to target location.
    RayDifferential nee_ray = mi.SpawnRayTo(vis.P1());
    for(int bounces = 0;; bounces++) {
      SurfaceInteraction isect;
      bool foundIntersection = scene.Intersect(nee_ray, &isect);
        
      // Move forward in the medium (if inside of a medium).
      MediumInteraction nee_mi;
      Float M = 0.0f;
      if (nee_ray.medium) {
        Spectrum majorants = nee_ray.medium->GetMajorant();
        int sampleChannel;
        for (int i=0; i < Spectrum::nSamples; i++){
          if (M < majorants[i]) {
            M = majorants[i];
            sampleChannel = i;
          }
        }
       
        if (channel >= 0) {
          M = majorants[channel];
          sampleChannel = channel;
        }

        ray.medium->SampleChannel(nee_ray, sampler, arena, &nee_mi, sampleChannel);
      }

      // Determine wheteher to continue tracking.
      if (nee_mi.IsValid()) {
				Spectrum scatter, absorption, null;
        ray.medium->GetCoefficients(nee_mi.p, scatter, absorption, null);

        if (null.IsBlack()) return Spectrum(0.0f);

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
  
        // If sample is not a null event, we absorbed or scattered (return 0).
				if (sampler.Get1D() > pNull)  return Spectrum(0.0f);

				W *= (Spectrum(M) - scatter - absorption)  / (M * pNull);
        W_NEE *= (Spectrum(M) - scatter - absorption)  / (M * pNull);

				nee_ray = nee_mi.SpawnRayTo(vis.P1());
      } else {
        if (!foundIntersection) return L;
        if (isect.primitive->GetMaterial() != nullptr) return Spectrum(0.0f);
        nee_ray = isect.SpawnRayTo(vis.P1());   
      }
    }
}

Spectrum SampleSpectralNextEvent(
    const Scene &scene, SurfaceInteraction isect, 
    Sampler &sampler, Ray ray, MemoryArena &arena, 
    Spectrum &W, Spectrum &W_NEE, int channel) {
    
    // Sample a light source  
    Vector3f wi;
    VisibilityTester vis;
    Float pdf_ems;
    int nLights = int(scene.lights.size());
    const std::shared_ptr<Light>& light = scene.lights[
      (int)(sampler.Get1D() * nLights)
    ];
    Spectrum L = light->Sample_Li(isect, sampler.Get2D(), &wi, &pdf_ems, &vis);
    pdf_ems /= nLights;

    BxDFType flags;
    Spectrum bsdfFunc = isect.bsdf->f(-ray.d, wi, flags)
                        * AbsDot(wi, isect.shading.n);
    Spectrum bsdfPdf = isect.bsdf->Pdf(-ray.d, wi, flags);

		if (pdf_ems == 0.0f || bsdfFunc.IsBlack()) return Spectrum(0.0f);

    W_NEE *= (bsdfFunc / pdf_ems); 
    W *= (bsdfFunc / bsdfPdf);

    // Ray pointing to target location.
    RayDifferential nee_ray = isect.SpawnRayTo(vis.P1());
    for(int bounces = 0;; bounces++) {
      SurfaceInteraction nee_isect;
      bool foundIntersection = scene.Intersect(nee_ray, &nee_isect);
        
      // Move forward in the medium (if inside of a medium).
      MediumInteraction nee_mi;
      Float M = 0.0f;
      if (nee_ray.medium) {
        Spectrum majorants = nee_ray.medium->GetMajorant();
        int sampleChannel;
        for (int i=0; i < Spectrum::nSamples; i++){
          if (M < majorants[i]) {
            M = majorants[i];
            sampleChannel = i;
          }
        }
       
        if (channel >= 0) {
          M = majorants[channel];
          sampleChannel = channel;
        }

        ray.medium->SampleChannel(nee_ray, sampler, arena, &nee_mi, sampleChannel);
      }

      // Determine wheteher to continue tracking.
      if (nee_mi.IsValid()) {
				Spectrum scatter, absorption, null;
        ray.medium->GetCoefficients(nee_mi.p, scatter, absorption, null);

        if (null.IsBlack()) return Spectrum(0.0f);

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
  
        // If sample is not a null event, we absorbed or scattered (return 0).
				if (sampler.Get1D() > pNull)  return Spectrum(0.0f);

				W *= (Spectrum(M) - scatter - absorption)  / (M * pNull);
        W_NEE *= (Spectrum(M) - scatter - absorption)  / (M * pNull);

				nee_ray = nee_mi.SpawnRayTo(vis.P1());
      } else {
        if (!foundIntersection) return L;
        if (nee_isect.primitive->GetMaterial() != nullptr) return Spectrum(0.0f);
        nee_ray = nee_isect.SpawnRayTo(vis.P1());   
      }
    }
}

Spectrum SampleNextEvent(
    const Scene &scene, SurfaceInteraction isect, 
    Sampler &sampler, Ray ray, int hero, MemoryArena &arena, 
    std::vector<double> &P_OVER_F, std::vector<double> &P_OVER_F_NEE,
    bool equalizeMajorant,
    bool deltaTrack) {
 
		BxDFType bsdfFlags = BSDF_ALL;
   
    // For each lambda there are two variations, NEE and NP.
    P_OVER_F = std::vector<double>(9, 1.0f);
    P_OVER_F_NEE = std::vector<double>(9, 1.0f);

    // uniform random selection of light in scene.
    Vector3f wi;
    VisibilityTester vis;
    Float pdf_ems;
    int nLights = int(scene.lights.size());
    const std::shared_ptr<Light>& light = scene.lights[
      (int)(sampler.Get1D() * nLights)
    ];
    Spectrum L = light->Sample_Li(isect, sampler.Get2D(), &wi, &pdf_ems, &vis);
    pdf_ems /= nLights;

    Spectrum bsdfFunc = isect.bsdf->f(isect.wo, wi, bsdfFlags) * AbsDot(wi, isect.shading.n);
    Spectrum bsdfPdf = isect.bsdf->Pdf(isect.wo, wi, bsdfFlags);

		if (pdf_ems == 0.0f || bsdfFunc.IsBlack()) return Spectrum(0.0f);

    // Update path throughput with phase and with directional pdf. 
    UpdateWeights(P_OVER_F, bsdfFunc, bsdfPdf);
    UpdateWeights(P_OVER_F_NEE, bsdfFunc, pdf_ems);

    // Ray pointing to target location.
    RayDifferential nee_ray = isect.SpawnRayTo(vis.P1());
    for(int bounces = 0;; bounces++) {
      SurfaceInteraction nee_isect;
      bool foundIntersection = scene.Intersect(nee_ray, &nee_isect);
        
      // Move forward in the medium (if inside of a medium).
      MediumInteraction nee_mi;
      Spectrum Tr(1.0f), majorants;
      if (nee_ray.medium) {
        majorants = nee_ray.medium->GetMajorant();
        if (equalizeMajorant) {
          int sampleChannel;
          majorants = HomogenizeMajorants(majorants, sampleChannel);
          Tr = SampleHomogenizedMajorants(nee_ray, sampler, arena, &nee_mi, 
                                          sampleChannel);
        } else {
          Tr = nee_ray.medium->SampleChannel(nee_ray, sampler, arena, 
                                             &nee_mi, hero);
        }
        Spectrum freeFlightPdf = nee_mi.IsValid() ? majorants * Tr : Tr;
        UpdateWeights(P_OVER_F, Tr, freeFlightPdf);
        UpdateWeights(P_OVER_F_NEE, Tr, freeFlightPdf);
      }

      // Determine whether to continue tracking.
      if (nee_mi.IsValid()) {
        // null vertex with forward scattering.
        // Get the scattering and absorption coefficients. 
        Spectrum scatter, absorption, null;
        nee_ray.medium->GetCoefficients(nee_mi.p, scatter, absorption, null);
        
        if (null.IsBlack()) return Spectrum(0.0f); 

        if (deltaTrack) {
          Spectrum pNull =  null[hero] / majorants[hero];
          if (sampler.Get1D() > null[hero] / majorants[hero])return Spectrum(0.0f);
          UpdateWeights(P_OVER_F, null, null / majorants);
          UpdateWeights(P_OVER_F_NEE, null, null / majorants);
        } else {
          UpdateWeights(P_OVER_F, null, null / majorants);
          UpdateWeights(P_OVER_F_NEE, null, Spectrum(1.0f));
        }
        nee_ray = nee_mi.SpawnRayTo(vis.P1());
      } else {
        if (!foundIntersection) return L;
        if (nee_isect.primitive->GetMaterial() != nullptr) return Spectrum(0.0f);
        nee_ray = nee_isect.SpawnRayTo(vis.P1());   
      }
    }
}

Spectrum SampleNextEvent(
    const Scene &scene, MediumInteraction mi, 
    Sampler &sampler, Ray ray, int hero, MemoryArena &arena, 
    std::vector<double> &P_OVER_F, std::vector<double> &P_OVER_F_NEE,
    bool equalizeMajorant,
    bool deltaTrack) {
    
    // For each lambda there are two variations, NEE and NP.
    P_OVER_F = std::vector<double>(9, 1.0f);
    P_OVER_F_NEE = std::vector<double>(9, 1.0f);

    // uniform random selection of light in scene.
    Vector3f wi;
    VisibilityTester vis;
    Float pdf_ems;
    int nLights = int(scene.lights.size());
    const std::shared_ptr<Light>& light = scene.lights[
      (int)(sampler.Get1D() * nLights)
    ];
    Spectrum L = light->Sample_Li(mi, sampler.Get2D(), &wi, &pdf_ems, &vis);
    pdf_ems /= nLights;

    // Update path throughput with phase and with directional pdf. 
    UpdateWeights(P_OVER_F, 
                  mi.phase->p(-ray.d, wi), 
                  mi.phase->p(-ray.d, wi));
  
    UpdateWeights(P_OVER_F_NEE, 
                  mi.phase->p(-ray.d, wi), 
                  pdf_ems);

    // Ray pointing to target location.
    RayDifferential nee_ray = mi.SpawnRayTo(vis.P1());
    for(int bounces = 0;; bounces++) {
      SurfaceInteraction isect;
      bool foundIntersection = scene.Intersect(nee_ray, &isect);
        
      // Move forward in the medium (if inside of a medium).
      MediumInteraction nee_mi;
      Spectrum Tr(1.0f), majorants;
      if (nee_ray.medium) {
        majorants = nee_ray.medium->GetMajorant();
        if (equalizeMajorant) {
          int sampleChannel;
          majorants = HomogenizeMajorants(majorants, sampleChannel);
          Tr = SampleHomogenizedMajorants(nee_ray, sampler, arena, &nee_mi, 
                                          sampleChannel);
        } else {
          Tr = nee_ray.medium->SampleChannel(nee_ray, sampler, arena, 
                                             &nee_mi, hero);
        }
        Spectrum freeFlightPdf = nee_mi.IsValid() ? majorants * Tr : Tr;
        UpdateWeights(P_OVER_F, Tr, freeFlightPdf);
        UpdateWeights(P_OVER_F_NEE, Tr, freeFlightPdf);
      }

      // Determine wheteher to continue tracking.
      if (nee_mi.IsValid()) {
        // null vertex with forward scattering.
        // Get the scattering and absorption coefficients. 
        Spectrum scatter, absorption, null;
        nee_ray.medium->GetCoefficients(nee_mi.p, scatter, absorption, null);
        
        if (null.IsBlack()) return Spectrum(0.0f); 

        if (deltaTrack) {
          Spectrum pNull =  null[hero] / majorants[hero];
          if (sampler.Get1D() > null[hero] / majorants[hero])return Spectrum(0.0f);
          UpdateWeights(P_OVER_F, null, null / majorants);
          UpdateWeights(P_OVER_F_NEE, null, null / majorants);
        } else {
          UpdateWeights(P_OVER_F, null, null / majorants);
          UpdateWeights(P_OVER_F_NEE, null, Spectrum(1.0f));
        }
        nee_ray = nee_mi.SpawnRayTo(vis.P1());
      } else {
        if (!foundIntersection) return L;
        if (isect.primitive->GetMaterial() != nullptr) return Spectrum(0.0f);
        nee_ray = isect.SpawnRayTo(vis.P1());   
      }
    }
}

bool CalculateNullPath(const Scene &scene, 
                       Sampler &sampler, 
                       MemoryArena &arena,
                       int hero,
                       Interaction &start,
                       Interaction &stop,
                       std::vector<double> &P_OVER_F,
                       std::vector<double> &P_OVER_F_FORCED_NULL,
                       std::string &path,
                       bool equalizeMajorant,
                       bool printPath) { 
  // For each lambda there are two variations, NEE and NP.
  P_OVER_F = std::vector<double>(9, 1.0f);
  P_OVER_F_FORCED_NULL = std::vector<double>(9, 1.0f);

  // Ray pointing to target location.
  RayDifferential ray = start.SpawnRayTo(stop);

  assert(start.p != stop.p);

  std::stringstream symbPath;
  for (int bounces = 0;; bounces++) {
    SurfaceInteraction isect;
    bool foundIntersection = scene.Intersect(ray, &isect);

    MediumInteraction mi;
    Spectrum Tr(1.0f), majorants; 
    if (ray.medium) {
      majorants = ray.medium->GetMajorant();
      if (equalizeMajorant) {
        int sampleChannel;
        majorants = HomogenizeMajorants(majorants, sampleChannel);
        Tr = SampleHomogenizedMajorants(ray, sampler, arena, &mi, 
                                        sampleChannel);
      } else {
        Tr = ray.medium->SampleChannel(ray, sampler, arena, 
                                       &mi, hero);
      }
      Spectrum freeFlightPdf = mi.IsValid() ? majorants * Tr : Tr;
      UpdateWeights(P_OVER_F, Tr, freeFlightPdf);
      UpdateWeights(P_OVER_F_FORCED_NULL, Tr, freeFlightPdf);
    }

    // Determine wheteher to continue tracking.
    if (mi.IsValid()) {
      // null vertex with forward scattering.
      // Get the scattering and absorption coefficients. 
      Spectrum scatter, absorption, null;
      ray.medium->GetCoefficients(mi.p, scatter, absorption, null);
      
      if (null.IsBlack()) return false; 

      UpdateWeights(P_OVER_F, null, null / majorants);
      UpdateWeights(P_OVER_F_FORCED_NULL, null, Spectrum(1.0f));
     
      ray = mi.SpawnRayTo(stop);   
    } else {
      if (!foundIntersection) return true;
      else if (isect.primitive->GetMaterial() != nullptr) return false;
      ray = isect.SpawnRayTo(stop);
    }
  }
}

Spectrum SampleNEE(const Scene &scene,
    Sampler &sampler, 
    MediumInteraction mi,
    Ray ray, 
    Ray lastRealRay,
    int hero,
    MemoryArena &arena, 
    std::vector<double> &P_OVER_F_NEE,
    std::vector<double> &P_OVER_F_EA,
    Spectrum scatterAtHitPoint,
    bool equalizeMajorant,
    bool printPath) {


  P_OVER_F_NEE = std::vector<double>(9, 1.0f); 
  P_OVER_F_EA = std::vector<double>(9, 1.0f); 
  
  // Choose a light source (not from mixture density).
  // =========================================================//
  // uniform random selection of light in scene.
  int nLights = int(scene.lights.size());
	const std::shared_ptr<Light>& light = scene.lights[
		std::min((int)(sampler.Get1D() * nLights), nLights - 1)];

	// Currently EA only supports area lights, so have to respect that for NEE. 
	const std::shared_ptr<AreaLight> areaLight 
		= std::dynamic_pointer_cast<AreaLight>(light);		
	assert(areaLight);

  Vector3f wi; 
  Float lightDirPdf;
  VisibilityTester vis;
  Spectrum L;
  L = light->Sample_Li(mi, sampler.Get2D(), &wi, &lightDirPdf, &vis);
  Float scatterFunc = mi.phase->p(-ray.d, wi);

  // Construct a path of all null events from scatter position to light.
  // =========================================================//
  std::vector<double> P_OVER_F_SCATTER_TO_LIGHT, P_OVER_F_NEE_SCATTER_TO_LIGHT;
  std::string symbSubPath;
  Interaction lightPos = vis.P1();
  bool reachedLight = CalculateNullPath(scene, sampler, arena, hero,
                        mi, 
                        lightPos,
                        P_OVER_F_SCATTER_TO_LIGHT,
                        P_OVER_F_NEE_SCATTER_TO_LIGHT,
                        symbSubPath,
                        false,
                        printPath);
  if (!reachedLight) {
    if (printPath) std::cout << "Path didn't make it from scatter to light!" << std::endl;
    return Spectrum(0.0f);
  }
  
  // Determine the probability of the shape of this EA path.
  // =========================================================//
  Float sampledDistance = (mi.p - lastRealRay.o).Length();
  Float aRaySurfaceDist = lastRealRay.tMax * lastRealRay.d.Length();
  Float distPdf = DistEquiangularPDF(scene, lightPos.p, lastRealRay, 
                               aRaySurfaceDist, sampledDistance);
  Float lightPdf = lightDirPdf / nLights;

  // NEE
  UpdateWeights(P_OVER_F_NEE, scatterFunc, lightPdf);
  MultiplyVec(P_OVER_F_NEE, P_OVER_F_NEE_SCATTER_TO_LIGHT);

  // EA
  MultiplyVec(P_OVER_F_EA, P_OVER_F_NEE_SCATTER_TO_LIGHT); 
  UpdateWeights(P_OVER_F_EA, scatterAtHitPoint * scatterFunc, distPdf * lightPdf); 

  return L; 
}

Spectrum SampleEA(
    const Scene &scene,
    Sampler &sampler, 
    Ray ray, 
    int hero,
    MemoryArena &arena, 
    std::vector<double> &P_OVER_F_NEE,
    std::vector<double> &P_OVER_F_EA,
    bool equalizeMajorant,
    bool printPath) {
  std::stringstream symbPath;
  symbPath << "\\frac{f(\\overline{x})}{p(\\overline{x})}=\n";

  P_OVER_F_NEE = std::vector<double>(9, 1.0f);
  P_OVER_F_EA = std::vector<double>(9, 1.0f);

  Ray ea_ray = ray;

  // Choose a distance
  // ===========================================================================

  // Choose a light source for distance sample
	int nLights = int(scene.lights.size());
	const std::shared_ptr<Light>& light = scene.lights[
		std::min((int)(sampler.Get1D() * nLights), nLights - 1)];

	// Currently EA only supports area lights. 
	const std::shared_ptr<AreaLight> areaLight 
		= std::dynamic_pointer_cast<AreaLight>(light);		

	assert(areaLight);

	Float l;
  std::shared_ptr<Shape> shape = areaLight->GetShape();
  Point3f lightPoint = (*shape->ObjectToWorld)(Point3f(0, 0, 0));

  // Intersect ray with the scene.
  SurfaceInteraction isect;
  bool foundIntersection = scene.Intersect(ea_ray, &isect);
  if (!ea_ray.medium) return Spectrum(0.0f);

  Float aRaySurfaceDist = (isect.p - ea_ray.o).Length(); 
  MediumInteraction mi(ea_ray.o, -ea_ray.d, ea_ray.time, ea_ray.medium,
    ARENA_ALLOC(arena, HenyeyGreenstein)(ea_ray.medium->GetG()));

  // Sample a distance in the direction of the ray, relative to sampled light pos. 
	Float distPdf;
	MediumInteraction scatterPos = SampleDistEquiangular(scene, arena, 
                                                       lightPoint, ea_ray, 
                                                       aRaySurfaceDist, 
                                                       sampler, distPdf);

  // Make sure the scattering event is not happening in a zero density region.
  Spectrum majorants = ea_ray.medium->GetMajorant();
  Spectrum scatter = ea_ray.medium->GetScattering(scatterPos.p);

  if (scatter.IsBlack()) return Spectrum(0.0f);

  // Choose a direction towards a light.
  // ===========================================================================
  VisibilityTester lightVis;
  Vector3f dirToLight;
  Float lightDirPdf;
  Spectrum L = light->Sample_Li(scatterPos, sampler.Get2D(), &dirToLight, &lightDirPdf, 
                                &lightVis);
 
  Float scatterFunc = mi.phase->p(-ea_ray.d, dirToLight);
    // Calculate transmittance between origin and scatter position.
  // =============================================================================
  std::vector<double> P_OVER_F_REF_TO_SCATTER, P_OVER_F_EA_REF_TO_SCATTER;
  std::string symbSubPath1;
  bool reachedScatter = CalculateNullPath(scene, sampler, arena, hero, 
                       mi, scatterPos,
                       P_OVER_F_REF_TO_SCATTER,
                       P_OVER_F_EA_REF_TO_SCATTER,
                       symbSubPath1,
                       false,
                       printPath);
 
  if (!reachedScatter) {
    if (printPath) std::cout << "Path didn't make it from origin to scatter!" << std::endl;
    return Spectrum(0.0f);
  }
  // Calculate transmittance between scatter position and light.
  // =============================================================================
  std::vector<double> P_OVER_F_SCATTER_TO_LIGHT, P_OVER_F_EA_SCATTER_TO_LIGHT;
  std::string symbSubPath2;
  Interaction lightPos = lightVis.P1();
  bool reachedLight = CalculateNullPath(scene, sampler, arena, hero,
                        scatterPos, lightPos,
                        P_OVER_F_SCATTER_TO_LIGHT,
                        P_OVER_F_EA_SCATTER_TO_LIGHT,
                        symbSubPath2,
                        false,
                        printPath);
  if (!reachedLight) {
    if (printPath) std::cout << "Path didn't make it from scatter to light!" << std::endl;
    return Spectrum(0.0f);
  }

  Float lightPdf = lightDirPdf / nLights;

  // NEE PATH
	UpdateWeights(P_OVER_F_NEE, scatterFunc, lightPdf);
	MultiplyVec(P_OVER_F_NEE, P_OVER_F_EA_SCATTER_TO_LIGHT);
	
  // EA PATH
  MultiplyVec(P_OVER_F_EA, P_OVER_F_EA_REF_TO_SCATTER); // Tr to scatter.
  MultiplyVec(P_OVER_F_EA, P_OVER_F_EA_SCATTER_TO_LIGHT); // Tr to light.
  UpdateWeights(P_OVER_F_EA, scatter * scatterFunc, distPdf * lightPdf);

  if (printPath) {
    Point3f p0 = ea_ray.o;                 // The point where the EA sampling starts.
    Point3f p1 = scatterPos.p;          // The scatter point.
    Point3f p2 = lightPoint;            // The sampled light position.
    Point3f p3 = lightVis.P1().p;            // Intersection with light.
    Float dist = mi.SpawnRayTo(scatterPos.p).d.Length(); // The sampled distance.

    // Float distPdf                       The probability of the sampled distance.
    symbPath << "points:" << std::endl;
    symbPath << "p0=(" << p0[0] << ", " << p0[1] << ", " << p0[2] << ")" << std::endl;
    symbPath << "p1=(" << p1[0] << ", " << p1[1] << ", " << p1[2] << ")" << std::endl;
    symbPath << "p2=(" << p2[0] << ", " << p2[1] << ", " << p2[2] << ")" << std::endl;
    symbPath << "p3=(" << p3[0] << ", " << p3[1] << ", " << p3[2] << ")" << std::endl;
    symbPath << "distances:" << std::endl;
    symbPath << "d(p0,p1)=" << dist << std::endl;
    symbPath << "d(p1,p3)=" << (lightPos.p - scatterPos.p).Length() << std::endl;
    symbPath << "pdfs:" << std::endl;
    symbPath << "pdf(light)" << lightPdf << std::endl;
    symbPath << "pdf(d(p0,p1))=" << distPdf << std::endl;
    symbPath << "scatterFunc=" << scatterFunc << std::endl;
    symbPath << "L==(" << L[0] << ", " << L[1] << ", " << L[2] << ")" << std::endl;
    symbPath << "path(p0,p1)=\n" << symbSubPath1 << std::endl;;
    symbPath << "path(p1,p3)=\n" << symbSubPath2 << std::endl;
    std::cout << symbPath.str() << std::endl;
  }

  return L; 
}

Spectrum HomogenizeMajorants(Spectrum &majorants, int &maxChannel) {
  Float maxMajorant = 0.0f;
  Float M = 0.0f;
  for (int i = 0 ; i < Spectrum::nSamples; i++){
    if (M < majorants[i]) {
      M = majorants[i];
      maxChannel = i;
    }
  }
  majorants = Spectrum(M);
  return majorants;
}

Spectrum SampleHomogenizedMajorants(const Ray &ray, Sampler &sampler, 
                                    MemoryArena &arena, MediumInteraction *mi, 
                                    int channel) {
  Spectrum Tr;
  Tr = ray.medium->SampleChannel(ray, sampler, arena, mi, channel);
  return Spectrum(Tr[channel]);
}

void UpdateWeights(std::vector<double> &P_OVER_F, Spectrum F, Spectrum P) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      P_OVER_F[i * 3 + j] *= F[i] == 0.0f ? 0.0f : P[j] / F[i];
    }
  }
}

std::vector<double> SetWeights(std::vector<double> &vec, double val) {
  for (int i = 0; i < vec.size(); i++) vec[i] = val;
  return vec;
}

std::vector<double> SetWeights(std::vector<double> &vec, std::vector<double> val) {
  for (int i = 0; i < vec.size(); i++) vec[i] = val[i];
  return vec;
}

std::vector<double> MultiplyVec(std::vector<double> &vec, double val) {
  for (int i=0; i < vec.size(); i++) vec[i] *=val;
  return vec;
} 

std::vector<double> MultiplyVec(std::vector<double> &vec, Spectrum s) {
  for (int i=0; i < vec.size(); i++) vec[i] *= s[i];
  return vec;
} 

std::vector<double> MultiplyVec(std::vector<double> &vec, 
                                std::vector<double> val) {
  for (int i=0; i < vec.size(); i++) vec[i] *= val[i];
  return vec;
} 

std::vector<double> DivideVec(std::vector<double> &vec, double val) {
  for (int i=0; i < vec.size(); i++) vec[i] /=val;
  return vec;
} 

std::vector<double> DivideVecVal(const std::vector<double> &vec, double val) {
  std::vector<double> newVec(vec.size());
  for (int i=0; i < vec.size(); i++) newVec[i] = vec[i] / val;
  return newVec;
} 

double AverageVec(const std::vector<double> &vec) {
  double avg = 0.0f;
  for (int i=0; i < vec.size(); i++) avg += vec[i];
  return avg / vec.size();
}

double MaxVal(const std::vector<double> &vec) {
  double maxVal = 0.0f;
  for (int i = 0; i < vec.size(); i++) maxVal = maxVal > vec[i] ? maxVal : vec[i];
  return maxVal;
}

bool AllZerosVec(const std::vector<double> &vec) {
  bool isAllZeros = true;
  for (int i=0; i < vec.size(); i++) isAllZeros = isAllZeros && vec[i] == 0;
  return isAllZeros;
}

std::vector<double> SpecToVec(Spectrum spec) {
  std::vector<double> vec;
  for (int i = 0; i < Spectrum::nSamples; i++) vec.push_back(spec[i]);
  return vec;
}

std::vector<double> SpecToVec(Spectrum spec, 
                              std::vector<double> &vec) {
  for (int i = 0; i < Spectrum::nSamples; i++) vec[i] = spec[i];
  return vec;
}

std::vector<double> FloatToVec(Float f, std::vector<double> &vec) {
  for (int i = 0; i < vec.size(); i++) {
    vec[i] = f;
  }
  return vec;
}

Spectrum VecToSpec(const std::vector<double> &vec) {
  Spectrum s;
  for (int i=0; i < vec.size(); i++) s[i] = vec[i];
  return s;
}

void PrintVec(const std::vector<double> &vec) {
  for (int i = 0; i < vec.size(); i++) {
    std::cout << vec[i] << " ";
  }
  std::cout << std::endl;
}

Spectrum 
FlagFilter(Spectrum s, 
           bool ignoreNoScatter,
					 int viewNumScatterEvents,
           int numScatterEvents) {
  if (ignoreNoScatter) {
    s = numScatterEvents == 0 ? Spectrum(0.0f) : s;
  }

	if (viewNumScatterEvents >= 0) {
		s = viewNumScatterEvents == numScatterEvents ? s : Spectrum(0.0f);	
	}

  return s;
}

// Can be called after taking the most recent scatter event into account.
bool ExitEarly(bool singleScatter, int viewNumScatterEvents, int numScatterEvents) {
  if (singleScatter && numScatterEvents > 1) return true;
  if (viewNumScatterEvents > 0 && numScatterEvents > viewNumScatterEvents) return true;
  return false;
}

}  // namespace pbrt
