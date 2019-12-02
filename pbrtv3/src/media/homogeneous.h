
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

#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_MEDIA_HOMOGENEOUS_H
#define PBRT_MEDIA_HOMOGENEOUS_H

// media/homogeneous.h*
#include "medium.h"
#include "imageio.h"
#include <cmath>
#include <cctype>

namespace pbrt {

struct MarbleScanlineOpts {
  std::string filename;
  int colorChannel;
  int x;

  MarbleScanlineOpts(): x(0), filename(""), colorChannel(0) {};
  MarbleScanlineOpts(std::string filename, int x, int colorChannel)
    : filename(filename)
    , x(x)
    , colorChannel(colorChannel) {};
};

struct MarbleScanline {
  std::vector<Float> scanline;
  int length;

  MarbleScanline(MarbleScanlineOpts opts) {
    if (opts.filename.empty()) return;
   
    Point2i res;
    std::unique_ptr<RGBSpectrum[]> img(ReadImageEXR(opts.filename, &res.x, &res.y));
    length = res.y;
    Float maxVal = 0.0f;
    Float minVal = std::numeric_limits<Float>::infinity();
    for (int y = 0; y  < length; y++) {
      Float val = img[opts.x + (y * res.x)][opts.colorChannel];
      scanline.push_back(val);
      if (val > maxVal) maxVal = val;
      if (val < minVal) minVal = val;
    }

    for (int y = 0; y  < length; y++) {
      scanline[y] = (scanline[y] - minVal + 0.01) / (maxVal - minVal + .01);
    }
  }

  Float get(Float y) const {
    y = std::abs(y);

    // Get the pixel that y falls into. 
    int pixel_index = round(y);
    Float pixel_offset = y - pixel_index;

    int lowerPixel, upperPixel;
    bool inLowerHalf = std::abs(pixel_offset) < 0.5f;
    lowerPixel = inLowerHalf ? pixel_index - 1 : pixel_index;
    upperPixel = inLowerHalf ? pixel_index : pixel_index + 1;

    // Determine if in an odd or even segment.
    bool lowerIsOdd = (int)(round(lowerPixel / length)) % 2 != 0;
    bool upperIsOdd = (int)(round(upperPixel / length)) % 2 != 0;

    // Map to the original domain.
    y = pixel_offset + pixel_index % length;
    lowerPixel %= length;
    upperPixel %= length;

    // Flip offset if in an odd segment.
    Float savedLowerPixel = lowerPixel;
    y = lowerIsOdd ? length - y : y;
    lowerPixel = lowerIsOdd ? length - upperPixel : lowerPixel;
    upperPixel = lowerIsOdd ? length - savedLowerPixel : upperPixel;

    // Mixture weights
    Float lowerPixelCenter = lowerPixel + 0.5;
    Float distanceToLower = y - lowerPixelCenter;
    Float lowerWeight = (1.0f - distanceToLower);

    return lowerWeight * scanline[lowerPixel] + (1.0f - lowerWeight) * scanline[upperPixel];
  }
};

// HomogeneousMedium Declarations
class HomogeneousMedium : public Medium {
  public:
    // HomogeneousMedium Public Methods
    HomogeneousMedium(const Spectrum &sigma_a, const Spectrum &sigma_s, 
                      const Spectrum &sigma_n,
                      Float g, bool enableDualDensity = false, 
                      bool enableMarble = false,
                      Point3f noiseScale = Point3f(1.0f,1.0f,1.0f),
                      Float noiseScale2 = 1.0f,
                      Float a = 1.0f,
                      Float w = 1.0f,
                      int m_num_octaves = 1,
                      const Spectrum &sigma_a_2 = Spectrum(0), 
                      const Spectrum &sigma_s_2 = Spectrum(0),
                      const Spectrum &sigma_n_2 = Spectrum(0),
                      const MarbleScanlineOpts marbleOpts = MarbleScanlineOpts());
  
  Spectrum Tr(const Ray &ray, Sampler &sampler, uint32_t flags = 0xFFFFFFFF, 
              TransportMode mode = TransportMode::Radiance) const;

  Spectrum Tr(const RayDifferential& ray, Sampler& sampler, 
              uint32_t flags = 0xFFFFFFFF, 
              TransportMode mode = TransportMode::Radiance) const;
	
  Spectrum Sample(const Ray &ray, Sampler &sampler, 
                  MemoryArena &arena, MediumInteraction *mi, 
                  uint32_t flags = 0xFFFFFFFF, 
                  TransportMode mode = TransportMode::Radiance) const;
  	
  Spectrum SampleChannel(const Ray &_ray, Sampler &sampler, MemoryArena &arena, 
                         MediumInteraction *mi, 
                         int channel, uint32_t flags = 0xFFFFFFFF, 
                         TransportMode mode = TransportMode::Radiance) const;
	
  Spectrum Sample(const RayDifferential &ray, Sampler &sampler, 
                  MemoryArena &arena, MediumInteraction *mi, 
                  uint32_t flags = 0xFFFFFFFF, 
                  TransportMode mode = TransportMode::Radiance) const;

  Spectrum GetAbsorption(const Point3f &p) const;
  Spectrum GetScattering(const Point3f &p) const;
  void GetCoefficients(const Point3f &p, Spectrum &s, Spectrum &a, Spectrum &n) const;
  Spectrum GetMajorant() const;

  Float GetG() const { return g; }

  bool SetMinAndMaxPosition(Ray &ray) const { return true; };

  private:
    // HomogeneousMedium Private Data
    const Spectrum sigma_a, sigma_s, sigma_t;
    const Float g;

    // Dual densities and marble:
    bool m_enableDualDensity;
    bool m_enabeMarble;
    Point3f m_noiseScale;
    Float m_noiseScale2;
    Float m_a;
    Float m_w;
    Float m_num_octaves;
    const Spectrum sigma_a_2, sigma_s_2, sigma_t_2;
    Spectrum max_sigma_t;
    const MarbleScanline marbleScanline;

    Float GetMarble(const Point3f &p) const;
};

}  // namespace pbrt

#endif  // PBRT_MEDIA_HOMOGENEOUS_H
