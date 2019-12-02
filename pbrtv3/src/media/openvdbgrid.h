#ifndef __OPENVDBGRIDMEDIUM_H__
#define __OPENVDBGRIDMEDIUM_H__

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>

#include "medium.h"
#include "transform.h"
#include "stats.h"
#include "texture.h"

namespace pbrt
{
  typedef openvdb::GridBase Grid;
  typedef openvdb::FloatGrid MonoGrid;
  typedef openvdb::Vec3SGrid RGBGrid;
  typedef openvdb::Vec3DGrid RGBGridf;
  typedef openvdb::tools::GridSampler<MonoGrid, openvdb::tools::BoxSampler> 
    MonoGridSampler;
          

  class OpenvdbGridMedium : public Medium
  {
  public:
    OpenvdbGridMedium(const std::string& filename, 
                      const Spectrum& sigma_a, const Spectrum& sigma_s, 
                      const Spectrum& sigma_n,
                      const Spectrum &sigma_a_2, const Spectrum& sigma_s_2,
                      const Spectrum& sigma_n_2,
                      const Float g, const Float scale, 
                      const Float gamma,
                      const Transform &mediumToWorld,
                      const bool useTemp,
                      const bool usePerlin,
                      const Point3f noiseScale,
                      const bool useLinearRamp,
                      const Float rampScale,
                      const int rampDimension, 
                      const Float rampOrigin);

    Point3i getDims() const;

    Spectrum Density(const Point3f &p) const;
    
    Float Temperature(const Point3f &p) const;

    Spectrum Sample(const Ray &ray, Sampler &sampler, 
                 MemoryArena &arena,
                 MediumInteraction *mi) const;
        
    Spectrum SampleChannel(const Ray &rWorld, Sampler &sampler,   
                    MemoryArena &arena, MediumInteraction *mi, int channel) const;

    Spectrum Tr(const Ray& ray, Sampler& sampler) const;

    Spectrum GetAbsorption(const Point3f &p) const;

    Spectrum GetScattering(const Point3f &p) const;

    void GetCoefficients(const Point3f &p, Spectrum &s, Spectrum &a, Spectrum &n) const;

    Spectrum GetMajorant() const;

    Spectrum InterpolateWithTemperature(Point3f p, Spectrum x, Spectrum y) const;

    Spectrum InterpolateWithPerlin(Point3f p, Spectrum x, Spectrum y) const;

    Float GetG() const { return m_g; }

    bool SetMinAndMaxPosition(Ray &ray) const;

  private:
    bool open(const std::string& filepath);
    
    openvdb::io::File* m_file;
    MonoGrid::Ptr m_density_grid;

    std::string m_filename;
    const Float m_g;
    const Float m_gamma;
    const Transform m_worldToMedium; // assume same for all grids.
    const Spectrum m_sigma_a, m_sigma_s, m_sigma_n;
    const Spectrum m_sigma_a_2, m_sigma_s_2, m_sigma_n_2;
    Spectrum m_sigma_t;

    Spectrum m_maxDensity;
    Point3i m_n;
    Point3i m_min;

    // temperature grid and color values for interpolation.
    const bool m_useTemp;
    Float m_maxTemperature;
    Point3i m_temp_n;
    Point3i m_temp_min;
    MonoGrid::Ptr m_temperature_grid;

    //linear ramp
    bool m_useLinearRamp;
    Float m_rampScale;
    int m_rampDimension;
    Float m_rampOrigin;

    // perlin noise
    bool m_usePerlin;
    Point3f m_noiseScale;
    
    void UpdateAndComputeMaxDensity();
    void ComputeMaxTemperature();
  };
}
#endif
