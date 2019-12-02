#include "openvdbgrid.h"
#include "medium.h"
#include "interaction.h"
#include "sampler.h"
#include <openvdb/tools/Interpolation.h>

namespace pbrt
{

  STAT_COUNTER("Integrator/openvdb lookups", lookups);

  OpenvdbGridMedium::OpenvdbGridMedium(const std::string& filename, 
                                      const Spectrum& sigma_a, 
                                      const Spectrum& sigma_s, 
                                      const Spectrum& sigma_n,
                                      const Spectrum& sigma_a_2, 
                                      const Spectrum& sigma_s_2, 
                                      const Spectrum& sigma_n_2,
                                      const Float g, const Float scale,
                                      const Float gamma,
                                      const Transform& mediumToWorld,
                                      const bool useTemp,
                                      const bool usePerlin,
                                      const Point3f noiseScale,
                                      const bool useLinearRamp,
                                      const Float rampScale,
                                      const int rampDimension,
                                      const Float rampOrigin)
  : m_filename(filename)
  , m_sigma_a(sigma_a*scale)
  , m_sigma_s(sigma_s*scale)
  , m_sigma_n(sigma_n*scale)
  , m_sigma_t(sigma_a*scale + sigma_s*scale + sigma_n*scale)
  , m_sigma_a_2(sigma_a_2*scale)
  , m_sigma_s_2(sigma_s_2*scale)
  , m_sigma_n_2(sigma_n_2*scale)
  , m_g(g)
  , m_gamma(gamma)
  , m_worldToMedium(Inverse(mediumToWorld))
  , m_useTemp(useTemp)
  , m_usePerlin(usePerlin)
  , m_noiseScale(noiseScale)
  , m_useLinearRamp(useLinearRamp)
  , m_rampScale(rampScale)
  , m_rampDimension(rampDimension)
  , m_rampOrigin(rampOrigin)
  {
    openvdb::initialize();
    if (open(m_filename))
    {
      Grid::Ptr baseGrid;
      // get the nx, ny, nz values of the grid
      // apparently openvdb files can have multiple grids inside them
      // right now we support only one grid
      for (openvdb::io::File::NameIterator nameIter = m_file->beginName();
         nameIter != m_file->endName(); ++nameIter)
      {
        std::cout << nameIter.gridName() << std::endl;
        baseGrid = m_file->readGrid(nameIter.gridName());

        if (nameIter.gridName() == "density" && baseGrid->isType<MonoGrid>()) {
          std::cout << "Density Grid Found" << std::endl;
          // Dimensions of grid.
          openvdb::Coord coord = baseGrid->evalActiveVoxelDim();
          m_n = Point3i(coord.x(), coord.y(), coord.z());

          // Minimum Index.
          openvdb::Coord min = baseGrid->evalActiveVoxelBoundingBox().min();
          m_min = Point3i(min.x(), min.y(), min.z());

          m_density_grid = openvdb::gridPtrCast<MonoGrid>(baseGrid);
        } else if (useTemp && nameIter.gridName() == "temperature") {
          std::cout << "Temperature Grid Found" << std::endl;
        
          // Dimensions of grid. 
          openvdb::Coord coord = baseGrid->evalActiveVoxelDim(); 
          m_temp_n = Point3i(coord.x(), coord.y(), coord.z());
         
          // Minimum Index.
          openvdb::Coord min = baseGrid->evalActiveVoxelBoundingBox().min();
          m_temp_min = Point3i(min.x(), min.y(), min.z());

          m_temperature_grid = openvdb::gridPtrCast<MonoGrid>(baseGrid);
        } else {
          continue;
        }

        // Print out all of the meta data associated with the density grid
        std::cout << "Meta Data: " << std::endl;
        for (openvdb::MetaMap::MetaIterator iter = baseGrid->beginMeta();
          iter != baseGrid->endMeta(); ++iter)
        {
          const std::string& name = iter->first;
          openvdb::Metadata::Ptr value = iter->second;
          std::string valueAsString = value->str();
          std::cout << "\t" << name << " = " << valueAsString << std::endl;
        }
      }
    } else {
      Error("Unable to open the openvdb file.\n");
      exit(-1);
    }

    UpdateAndComputeMaxDensity();
    
    if (m_useTemp | m_usePerlin) {
      if (m_useTemp) ComputeMaxTemperature();
      Spectrum sigma_t_2 = (sigma_a_2 + sigma_s_2 + sigma_n_2) * scale;
      for (int i = 0; i < Spectrum::nSamples; i++) {
        m_sigma_t[i] = std::max(m_sigma_t[i], sigma_t_2[i]);
      }
    }
  }

  bool OpenvdbGridMedium::open(const std::string& filename)
  {
    m_file = new openvdb::io::File(filename);
    if (m_file)
    {
      m_file->open();
      return true;
    }
    else
      return false;
  }

  void OpenvdbGridMedium::UpdateAndComputeMaxDensity() {
    if (!m_density_grid) return;
    
    float max = 0.0f;
    for (MonoGrid::ValueOnIter iter = m_density_grid->tree().beginValueOn();
        iter.test(); ++iter) {
      iter.setValue(std::pow(iter.getValue(), m_gamma));
      max = std::max(max, iter.getValue());
    }
    m_maxDensity = Spectrum(max);

    std::cout << "Max density calculated: " << m_maxDensity << std::endl;
  }

  void OpenvdbGridMedium::ComputeMaxTemperature() {
    if (!m_temperature_grid) return;
    m_maxTemperature = 0.0f;
    float max = 0.0f;
    for (MonoGrid::ValueOnCIter iter = m_temperature_grid->tree().cbeginValueOn();
        iter.test(); ++iter) {
      max = std::max(max, *iter);
    }
    m_maxTemperature = max;
    std::cout << "Max temperature calculated: " << m_maxTemperature << std::endl;
  }

  Spectrum OpenvdbGridMedium::Density(const Point3f& p) const
  {
    lookups++;
    openvdb::Vec3f pInd(p.x * m_n.x + m_min.x, 
                        p.y * m_n.y + m_min.y, 
                        p.z * m_n.z + m_min.z);
    MonoGridSampler sampler(m_density_grid->tree(), 
                            m_density_grid->transform());

    double density = sampler.isSample(pInd);

    if (m_useLinearRamp) {
      double minThreshold = 0.01f;
      if (density <= 0.0f) return Spectrum(density);
      double updatedDensity = density + m_rampScale * (p[m_rampDimension] - m_rampOrigin);
      return Spectrum(std::max(std::min(updatedDensity, Float(m_maxDensity[0])), minThreshold));
    }
    return Spectrum(density);
  }

  Float OpenvdbGridMedium::Temperature(const Point3f& p) const
  {
    openvdb::Vec3f pInd(p.x * m_temp_n.x + m_temp_min.x, 
                        p.y * m_temp_n.y + m_temp_min.y, 
                        p.z * m_temp_n.z + m_temp_min.z);
    MonoGridSampler sampler(m_temperature_grid->tree(), 
                            m_temperature_grid->transform());
    return sampler.isSample(pInd);
  }

  Point3i OpenvdbGridMedium::getDims() const { return m_n; }

  Spectrum OpenvdbGridMedium::Sample(const Ray &ray, Sampler &sampler,   
         MemoryArena &arena, MediumInteraction *mi) const {
    std::cout << "Sample not implemented for openvdb grid." << std::endl;
    return Spectrum(0.0f);
  }
  
  Spectrum OpenvdbGridMedium::SampleChannel(const Ray &rWorld, Sampler &sampler, 
         MemoryArena &arena, MediumInteraction *mi, int channel) const {
    ProfilePhase _(Prof::MediumSample);
    Ray normRay = Ray(rWorld.o, Normalize(rWorld.d), rWorld.d.Length() * rWorld.tMax);
    Ray ray = m_worldToMedium(normRay);
    const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
    Float tMin, tMax;
    if (!b.IntersectP(ray, &tMin, &tMax)) return Spectrum(1.f);
    Float t = tMin - std::log(1 - sampler.Get1D()) / GetMajorant()[channel];
    bool sampledMedium = t < tMax;
    t = std::min(t, tMax);
    if (sampledMedium)
      *mi = MediumInteraction(normRay(t), -normRay.d, normRay.time, 
                              this, ARENA_ALLOC(arena, HenyeyGreenstein)(m_g));
    return Exp(-GetMajorant() * (t-tMin));
  }

  Spectrum OpenvdbGridMedium::Tr(const Ray& ray, Sampler& sampler) const {
    std::cerr << "Tr not implemented for openvb grid" << std::endl;
    return Spectrum(0.0f);
  }  

  Spectrum OpenvdbGridMedium::GetAbsorption(const Point3f &p) const {
    if (m_useTemp) {
      Spectrum lerped_sigma_a = InterpolateWithTemperature(p, m_sigma_a, 
                                                              m_sigma_a_2);
      return Density(m_worldToMedium(p)) * lerped_sigma_a;
    } else if (m_usePerlin) {
      Spectrum lerped_sigma_a = InterpolateWithPerlin(p, m_sigma_a, m_sigma_a_2);
      return Density(m_worldToMedium(p)) * lerped_sigma_a;
    } else {
      return Density(m_worldToMedium(p)) * m_sigma_a;
    }
  };

  Spectrum OpenvdbGridMedium::GetScattering(const Point3f &p) const {
    if (m_useTemp) {
      Spectrum lerped_sigma_s = InterpolateWithTemperature(p, m_sigma_s, 
                                                              m_sigma_s_2);
      return Density(m_worldToMedium(p)) * lerped_sigma_s;
    } else if (m_usePerlin) {
      Spectrum lerped_sigma_s = InterpolateWithPerlin(p, m_sigma_s, m_sigma_s_2);
      return Density(m_worldToMedium(p)) * lerped_sigma_s;
    } else {
      return Density(m_worldToMedium(p)) * m_sigma_s;
    }
  };

  void OpenvdbGridMedium::GetCoefficients(const Point3f &p, Spectrum &s, Spectrum &a, Spectrum &n) const {
    Spectrum lerped_s, lerped_a;
    if (m_useTemp) {
      lerped_s = InterpolateWithTemperature(p, m_sigma_s, m_sigma_s_2);
      lerped_a = InterpolateWithTemperature(p, m_sigma_a, m_sigma_a_2);
    } else if (m_usePerlin) {
      lerped_s = InterpolateWithPerlin(p, m_sigma_s, m_sigma_s_2);
      lerped_a = InterpolateWithPerlin(p, m_sigma_a, m_sigma_a_2);
    } else {
      lerped_s = m_sigma_s;
      lerped_a = m_sigma_a;
    }
    Spectrum density = Density(m_worldToMedium(p)); 
    s = density * lerped_s;
    a = density * lerped_a;
    n = GetMajorant() - (s + a);
  };

  Spectrum OpenvdbGridMedium::GetMajorant() const {
    return m_maxDensity * m_sigma_t;
  };

  bool OpenvdbGridMedium::SetMinAndMaxPosition(Ray &ray) const {
      Ray normRay = Ray(ray.o, Normalize(ray.d), ray.d.Length() * ray.tMax);
      Ray r = m_worldToMedium(normRay);
      const Bounds3f b(Point3f(0, 0, 0), Point3f(1, 1, 1));
      Float tMin, tMax;
      if (!b.IntersectP(r, &tMin, &tMax)) return false;
      ray.o = normRay(tMin);
      ray.tMax = tMax; 
      return true;
  };
  
  Spectrum OpenvdbGridMedium::InterpolateWithTemperature(Point3f p, 
                                                         Spectrum x, 
                                                         Spectrum y) const {
    Float temp = Temperature(m_worldToMedium(p));
    Float m = (m_maxTemperature - temp)  / m_maxTemperature;
    return m * x + (1.0f - m) * y;
  }

  Spectrum OpenvdbGridMedium::InterpolateWithPerlin(Point3f p, Spectrum x, Spectrum y) const {
    Float noise = Noise(p * m_noiseScale);
    Float m = (noise + 1.0f) / 2.0f;
    return m * x + (1.0f - m) * y;
  }
}
