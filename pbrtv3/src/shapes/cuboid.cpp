#include "shapes/cuboid.h"
#include "paramset.h"
#include "efloat.h"
#include "stats.h"

namespace pbrt {
	Cuboid::Cuboid(const Transform *o2w, const Transform *w2o, bool reverseOrientation)
		: Shape(o2w, w2o, reverseOrientation)
	{
		m_box.pMin = (*o2w)(Point3f(-1.0f));
		m_box.pMax = (*o2w)(Point3f(1.0f));
	}

	Bounds3f Cuboid::ObjectBound() const
	{
		return m_box;
	}

	bool Cuboid::Intersect(const Ray & ray, Float * tHit, SurfaceInteraction * isect, bool testAlphaTexture) const
	{
		return false;
	}

	bool Cuboid::IntersectP(const Ray & ray, bool testAlphaTexture) const
	{
		return m_box.IntersectP(ray);
	}

	Float Cuboid::Area() const
	{
		return m_box.SurfaceArea();
	}

	Interaction Cuboid::Sample(const Point2f & u, Float * pdf) const
	{
		return Interaction();
	}

	std::shared_ptr<Cuboid> CreateCuboidShape(const Transform *o2w,
		const Transform *w2o,
		bool reverseOrientation,
		const ParamSet &params) {
		return std::make_shared<Cuboid>(o2w, w2o, reverseOrientation);
	}
	
}