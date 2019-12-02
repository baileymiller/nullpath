#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_TRANSMITTANCE_H
#define PBRT_CORE_TRANSMITTANCE_H

#include "pbrt.h"

namespace pbrt
{
	// Abstract base class transmittance
	class Transmittance
	{
	public:
		// Transmittance and corresponding pdf routines
		virtual Spectrum TrMedium(const Float opticalDepth) const = 0;
		virtual Spectrum TrSurface(const Float opticalDepth) const = 0;
		
        // technically in a path tracer, this is the probability, but we'll call it PDF for now.
        virtual Spectrum mediumPDF(const Float opticalDepth) const = 0;
        virtual Spectrum surfacePDF(const Float opticalDepth) const = 0;

		// sampling routine for a particular transmittance function
		virtual Float sample(const Float xi) const = 0;

        // Use a different sampling strategy based on Tr transmittance
        virtual Float sample2(const Float xi) const = 0;
        virtual Spectrum sample2PDF(const Float tau) const = 0;
        virtual Spectrum sample2Prob(const Float tau) const = 0;
	};
}

#endif
