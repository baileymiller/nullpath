#ifndef TRANSMITTANCE_EXPONENTIAL_H
#define TRANSMITTANCE_EXPONENTIAL_H

#include "transmittance.h"

namespace pbrt
{
	class TransmittanceExponential : public Transmittance
	{
	public:
		TransmittanceExponential();
		Spectrum TrMedium(const Float opticalDepth) const;
		Spectrum TrSurface(const Float opticalDepth) const;
		Spectrum mediumPDF(const Float opticalDepth) const;
		Spectrum surfacePDF(const Float opticalDepth) const;
		Float sample(const Float xi) const;
        Float sample2(const Float xi) const;
        Spectrum sample2PDF(const Float tau) const;
        Spectrum sample2Prob(const Float tau) const;
	};
}

#endif