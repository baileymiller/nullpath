#include "transmittance_exponential.h"
#include "spectrum.h"

namespace pbrt
{
	TransmittanceExponential::TransmittanceExponential()
	{
	}
	
	Spectrum TransmittanceExponential::TrMedium(const Float opticalDepth) const
	{
		return Spectrum(expf(-opticalDepth));
	}

	Spectrum TransmittanceExponential::TrSurface(const Float opticalDepth) const
	{
		return Spectrum(expf(-opticalDepth));
	}

	// the pdf is \sigma_t * exp(-\sigmaT * distance)
	Spectrum TransmittanceExponential::mediumPDF(const Float opticalDepth) const
	{
        return TrMedium(opticalDepth);
	}

	// the tracklength pdf is exp(-\sigma_t * distance)
	Spectrum TransmittanceExponential::surfacePDF(const Float opticalDepth) const
	{
        return TrSurface(opticalDepth);
	}

	Float TransmittanceExponential::sample(const Float u0) const
	{
		return -std::log(1.0f - u0);
	}

    Float TransmittanceExponential::sample2(const Float u0) const
    {
        return -std::log(1.0f - u0);
    }
    Spectrum TransmittanceExponential::sample2PDF(const Float tau) const
    {
        return TrSurface(tau);
    }
    Spectrum TransmittanceExponential::sample2Prob(const Float tau) const
    {
        return TrSurface(tau);
    }
}