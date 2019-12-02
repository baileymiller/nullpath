#ifndef __MEDIA_COMMON_H__
#define __MEDIA_COMMON_H__

#include "medium.h"
#include "transform.h"
#include "stats.h"

// common functions for different media types
// adding this file so as to keep the interface clean

namespace pbrt
{
	// Explicit tabulated transmittance function for sampling and evaluating.
	enum Normalization { NO_NORMALIZATION, MEAN_NORMALIZATION, MEDIAN_NORMALIZATION, UNIT_MFP_NORMALIZATION };

	// Different sampling methods for different types of methods
	enum SamplingMethod
	{
        SAMPLE_EXP,									// Sample an exponential media using a curve
        SAMPLE_EXP_WOODCOCK,                        // Sample an exponential media using regular woodcock tracking
		SAMPLE_NONEXP_FBM_RAYMARCH,					// Jan's method - sample an optical depth from an underlying curve and accumulate density till you reach that
		SAMPLE_AVERAGE_TRFUNCTION,					// sample average tr function on a per face basis
		SAMPLE_GAMMA_NONEXP_DAVIS_RAYMARCH,			// Jan's method - ray marched solution - use Davis equations
		SAMPLE_GAMMA_NONEXP_DAVIS_WOODCOCK,			// Jan's method - woodcock routines - use Davis equations
		SAMPLE_GAMMA_NONEXP_DAVIS_REGULARTRACKING,	// Jan's method - regular tracking - use Davis equations
		SAMPLE_GAMMA_PRODUCT_NAIVE_MAJORANT,		// woodcock tracking uses the \sigma_t scaling factor as the constant majorant and exponential steps
		SAMPLE_GAMMA_PRODUCT_TIGHT_MAJORANT,		// woodcock tracking uses a local majorant that is a product of implict function and constant majorant
		SAMPLE_GAMMA_PRODUCT_GAMMA_MAJORANT,		// woodcock tracking uses gamma steps and gammam function as the majorant.
		SAMPLE_GAMMA_PRODUCT_RAYMARCH,				// ray march to accumulate density provided by exp density
		SAMPLE_UNKNOWN_METHOD						// unknow method
	};

    // Different MIP methods
    enum MipMethod
    {
        MIP_NAIVE,                                  // use naive downsampling of sigmaT values
        MIP_EFFECTIVE,                              // use effective downsampling of sigmaT values
		MIP_EFFECTIVE_HISTOGRAM,					// use effective sigmaT by employing our effective histogram method
		MIP_AVERAGE_TR								// average TR values computed at various mip resolutions
    };

	class TabulatedTransmittanceFunction
	{
	public:
		TabulatedTransmittanceFunction(Normalization normalization = NO_NORMALIZATION);
		void insert(Float sampledDistance);
		Float sample(Float epsilon) const;
		void prepare();
		Float getDistance(size_t index) const { if (index < m_distances.size()) return m_distances[index]; else return 0.0f; }
		Float getVals(size_t index) const { if (index < m_vals.size()) return m_vals[index]; else return 0.0f; }
		Float getTransmittance(Float accumulatedOpticalThickness) const;
		size_t getSize() const { return m_distances.size(); }
		void writeDistances(const std::string& filename) const;
	private:
		void normalizeMeanDistances();
		void normalizeMedianValue();

		std::vector<Float> m_distances;
		std::vector<Float> m_vals;

		std::vector<Float> m_distances1;	// used for transmittance computation
		std::vector<Float> m_vals1;			// used for transmittance computation
		Normalization m_normalization;
	};

	// A convenient class to compute running averages, stds and variances
	// code from: https://www.johndcook.com/blog/standard_deviation/
	class RunningStat
	{
	public:
		RunningStat() : m_n(0) {}

		void clear()
		{
			m_n = 0;
		}

		void push(Float x)
		{
			m_n++;

			// See Knuth TAOCP vol 2, 3rd edition, page 232
			if (m_n == 1)
			{
				m_oldM = m_newM = x;
				m_oldS = 0.0;
			}
			else
			{
				m_newM = m_oldM + (x - m_oldM) / m_n;
				m_newS = m_oldS + (x - m_oldM)*(x - m_newM);

				// set up for next iteration
				m_oldM = m_newM;
				m_oldS = m_newS;
			}
		}

		int numDataValues() const
		{
			return m_n;
		}

		Float mean() const
		{
			return (m_n > 0) ? m_newM : 0.0;
		}

		Float variance() const
		{
			return ((m_n > 1) ? m_newS / (m_n - 1) : 0.0);
		}

		Float standardDeviation() const
		{
			return sqrt(variance());
		}

	private:
		int m_n;
		Float m_oldM, m_newM, m_oldS, m_newS;
	};

	// given a vector, gives an index into a probable table of values
	// based on the maximum compnnent of the ray-direction
	int maxDirectionComponent(const Vector3f& v);

	// some utility methods for printing debugging messages
	std::string getNormalizationString(Normalization n);
	std::string getSamplingMethodString(SamplingMethod s);
	std::string getMipMethodString(MipMethod m);
	std::string getTransmittanceModeString(TransmittanceMode m);
	std::string getTransmittanceFunctionTypeString(TransmittanceFunctionType t);
}

#endif


