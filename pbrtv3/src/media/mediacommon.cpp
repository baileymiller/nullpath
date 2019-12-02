#include "mediacommon.h"
#include <fstream>

namespace pbrt
{
	// Tabulated Transmittance Function 
	TabulatedTransmittanceFunction::TabulatedTransmittanceFunction(Normalization normalization)
		: m_normalization(normalization)
	{
		m_distances.push_back(0.0);
	}

	void TabulatedTransmittanceFunction::insert(Float sampledDistance)
	{
		m_distances.push_back(sampledDistance);
	}

	Float TabulatedTransmittanceFunction::sample(Float epsilon) const
	{
		// find the value in the sorted function
		auto l = std::upper_bound(m_vals.begin(), m_vals.end(), epsilon);
		auto index = l - m_vals.begin();

		// get the next and prev locations within the function list
		auto rangeCurr = m_distances[index];
		auto rangePrev = m_distances[index - 1];
		auto funCurr = m_vals[index];
		auto funPrev = m_vals[index - 1];

		auto lerpFactor = abs((funPrev - epsilon) / (funCurr - funPrev));
		return (1.0f - lerpFactor) * rangePrev + lerpFactor * rangeCurr;
	}

	void TabulatedTransmittanceFunction::prepare()
	{
		// Sort the distances
		m_distances.push_back(std::numeric_limits<float>::max());
		std::sort(m_distances.begin(), m_distances.end());

		// store the sorted distances
		m_distances1 = m_distances;

		// compute the transmittance function
		auto total = m_distances.size();
		auto current = total;
		Float invTotal = 1.0f / Float(total);
		for (auto& d : m_distances)
		{
			Float ratio = Float(current) / Float(total);
			m_vals.push_back(ratio);
			current--;
		}

		// change the last value to zero
		auto location = m_vals.size();
		m_vals[location - 1] = 0.0f;

		// keep the vals the same way they are
		m_vals1 = m_vals;

		// sort the vals
		std::sort(m_vals.begin(), m_vals.end());
		std::reverse(m_distances.begin(), m_distances.end());

		switch (m_normalization)
		{
		case pbrt::NO_NORMALIZATION:
			std::cout << "No Normalization performed" << std::endl;
			break;
		case pbrt::MEAN_NORMALIZATION:
			normalizeMeanDistances();
			std::cout << "Performing mean normalization" << std::endl;
			break;
		case pbrt::MEDIAN_NORMALIZATION:
			normalizeMedianValue();
			std::cout << "Performing median normalization" << std::endl;
			break;
		default:
			break;
		}
	}

	Float TabulatedTransmittanceFunction::getTransmittance(Float accumulatedOpticalThickness) const
	{
		// find the value in the sorted function
		auto l = std::upper_bound(m_distances1.begin(), m_distances1.end(), accumulatedOpticalThickness);
		auto index = l - m_distances1.begin();

		// get the next and prev locations within the function list
		auto rangeCurr = m_distances1[index];
		auto rangePrev = m_distances1[index - 1];
		auto funCurr = m_vals1[index];
		auto funPrev = m_vals1[index - 1];

		auto lerpFactor = abs((rangePrev - accumulatedOpticalThickness) / (rangePrev - rangeCurr));
		return (1.0f - lerpFactor) * funPrev + lerpFactor * funCurr;
	}

	void TabulatedTransmittanceFunction::normalizeMeanDistances()
	{
		Float totalDistances = 0.0;
		for (auto& d : m_distances)
		{
			if (d == std::numeric_limits<float>::max()) continue;
			totalDistances += d;
		}

		Float meanDistance = totalDistances / Float(m_distances.size());

		for (auto& d : m_distances)
		{
			if (d == std::numeric_limits<float>::max()) continue;
			d /= meanDistance;
		}

		for (auto& d : m_distances1)
		{
			if (d == std::numeric_limits<float>::max()) continue;
			d /= meanDistance;
		}
	}

	void TabulatedTransmittanceFunction::normalizeMedianValue()
	{
		Float medianValueDistance = sample(0.5);
		const Float expMedianValueDistance = 0.69314718056f;
		const Float factor = (1.0 / medianValueDistance) * expMedianValueDistance;
		for (auto& d : m_distances)
		{
			if (d == std::numeric_limits<float>::max()) continue;
			d *= factor;
		}
		for (auto& d : m_distances1)
		{
			if (d == std::numeric_limits<float>::max()) continue;
			d *= factor;
		}
	}

	void TabulatedTransmittanceFunction::writeDistances(const std::string& filename) const
	{
		std::cout << "Writing out the tabulated function : " << filename << std::endl;
		std::ofstream out(filename.c_str());
		for (auto& d : m_distances1)
		{
			if (d == std::numeric_limits<float>::max()) continue;
			out << d << ",";
		}
		out.close();
	}

	int maxDirectionComponent(const Vector3f & v)
	{
		if (abs(v.x) >= abs(v.y) && abs(v.x) >= abs(v.z))
		{
			if (v.x >= 0.0) return 0;
			else return 1;
		}
		else if (abs(v.y) >= abs(v.x) && abs(v.y) >= abs(v.z))
		{
			if (v.y >= 0.0) return 2;
			else return 3;
		}
		else if (abs(v.z) >= abs(v.x) && abs(v.z) >= abs(v.y))
		{
			if (v.z >= 0.0) return 4;
			else return 5;
		}
		else
		{
			std::cout << "Something wrong : v : [ " << v.x << "," << v.y << "," << v.z << "]" << std::endl;
			return 0;
		}
	}

	std::string getNormalizationString(Normalization n)
	{
		switch (n)
		{
		case pbrt::NO_NORMALIZATION:
			return "NO_NORMALIZATION";
			break;
		case pbrt::MEAN_NORMALIZATION:
			return "MEAN_NORMALIZATION";
			break;
		case pbrt::MEDIAN_NORMALIZATION:
			return "MEDIAN_NORMALIZATION";
			break;
		case pbrt::UNIT_MFP_NORMALIZATION:
			return "UNIT_MFP_NORMALIZATION";
			break;
		default:
			return "UNKNOWN";
			break;
		}
	}
	
	std::string getSamplingMethodString(SamplingMethod s)
	{
		switch (s)
		{
		case pbrt::SAMPLE_EXP:
			return "SAMPLE_EXP";
			break;
		case pbrt::SAMPLE_EXP_WOODCOCK:
			return "SAMPLE_EXP_WOODCOCK";
			break;
		case pbrt::SAMPLE_NONEXP_FBM_RAYMARCH:
			return "SAMPLE_NONEXP_FBM_RAYMARCH";
			break;
		case pbrt::SAMPLE_AVERAGE_TRFUNCTION:
			return "SAMPLE_AVERAGE_TRFUNCTION";
			break;
		case pbrt::SAMPLE_GAMMA_NONEXP_DAVIS_RAYMARCH:
			return "SAMPLE_GAMMA_NONEXP_DAVIS_RAYMARCH";
			break;
		case pbrt::SAMPLE_GAMMA_NONEXP_DAVIS_WOODCOCK:
			return "SAMPLE_GAMMA_NONEXP_DAVIS_WOODCOCK";
			break;
		case pbrt::SAMPLE_GAMMA_NONEXP_DAVIS_REGULARTRACKING:
			return "SAMPLE_GAMMA_NONEXP_DAVIS_REGULARTRACKING";
			break;
		case pbrt::SAMPLE_GAMMA_PRODUCT_NAIVE_MAJORANT:
			return "SAMPLE_GAMMA_PRODUCT_NAIVE_MAJORANT";
			break;
		case pbrt::SAMPLE_GAMMA_PRODUCT_TIGHT_MAJORANT:
			return "SAMPLE_GAMMA_PRODUCT_TIGHT_MAJORANT";
			break;
		case pbrt::SAMPLE_GAMMA_PRODUCT_GAMMA_MAJORANT:
			return "SAMPLE_GAMMA_PRODUCT_GAMMA_MAJORANT";
			break;
		case pbrt::SAMPLE_GAMMA_PRODUCT_RAYMARCH:
			return "SAMPLE_GAMMA_PRODUCT_RAYMARCH";
			break;
		case pbrt::SAMPLE_UNKNOWN_METHOD:
			return "SAMPLE_UNKNOWN_METHOD";
			break;
		default:
			return "SAMPLE_UNKNOWN_METHOD";
			break;
		}
	}

	std::string getMipMethodString(MipMethod m)
	{
		switch (m)
		{
		case pbrt::MIP_NAIVE:
			return "MIP_NAIVE";
			break;
		case pbrt::MIP_EFFECTIVE:
			return "MIP_EFFECTIVE";
			break;
		case pbrt::MIP_EFFECTIVE_HISTOGRAM:
			return "MIP_EFFECTIVE_HISTOGRAM";
			break;
		case pbrt::MIP_AVERAGE_TR:
			return "MIP_AVERAGE_TR";
			break;
		default:
			return "UNKNOWN";
			break;
		}
	}
	
	std::string getTransmittanceModeString(TransmittanceMode m)
	{
		switch (m)
		{
		case pbrt::NONEXP_SSC_MSC_SMC_MMC:
			return "NONEXP_SSC_MSC_SMC_MMC";
			break;
		case pbrt::NONEXP_SST_MST_SMC_MMC:
			return "NONEXP_SST_MST_SMC_MMC";
			break;
		case pbrt::NONEXP_SST_MSC_SMT_MMC:
			return "NONEXP_SST_MSC_SMT_MMC";
			break;
		case pbrt::MODE_UNKNOWN:
			return "MODE_UNKNOWN";
			break;
		default:
			return "MODE_UNKNOWN";
			break;
		}
	}
	
	std::string getTransmittanceFunctionTypeString(TransmittanceFunctionType t)
	{
		switch (t)
		{
		case pbrt::TRANSMITTANCE_EXPONENTIAL:
			return "TRANSMITTANCE_EXPONENTIAL";
			break;
		case pbrt::TRANSMITTANCE_DAVIS:
			return "TRANSMITTANCE_DAVIS";
			break;
        case pbrt::TRANSMITTANCE_DAVIS_2PARAM:
            return "TRANSMITTANCE_DAVIS_2PARAM";
            break;
		case pbrt::TRANSMITTANCE_QUADRATIC_RAMP:
			return "TRANSMITTANCE_QUADRATIC_RAMP";
			break;
		case pbrt::TRANSMITTANCE_LINEAR_RAMP:
			return "TRANSMITTANCE_LINEAR_RAMP";
			break;
		case pbrt::TRANSMITTANCE_UNKNOWN:
			return "TRANSMITTANCE_UNKNOWN";
			break;
		default:
			return "TRANSMITTANCE_UNKNOWN";
			break;
		}
	}
}