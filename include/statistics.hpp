#ifndef CALCNUM_STATISTICS_HPP_0131247378
#define CALCNUM_STATISTICS_HPP_0131247378

#include "utils.hpp"

#include <cassert>
#include <cmath>

#include <vector>
#include <algorithm>
#include <numeric>

namespace calcnum{

	struct sd_mean{
		double standard_deviation;
		double mean;
	};

	inline sd_mean analyze_data(std::vector<double> data){
		const double sum = std::accumulate(data.begin(), data.end(), 0.0); // FIXME:
		const double mean = sum/data.size();
		double variance = 0;
		for(const auto d : data){
			variance += (d - mean)*(d - mean);
		}
		variance /= data.size();
		const double standard_deviation = std::sqrt(variance);
		return {standard_deviation, mean};
	}

	inline bool is_outlier(const sd_mean& info, double d){
		const double  min = (info.mean - info.standard_deviation*3);
		const double  max = (info.mean + info.standard_deviation*3);
		return d < min || d > max;
	}

}

#endif
