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


	struct variance_accumulator{
		double mean;
		kahan_sum operator()(const kahan_sum& a, double b){
			return a + (b - mean)*(b - mean);
		}
	};

	inline sd_mean analyze_data(std::vector<double> data){
		assert(data.size() > 0 && "needs to be > 0");
		std::sort(data.begin(), data.end());
		const auto sum_k = std::accumulate(data.begin(), data.end(), kahan_sum());
		auto sum = d(sum_k);
		const double mean = sum/d(data.size());

		variance_accumulator var_acc = {mean};
		const auto variance_k = std::accumulate(data.begin(), data.end(), kahan_sum(), var_acc);
		assert(d(variance_k) >= 0 && "needs to be >= 0");
		auto variance = d(variance_k)/d(data.size());
		const double standard_deviation = std::sqrt(variance);
		return {standard_deviation, mean};
	}

	inline bool is_outlier(const sd_mean& info, double d){
		const double min = (info.mean - info.standard_deviation*1.5);
		const double max = (info.mean + info.standard_deviation*1.5);
		return d < min || d > max;
	}

}

#endif
