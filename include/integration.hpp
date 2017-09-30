//          Copyright Federico Kircheis 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CALCNUM_INTEGRATION_HPP_5830907329
#define CALCNUM_INTEGRATION_HPP_5830907329

#include "intervals.hpp"
#include "utils.hpp"

#include "interpolation.hpp"

#include <functional>
#include <iterator>
#include <cassert>

namespace calcnum{

	inline double integrate_pto_medio(std::function<double(double)> f, const open_interval& interv, std::size_t subintervals = 1){
		assert(subintervals > 0);
		const double h = length(interv)/d(subintervals);

		double result = 0;
		for(std::size_t i = 1; i != subintervals+1; ++i){
			double x_i = std::fma(i, h, interv.a - h/2);
			result = std::fma(h, f(x_i), result);
		}
		return result;
	}

	inline double integrate_trapezoidal(std::function<double(double)> f, closed_interval interv, std::size_t subintervals = 1){
		assert(subintervals > 0);
		const double h = length(interv)/d(subintervals);

		double result = (f(interv.a) + f(interv.b))/2;
		for(std::size_t i = 1; i != subintervals; ++i){
			double x_i = std::fma(i, h, interv.a);
			result += f(x_i);
		}
		return h*result;
	}

	template<class intervaltype>
	double integrate_newton_cotes(std::function<double(double)> f, const intervaltype& interv, std::size_t degree, std::size_t subintervals = 1){
		// cal
		assert(degree > 0);
		assert(subintervals > 0);

		double par_sum = 0;
		equidistant_nodes_closed ext_nodes(closed_interval{interv.a, interv.b}, subintervals+1);
		for(auto it = std::begin(ext_nodes); it != std::end(ext_nodes)-1; ++it){
			intervaltype interv2 = {*it, *(it+1)};

			auto nodes = make_node_gen(interv2, degree);
			std::vector<double> data_x;
			data_x.reserve(degree);
			std::copy(std::begin(nodes), std::end(nodes), std::back_inserter(data_x));

			for(std::size_t i = 0; i != data_x.size(); ++i){
				auto weight = integrate(lagr_base_poly(std::begin(data_x), std::end(data_x), i), interv2);
				par_sum += f(data_x[i]) * weight;
			}

		}
		return par_sum;
	}

}

#endif
