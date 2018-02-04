//          Copyright Federico Kircheis 2017-2018
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CALCNUM_FUN_HPP_0790906003
#define CALCNUM_FUN_HPP_0790906003

#include "matrix.hpp"

#include <functional>

namespace calcnum{
	// conversion between functions with 1D-vector and scalar as parameter
	inline std::function<calcnum::vector<1>(double, calcnum::vector<1>)> vectorize(std::function<double(double,double)> f){
		return [f](double t, calcnum::vector<1> v){
			auto res = f(t, v[0]);
			return calcnum::vector<1>({res});
		};
	}
	inline std::function<calcnum::vector<1>(double, calcnum::vector<1>)> vectorize(std::function<calcnum::vector<1>(double, calcnum::vector<1>)> f){
		return f;
	}

	inline std::function<double(double,double)> linearize(std::function<calcnum::vector<1>(double, calcnum::vector<1>)> f){
		return [f](double t, double v){
			auto res = f(t, calcnum::vector<1>({v}));
			return res[0];
		};
	}
	inline std::function<double(double,double)> linearize(std::function<double(double,double)> f){
		return f;
	}

	using cauchy_time = double;
	using cauchy_var = double;
	using cauchy_eq = std::function<cauchy_var(cauchy_time,cauchy_var)>;
	struct cauchy_prob{
		cauchy_eq f;
		cauchy_time t0;
		cauchy_var y0;
	};

}

#endif
