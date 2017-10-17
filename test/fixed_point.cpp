//          Copyright Federico Kircheis 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <catch.hpp>


#include "fixed_point.hpp"
#include "statistics.hpp"

#include <functional>


namespace{
	double cos_to_fix(double x){
		return std::cos(x);
	}
	double sin_to_fix(double x){
		return std::sin(x);
	}

	std::vector<double> simple_fix_point_test(const std::function<double(double)>& f, double x0, double sol, std::size_t max_iter = 100){
		std::vector<double> err;
		calcnum::fixed_point it(f, x0);
		while(std::fabs(it->f - sol) > 0.0000000001 && max_iter > 0){
			--max_iter;
			++it;
			err.push_back(std::abs(sol-it->f));
		}
		return err;
	}
}
TEST_CASE("fixed_point", "[fixed_point]"){
	SECTION("cos(x)=x"){
		auto err = simple_fix_point_test(cos_to_fix, 1, 0.73908513321516);
		const auto conv = calcnum::calculate_convergency(err);
		auto res = calcnum::analyze_data(calcnum::clear_from_inf_nan(conv));
		INFO(res);
		REQUIRE_FALSE(calcnum::is_outlier(res, 1));
	}
	SECTION("sin(x)=x"){
		auto err = simple_fix_point_test(sin_to_fix, 1, 0);
		const auto conv = calcnum::calculate_convergency(err);
		auto res = calcnum::analyze_data(calcnum::clear_from_inf_nan(conv));
		INFO(res);
		REQUIRE_FALSE(calcnum::is_outlier(res, 1));
	}
}
