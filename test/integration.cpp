//          Copyright Federico Kircheis 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <catch.hpp>

#include "integration.hpp"
#include "utils.hpp"
#include "statistics.hpp"

#include <functional>

namespace  {

	double fun_to_integrate(double x){
		return x - std::sin(x);
	}
	const calcnum::closed_interval close_int_to_integrate{0,10};
	const calcnum::open_interval   open_int_to_integrate{0,10};
	const double integrated_value = 48.160928470923547547741136052175935166;

}
TEST_CASE("integrate_pto_medio"){
	std::vector<double> err;
	std::vector<double> step;
	for(std::size_t i = 1; i != 10; ++i){
		auto val = calcnum::integrate_pto_medio(fun_to_integrate, open_int_to_integrate, i * 100);
		auto new_err = std::fabs(val-integrated_value);
		err.push_back(new_err);
		step.push_back(length(open_int_to_integrate)/(i*100.0));
	}

	const auto conv = calcnum::calculate_convergency(err, step);
	auto res = calcnum::analyze_data(calcnum::clear_from_inf_nan(conv));
	REQUIRE_FALSE(calcnum::is_outlier(res, 2));

}

TEST_CASE("integrate_trapezi"){
	std::vector<double> err;
	std::vector<double> step;
	for(std::size_t i = 1; i != 10; ++i){
		auto val = calcnum::integrate_trapezoidal(fun_to_integrate, close_int_to_integrate, i * 100);
		auto new_err = std::fabs(val-integrated_value);
		err.push_back(new_err);
		step.push_back(length(open_int_to_integrate)/(i*100.0));
	}
	const auto conv = calcnum::calculate_convergency(err, step);
	auto res = calcnum::analyze_data(calcnum::clear_from_inf_nan(conv));
	INFO(res);
	REQUIRE_FALSE(calcnum::is_outlier(res, 2));

}

TEST_CASE("newton-cotes"){
	using calcnum::approx_equal;
	SECTION("pto_medio"){
		auto val_medio = integrate_pto_medio(fun_to_integrate, open_int_to_integrate);
		auto val_cotes = integrate_newton_cotes(fun_to_integrate, open_int_to_integrate, 1);
		REQUIRE(val_cotes == Approx(val_medio));
	}
	SECTION("trapezi"){
		auto val_trapezi = integrate_trapezoidal(fun_to_integrate, close_int_to_integrate);
		auto val_cotes = integrate_newton_cotes(fun_to_integrate, close_int_to_integrate, 2);
		REQUIRE(val_cotes == Approx(val_trapezi));
	}
	SECTION("cav_simps"){
		auto val_cotes = integrate_newton_cotes(fun_to_integrate, close_int_to_integrate, 3);
		auto val_cavsimps= 57.299530349;
		REQUIRE(val_cotes == Approx(val_cavsimps));
	}
}

TEST_CASE("newton-cotes composito"){
	using calcnum::approx_equal;
	SECTION("pto_medio"){
		for(std::size_t i = 1; i != 100; ++i){
			auto val_medio = integrate_pto_medio(fun_to_integrate, open_int_to_integrate, i);
			auto val_cotes = integrate_newton_cotes(fun_to_integrate, open_int_to_integrate, 1, i);
			REQUIRE(val_cotes == Approx(val_medio));
		}
	}
	SECTION("trapezi"){
		for(std::size_t i = 1; i != 100; ++i){
			auto val_trapezi = integrate_trapezoidal(fun_to_integrate, close_int_to_integrate, i);
			auto val_cotes = integrate_newton_cotes(fun_to_integrate, close_int_to_integrate, 2, i);
			REQUIRE(val_cotes == Approx(val_trapezi));
		}
	}
	SECTION("cav_simps"){
		auto val_cotes = integrate_newton_cotes(fun_to_integrate, close_int_to_integrate, 3,1);
		auto val_cavsimps= 57.299530349;
		REQUIRE(val_cotes == Approx(val_cavsimps));
	}
}
