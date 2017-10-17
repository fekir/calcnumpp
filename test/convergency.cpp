//          Copyright Federico Kircheis 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <catch.hpp>

#include "utils.hpp"

TEST_CASE("convergency (newton)"){
	using calcnum::approx_equal;

	std::vector<double> err = {0.1};
	auto multby = 2; // generate randomly at every step, but != err[0] (?)
	SECTION("linear"){
		for(int i = 0; i != 10; ++i){
			err.push_back(multby*err.back());
		}
		auto conv = calcnum::calculate_convergency(err);

		for(auto v : conv){
			REQUIRE(v == Approx(1));
		}
	}
	SECTION("quadratic"){
		for(int i = 0; i != 5; ++i){
			err.push_back(multby*(err.back()*err.back()));
		}
		auto conv = calcnum::calculate_convergency(err);

		for(auto v : conv){
			REQUIRE(v == Approx(2));
		}
	}
}

TEST_CASE("convergency (integration)"){
	using calcnum::approx_equal;

	std::vector<double> err;
	std::vector<double> step = {0.1};
	for(int i = 0; i != 10; ++i){
		step.push_back(0.1*step.back());
	}
	auto multby = 2; // generate randomly at every step
	SECTION("linear"){
		for(auto s : step){
			err.push_back(multby*s);
		}
		auto conv = calcnum::calculate_convergency(err, step);

		for(auto v : conv){
			REQUIRE(v == Approx(1));
		}
	}
	SECTION("quadratic"){
		for(auto s : step){
			err.push_back(multby*s*s);
		}
		auto conv = calcnum::calculate_convergency(err, step);

		for(auto v : conv){
			REQUIRE(v == Approx(2));
		}
	}
}
