//          Copyright Federico Kircheis 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <catch.hpp>

#include "bisec.hpp"


#include <functional>

namespace{
	// all test functions have a root in x = 5
	double fun_linear(double x){
		return x-5;
	}
	double fun_quad(double x){
		return x*x-25;
	}
	double fun_tri(double x){
		return x*x*x-125;
	}

	// simple test to see if the properties of the bisection algorithm are correct (intervals are at least halfed at every iteration)
	void simple_bisec_test(const std::function<double(double)>& f, calcnum::closed_interval interval){
		calcnum::bisect_iter it(f, interval);
		auto oldst = *it;
		while(length(it->interval)> 0.00001){
			++it;
			const auto& st = *it;
			REQUIRE(calcnum::length(oldst.interval)/2 >= calcnum::length(st.interval));
			oldst = st;
		}
	}
}

TEST_CASE("bisection", "[bisection]"){
	SECTION("normal case"){
		calcnum::closed_interval interval{0,11};
		SECTION("linear"){
			simple_bisec_test(fun_linear, interval);
		}
		SECTION("quad"){
			simple_bisec_test(fun_quad, interval);
		}
		SECTION("tri"){
			simple_bisec_test(fun_tri, interval);
		}
	}
	SECTION("exact root"){
		calcnum::closed_interval interval{0,10};
		SECTION("linear"){
			simple_bisec_test(fun_linear, interval);
		}
		SECTION("quad"){
			simple_bisec_test(fun_quad, interval);
		}
		SECTION("tri"){
			simple_bisec_test(fun_tri, interval);
		}
	}
}
