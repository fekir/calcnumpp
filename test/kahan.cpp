//          Copyright Federico Kircheis 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <catch.hpp>

#include "utils.hpp"


using namespace calcnum;


TEST_CASE("kahan", "[utils][kahan]"){
	SECTION("verify plus equivalence"){
		const auto err = 0.0001;
		kahan_sum h1;
		kahan_sum h2;
		kahan_sum h3;
		for(int i =0; i != 10; ++i){
			h1 = h1 + 1;
			h2 += 1;
			++h3;
		}
		REQUIRE( approx_equal(d(h1), 10, err));
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
		REQUIRE( d(h1) == d(h2));
		REQUIRE( d(h1) == d(h3));
#pragma GCC diagnostic pop
	}
	SECTION("verify minus equivalence"){
		const auto err = 0.0001;
		kahan_sum h1;
		kahan_sum h2;
		kahan_sum h3;
		for(int i = 0; i != 10; ++i){
			h1 = h1 - 1;
			h2 -= 1;
			--h3;
		}
		REQUIRE( approx_equal(d(h1), -10, err));
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
		REQUIRE( d(h1) == d(h2));
		REQUIRE( d(h1) == d(h3));
#pragma GCC diagnostic pop
	}

	// FIXME: find an example that shows that the sum is more precise/stable
}
