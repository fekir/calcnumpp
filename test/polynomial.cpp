//          Copyright Federico Kircheis 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <catch.hpp>

#include "utils.hpp"
#include "intervals.hpp"
#include "polynomial.hpp"

#include <functional>
#include <iterator>

using namespace calcnum;

TEST_CASE("polynomial"){

	SECTION("normalize"){
		const polynomial p({0, 0, 1, 0});
		REQUIRE(p.coeffs().size()==3);
	}

	const polynomial poly({-4, 14, -7, 1});

	SECTION("evaluate-hoerner"){
		auto res = poly(0);
		REQUIRE(res == Approx(-4));

		res = poly(1);
		REQUIRE(res == Approx(4));
		res = poly(2);
		REQUIRE(res == Approx(4));
		res = poly(4);
		REQUIRE(res == Approx(4));
		res = poly(0.341032918083006);
		REQUIRE(res == Approx(0));
	}

	const polynomial poly2({1, 2, 3});

	SECTION("sum"){
		const auto num_coeff = std::max(poly2.coeffs().size(), poly.coeffs().size());
		const auto poly3 = poly + poly2;
		REQUIRE(poly3.coeffs().size() <= num_coeff);
		REQUIRE(poly3.coeffs()[0] == Approx(-3));
		REQUIRE(poly3.coeffs()[1] == Approx(16));
		REQUIRE(poly3.coeffs()[2] == Approx(-4));
		REQUIRE(poly3.coeffs()[3] == Approx(1));
	}

	SECTION("multiply-scalar"){
		const auto scalar = 2;
		const auto num_coeff = poly2.coeffs().size() + poly.coeffs().size()-1;
		const auto poly3 = poly * scalar;
		REQUIRE(poly3.coeffs().size() == poly.coeffs().size());
		for(std::size_t i = 0; i != poly.coeffs().size();++i ){
			REQUIRE(scalar*poly.coeffs()[i] == Approx(poly3.coeffs()[i]));
		}
	}

	SECTION("multiply"){
		const auto num_coeff = poly2.coeffs().size() + poly.coeffs().size()-1;
		const auto poly3 = poly * poly2;
		REQUIRE(poly3.coeffs().size() == num_coeff);
		REQUIRE(poly3.coeffs()[0] == Approx(-4));
		REQUIRE(poly3.coeffs()[1] == Approx(6));
		REQUIRE(poly3.coeffs()[2] == Approx(9));
		REQUIRE(poly3.coeffs()[3] == Approx(29));
		REQUIRE(poly3.coeffs()[4] == Approx(-19));
		REQUIRE(poly3.coeffs()[5] == Approx(3));
	}

}


TEST_CASE("derive"){
	const polynomial poly({-4, 14, -7, 1});
	const polynomial dpoly = derive(poly);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
	REQUIRE(dpoly.coeffs() == std::vector<double>({14, -7.0 * 2, 1.0 *3}));
#pragma GCC diagnostic pop


	SECTION("integrate"){
		const auto ipoly = integrate(poly);
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
		REQUIRE(derive(ipoly).coeffs() == poly.coeffs());
#pragma GCC diagnostic pop


		SECTION("integrate on interval"){
			polynomial poly2({1,1});
			auto integral = integrate(poly2, closed_interval{0,4});
			REQUIRE(integral == Approx(12));
		}
	}
}
