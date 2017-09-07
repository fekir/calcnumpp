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
	const auto err = 0.0001;
	const polynomial poly({-4, 14, -7, 1});
	const polynomial poly2({1, 2, 3});

	SECTION("normalize"){
		const polynomial p({0, 0, 1, 0});
		REQUIRE(p.coeffs().size()==3);
	}

	SECTION("evaluate-hoerner"){
		auto res = poly(0);
		REQUIRE(approx_equal(res, -4, err));

		res = poly(1);
		REQUIRE(approx_equal(res, 4, err));
		res = poly(2);
		REQUIRE(approx_equal(res, 4, err));
		res = poly(4);
		REQUIRE(approx_equal(res, 4, err));
		res = poly(0.34103);
		REQUIRE(approx_equal(res, 0, err));
	}



	SECTION("sum"){
		const auto num_coeff = std::max(poly2.coeffs().size(), poly.coeffs().size());
		const auto poly3 = poly + poly2;
		REQUIRE(poly3.coeffs().size() <= num_coeff);
		REQUIRE(approx_equal(poly3.coeffs()[0], -3, err));
		REQUIRE(approx_equal(poly3.coeffs()[1], 16, err));
		REQUIRE(approx_equal(poly3.coeffs()[2], -4, err));
		REQUIRE(approx_equal(poly3.coeffs()[3], 1, err));
	}

	SECTION("multiply-scalar"){
		const auto scalar = 2;
		const auto num_coeff = poly2.coeffs().size() + poly.coeffs().size()-1;
		const auto poly3 = poly * scalar;
		REQUIRE(poly3.coeffs().size() == poly.coeffs().size());
		for(std::size_t i = 0; i != poly.coeffs().size();++i ){
			REQUIRE(approx_equal(scalar*poly.coeffs()[i], poly3.coeffs()[i], err));
		}
	}

	SECTION("multiply"){
		const auto num_coeff = poly2.coeffs().size() + poly.coeffs().size()-1;
		const auto poly3 = poly * poly2;
		REQUIRE(poly3.coeffs().size() == num_coeff);
		REQUIRE(approx_equal(poly3.coeffs()[0], -4, err));
		REQUIRE(approx_equal(poly3.coeffs()[1],  6, err));
		REQUIRE(approx_equal(poly3.coeffs()[2],  9, err));
		REQUIRE(approx_equal(poly3.coeffs()[3], 29, err));
		REQUIRE(approx_equal(poly3.coeffs()[4],-19, err));
		REQUIRE(approx_equal(poly3.coeffs()[5],  3, err));
	}

}


TEST_CASE("derive"){
	const polynomial poly({-4, 14, -7, 1});
	const polynomial dpoly = derive(poly);
	REQUIRE(dpoly.coeffs() == std::vector<double>({14, -7.0 * 2, 1.0 *3}));

	SECTION("integrate"){
		const auto ipoly = integrate(poly);
		REQUIRE(derive(ipoly).coeffs() == poly.coeffs());


		SECTION("integrate on interval"){
			polynomial poly2({1,1});
			auto integral = integrate(poly2, closed_interval{0,4});
			REQUIRE(approx_equal(integral, 12, 0.0001));
		}
	}
}
