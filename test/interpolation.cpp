//          Copyright Federico Kircheis 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <catch.hpp>

#include "utils.hpp"

#include "interpolation.hpp"

#include <functional>
#include <iterator>

using namespace calcnum;

TEST_CASE("lagrange"){
	const auto err = 0.0001;


	SECTION("exact"){
		std::vector<double> data_x = {1, 2, 3};
		std::vector<double> data_y = {1, 4, 9};

		const auto pol1 = data_y[0]*lagr_base_poly(std::begin(data_x), std::end(data_x), 0);
		{
			const auto& coeff1 = pol1.coeffs();
			REQUIRE(approx_equal(coeff1[0], 3, err));
			REQUIRE(approx_equal(coeff1[1], -2.5, err));
			REQUIRE(approx_equal(coeff1[2], 0.5, err));
		}
		const auto pol2 = data_y[1]*lagr_base_poly(std::begin(data_x), std::end(data_x), 1);
		{
			const auto& coeff2 = pol2.coeffs();
			REQUIRE(approx_equal(coeff2[0], -12, err));
			REQUIRE(approx_equal(coeff2[1], 16, err));
			REQUIRE(approx_equal(coeff2[2], -4, err));
		}
		const auto pol3 = data_y[2]*lagr_base_poly(std::begin(data_x), std::end(data_x), 2);
		{
			const auto& coeff3 = pol3.coeffs();
			REQUIRE(approx_equal(coeff3[0], 9, err));
			REQUIRE(approx_equal(coeff3[1], -13.5, err));
			REQUIRE(approx_equal(coeff3[2], 4.5, err));
		}
		auto pol = pol1 + pol2 + pol3;
		{
			const auto& coeff = pol.coeffs();
			REQUIRE(approx_equal(coeff[0], 0, err));
			REQUIRE(approx_equal(coeff[1], 0, err));
			REQUIRE(approx_equal(coeff[2], 1, err));
		}

		auto lagr_pol = lagr_poly(std::begin(data_x), std::end(data_x),std::begin(data_y), std::end(data_y));
		REQUIRE(pol.coeffs() == lagr_pol.coeffs());
	}

	SECTION("higher order"){

		std::vector<double> data_x = {1, 2, 3};
		std::vector<double> data_y = {1, 8, 27};

		auto pol = lagr_poly(std::begin(data_x), std::end(data_x),std::begin(data_y), std::end(data_y));
		REQUIRE(pol.coeffs() == std::vector<double>({6, -11, 6}));
	}

}

TEST_CASE("chebyshev"){
	std::size_t n = 10;
	std::vector<double> nodes;
	nodes.resize(n);

	open_interval interv{-1,1};
	chebyshev_nodes cb(interv, n);
	std::copy(std::begin(cb), std::end(cb), std::begin(nodes));

	REQUIRE(std::is_sorted(std::begin(nodes), std::end(nodes)));
	REQUIRE(*std::min_element(std::begin(cb), std::end(cb)) > interv.a);
	REQUIRE(*std::min_element(std::begin(cb), std::end(cb)) < interv.b);
	REQUIRE(approx_equal(nodes[0],-0.98768, 0.0001 ));
}

TEST_CASE("spline1"){
	const auto err = 0.0001;
	const std::vector<double> data_x = {-2, 0, 1, 2, 3};
	const std::vector<double> data_y = {-8, 0, 1, 8, 27};
	auto spline = create_interp_spline_deg1(std::begin(data_x), std::end(data_x),std::begin(data_y), std::end(data_y));
	for(std::size_t i = 0; i != data_x.size(); ++i){
		REQUIRE(approx_equal(spline(data_x[i]), data_y[i], err));
	}
}
