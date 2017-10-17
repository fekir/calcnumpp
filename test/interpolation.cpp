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

	SECTION("exact"){
		std::vector<double> data_x = {1, 2, 3};
		std::vector<double> data_y = {1, 4, 9};

		const auto pol1 = data_y[0]*lagr_base_poly(std::begin(data_x), std::end(data_x), 0);
		{
			const auto& coeff1 = pol1.coeffs();
			REQUIRE(coeff1[0] == Approx(3));
			REQUIRE(coeff1[1] == Approx(-2.5));
			REQUIRE(coeff1[2] == Approx(0.5));
		}
		const auto pol2 = data_y[1]*lagr_base_poly(std::begin(data_x), std::end(data_x), 1);
		{
			const auto& coeff2 = pol2.coeffs();
			REQUIRE(coeff2[0] == Approx(-12));
			REQUIRE(coeff2[1] == Approx(16));
			REQUIRE(coeff2[2] == Approx(-4));
		}
		const auto pol3 = data_y[2]*lagr_base_poly(std::begin(data_x), std::end(data_x), 2);
		{
			const auto& coeff3 = pol3.coeffs();
			REQUIRE(coeff3[0] == Approx(9));
			REQUIRE(coeff3[1] == Approx(-13.5));
			REQUIRE(coeff3[2] == Approx(4.5));
		}
		auto pol = pol1 + pol2 + pol3;
		{
			const auto& coeff = pol.coeffs();
			REQUIRE(coeff[0] == Approx(0));
			REQUIRE(coeff[1] == Approx(0));
			REQUIRE(coeff[2] == Approx(1));
		}

		auto lagr_pol = lagr_poly(std::begin(data_x), std::end(data_x),std::begin(data_y), std::end(data_y));

		REQUIRE(pol.coeffs().size() == lagr_pol.coeffs().size());
		for(std::size_t i = 0; i != pol.coeffs().size(); ++i){
			REQUIRE(pol.coeffs()[i] == Approx(lagr_pol.coeffs()[i]));
		}
	}

	SECTION("higher order"){

		std::vector<double> data_x = {1, 2, 3};
		std::vector<double> data_y = {1, 8, 27};

		auto pol = lagr_poly(std::begin(data_x), std::end(data_x),std::begin(data_y), std::end(data_y));
		REQUIRE(pol.coeffs()[0] == Approx(6));
		REQUIRE(pol.coeffs()[1] == Approx(-11));
		REQUIRE(pol.coeffs()[2] == Approx(6));
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
	REQUIRE(nodes[0] == Approx(-0.98768));
}

TEST_CASE("spline1"){
	const std::vector<double> data_x = {-2, 0, 1, 2, 3};
	const std::vector<double> data_y = {-8, 0, 1, 8, 27};
	auto spline = create_interp_spline_deg1(std::begin(data_x), std::end(data_x),std::begin(data_y), std::end(data_y));
	for(std::size_t i = 0; i != data_x.size(); ++i){
		REQUIRE(spline(data_x[i]) == Approx(data_y[i]));
	}
}
