//          Copyright Federico Kircheis 2017-2018
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <catch.hpp>

#include "euler.hpp"

#include "fixed_point.hpp"

#include <cmath>
struct ode_to_solve{
	std::function<double(double,double)> f;
	std::function<double(double)> sol;
	double t0;
	double tmax;
	double y0;
};

using namespace calcnum;
TEST_CASE("explicit euler"){
	cauchy_prob cp = {
	    [](double t, double y){return t*t *(1-3*y);},
	    cauchy_time(0),
	    cauchy_var(2)
	};
	auto sol = [](double t){ return (1+5*std::exp(-t*t*t))/3;};
	double tmax = 3;

	const double h = 0.01;
	calcnum::explicit_euler ee(cp, h);
	while(ee->t < tmax){
		REQUIRE(sol(ee->t) == Approx(ee->y).epsilon(h));
		++ee;
	}
}

TEST_CASE("implicit euler"){
/*
	cauchy_prob cp0 = {
		[](double t, double y){return 3*(y-t);},
		cauchy_time(0),
		cauchy_var(1.0/3)
	};
	auto sol = [](double t){ return t+1.0/3;},
	double tmax = 5;
	cauchy_prob cp1 = {
		[](double t, double y){(void)t;return -50*y;},
		cauchy_time(0),
		cauchy_var(1)
	};
	cont auto sol = [](double t){ return std::exp(exp(-50*t));},
	double tmax = 2;
*/
	cauchy_prob cp = {
	    [](double t, double y){return t*t *(1-3*y);},
	    cauchy_time(0),
	    cauchy_var(2)
	};
	double tmax = 3;
	auto sol = [](double t){ return (1+5*std::exp(-t*t*t))/3;};

	// callback used in the implicit step
	auto callback = [](cauchy_var yn, cauchy_time tn1, cauchy_eq f, double h){
		auto fun = [yn, tn1, f, h](cauchy_var y){
			return h* f(tn1,y) + yn;
		};
		calcnum::fixed_point fp(fun, yn);
		for(int i = 0; i != 100; ++i){
			++fp;
		}
		return fp->f;
	};

	const double h = 0.01;
	calcnum::implicit_euler ee(cp, h, callback);
	double err = 0;
	while(ee->t < tmax){
		REQUIRE(sol(ee->t) == Approx(ee->y).epsilon(h));
		err = std::fmax(err, std::fabs(ee->y - sol(ee->t)));
		(void)err;
		++ee;
	}

}
