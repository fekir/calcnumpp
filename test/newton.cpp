#include <catch.hpp>

#include "newton.hpp"


#include <functional>

namespace{
	// functions with a 0 in x=2
	double fun_new(double x){
		return x*x*x-3*x-2;
	}
	double dfun_new(double x){
		return 3*x*x-3;
	}
	double fun_sec(double x){
		return x+0.5-std::sin(x);
	}

	void simple_newton_test(const std::function<double(double)>& f, const std::function<double(double)>& df, double x0){
		calcnum::newton_iter it(f, df, x0);
		auto oldst = *it;
		while(std::fabs(it->fx) > 0.0000001){
			++it;
			const auto& st = *it;
			oldst = st;
		}
	}
	void simple_secanti_test(const std::function<double(double)>& f, double x0, double x1){
		calcnum::secanti_iter it(f, x0, x1);
		auto oldst = *it;
		while(std::fabs(it->fxn) > 0.00001){
			++it;
			const auto& st = *it;
			oldst = st;
		}
	}
}

TEST_CASE("newton"){
	simple_newton_test(fun_new, dfun_new, 3);
}



TEST_CASE("secanti"){
	simple_secanti_test(fun_sec,-1,-2);
}
