#include <catch.hpp>

#include "newton.hpp"


#include <functional>

namespace{
	// functions with a 0 in x=2
	double fun_new(double x){
		return x*x*x - 3*x - 2;
	}
	double dfun_new(double x){
		return 3*x*x - 3;
	}

	double fun_new_doubleroot(double x){
		return std::exp(x) - x*x - std::sin(x)-1;
	}
	double dfun_new_doubleroot(double x){
		return std::exp(x) - 2*x - std::cos(x);
	}

	std::vector<double> simple_newton_test(const std::function<double(double)>& f, const std::function<double(double)>& df, double x0, double sol){
		std::vector<double> err;
		std::vector<double> h;

		calcnum::newton_iter it(f, df, x0);
		auto oldst = *it;
		while(std::fabs(it->fx) > 0.0000000001){
			++it;
			err.push_back(std::abs(sol-it->x));
			h.push_back(std::abs(sol-oldst.x));
			const auto& st = *it;
			oldst = st;
		}

		return calcnum::calculate_convergency(err,h);
	}
	std::vector<double>  simple_secanti_test(const std::function<double(double)>& f, double x0, double x1, double sol){
		std::vector<double> err;
		std::vector<double> h;
		calcnum::secanti_iter it(f, x0, x1);
		auto oldst = *it;
		while(std::fabs(it->fxn) > 0.0000000001){
			++it;
			err.push_back(std::abs(sol-it->xn));
			h.push_back(std::abs(sol-oldst.xn));
			const auto& st = *it;
			oldst = st;
		}
		return calcnum::calculate_convergency(err,h);
	}
}

TEST_CASE("newton"){
    SECTION("simple_fun"){
        auto conv = simple_newton_test(fun_new, dfun_new, 7, 2);
        (void) conv; // FIXME: REQUIRE conv ~ 2
    }
    SECTION("complex_fun"){
        SECTION("quad_conv"){
            auto conv = simple_newton_test(fun_new_doubleroot, dfun_new_doubleroot, 0.1,0);
            (void) conv; // FIXME: REQUIRE conv ~ 1
        }
        SECTION("lin_conv"){
            auto conv = simple_newton_test(fun_new_doubleroot, dfun_new_doubleroot, 1.2,1.27970133100099630500239591776735167562639703793577);
            (void) conv; // FIXME: REQUIRE conv ~ 2
        }
    }
}



TEST_CASE("secanti"){
    auto conv = simple_secanti_test(fun_new,7,7.01, 2);
    (void) conv; // FIXME: REQUIRE conv ~ 2
}
