#include <catch.hpp>

#include "newton.hpp"
#include "statistics.hpp"

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

	std::vector<double> simple_newton_test(const std::function<double(double)>& f, const std::function<double(double)>& df, double x0, double sol, std::size_t max_iter = 100){
		std::vector<double> err;

		calcnum::newton_iter it(f, df, x0);
		auto oldst = *it;
		while(std::fabs(it->fx) > 0.0000000001){
			++it;
			err.push_back(std::abs(sol-it->x));
			const auto& st = *it;
			oldst = st;
		}
		return err;
	}
	std::vector<double> simple_secanti_test(const std::function<double(double)>& f, double x0, double x1, double sol, std::size_t max_iter = 100){
		std::vector<double> err;
		calcnum::secanti_iter it(f, x0, x1);
		auto oldst = *it;
		while(std::fabs(it->xn - sol) > 0.0000000001){
			++it;
			err.push_back(std::abs(sol-it->xn));
			const auto& st = *it;
			oldst = st;
		}
		return err;
	}
}

TEST_CASE("newton"){
    SECTION("simple_fun"){
        const auto err = simple_newton_test(fun_new, dfun_new, 7, 2);
        const auto conv = calcnum::calculate_convergency(err);
        auto res = calcnum::analyze_data(calcnum::clear_from_inf_nan(conv));
        REQUIRE_FALSE(calcnum::is_outlier(res, 2));
    }
    SECTION("complex_fun"){
        SECTION("quad_conv"){
            const auto err = simple_newton_test(fun_new_doubleroot, dfun_new_doubleroot, 0.1,0);
            const auto conv = calcnum::calculate_convergency(err);
            auto res = calcnum::analyze_data(calcnum::clear_from_inf_nan(conv));
            REQUIRE_FALSE(calcnum::is_outlier(res, 1));
        }
        SECTION("lin_conv"){
            auto err = simple_newton_test(fun_new_doubleroot, dfun_new_doubleroot, 1.2,1.27970133100099630500239591776735167562639703793577);
            const auto conv = calcnum::calculate_convergency(err);
            auto res = calcnum::analyze_data(calcnum::clear_from_inf_nan(conv));
            REQUIRE_FALSE(calcnum::is_outlier(res, 2));
        }
    }
}



TEST_CASE("secanti"){
    SECTION("simple_fun"){
        const auto err = simple_secanti_test(fun_new, 7, 7.01, 2);
        const auto conv = calcnum::calculate_convergency(err);
        auto res = calcnum::analyze_data(calcnum::clear_from_inf_nan(conv));

        REQUIRE_FALSE(calcnum::is_outlier(res, 1.618));
    }
    SECTION("complex_fun"){
        SECTION("quad_conv"){
            const auto err = simple_secanti_test(fun_new_doubleroot, 0.1, 0.2, 0);
            const auto conv = calcnum::calculate_convergency(err);
            auto res = calcnum::analyze_data(calcnum::clear_from_inf_nan(conv));
            REQUIRE_FALSE(calcnum::is_outlier(res, 1));
        }
        SECTION("lin_conv"){
            const auto err = simple_secanti_test(fun_new_doubleroot, 1.2, 1.1, 1.27970133100099630500239591776735167562639703793577);
            const auto conv = calcnum::calculate_convergency(err);
            auto res = calcnum::analyze_data(calcnum::clear_from_inf_nan(conv));
            REQUIRE_FALSE(calcnum::is_outlier(res, 1.618));
        }
    }
}
