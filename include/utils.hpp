#ifndef CALCNUM_UTILS_HPP_7922187677
#define CALCNUM_UTILS_HPP_7922187677

#include <cassert>
#include <cmath>

#include <vector>
#include <algorithm>
#include <numeric>

namespace calcnum{

	// operator== should, most of the times, be avoided between floating point values
	inline bool approx_equal(double a, double b, double err){
		assert(err>0 && "error for comparing values should be strict positive");
		return a-err <= b && a+err >= b;
	}

	// converting to double, more readable in complex formulas instead of explicit cast
	template<class T>
	double d(T v){
		return static_cast<double>(v);
	}


	/// Returns the sign (less zero, zero, greater zero) of a floating point number
	/// It does not work correctly with nan and does not distinguish between +0.0/-0.0, see copysign/signbit if you need the sign of those types too
	/// It otherwise should work with any double, +inf and -inf included
	enum sign : int {
		lesszero = -1, zero = 0, greaterzero =1
	};
	inline sign signum(double d){
		return (d<0) ? sign::lesszero : (d>0) ? sign::greaterzero : sign::zero;
	}

	inline std::vector<double> calculate_convergency(std::vector<double> err, std::vector<double> step){
		std::transform(err.begin(), err.end(), err.begin(), [](double d){return std::log(std::fabs(d));});
		std::transform(step.begin(), step.end(), step.begin(), [](double d){return std::log(std::fabs(d));});

		// FIXME: can maybe do in-place
		std::vector<double> diff_lerr;
		std::adjacent_difference(err.begin(), err.end(), std::back_inserter(diff_lerr));

		std::vector<double> diff_step;
		std::adjacent_difference(step.begin(), step.end(), std::back_inserter(diff_step));

		std::vector<double> conv;
		conv.reserve(diff_lerr.size()-1);
		std::transform(diff_lerr.begin()+1, diff_lerr.end(), diff_step.begin()+1, std::back_inserter(conv), [](double a ,double b){return a/b;});
		return conv;
	}

	/// Helper class for verifying invariant of an object
	template<class invar_obj>
	struct test_invariant {
		const invar_obj& obj;
		test_invariant(const invar_obj& obj_) : obj(obj_)  {
			assert(obj.invariant());
		}
		~test_invariant() {
			assert(obj.invariant());
		}
	};
}

#endif
