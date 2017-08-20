#ifndef CALCNUM_UTILS_HPP_7922187677
#define CALCNUM_UTILS_HPP_7922187677

#include <cassert>

namespace calcnum{

	// operator== should be avoided between floating point values
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
		return d<0 ? sign::lesszero : d>0 ? sign::greaterzero : sign::zero;
	}

}

#endif
