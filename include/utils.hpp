#pragma once

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

}

