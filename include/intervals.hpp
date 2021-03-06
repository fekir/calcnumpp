//          Copyright Federico Kircheis 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CALCNUM_INTERVAL_HPP_0355809029
#define CALCNUM_INTERVAL_HPP_0355809029

#include <cassert>

namespace calcnum{
	struct closed_interval{
		double a;
		double b;
	};

	inline bool invariant(const closed_interval& interval) {
		return interval.a <= interval.b;
	}

	inline double midpoint(const closed_interval& interval){
		assert(invariant(interval));
		return (interval.a + interval.b)/2;
	}

	inline double length(const closed_interval& interval){
		assert(invariant(interval));
		return interval.b - interval.a;
	}

	struct open_interval{
		double a;
		double b;
	};

	inline bool invariant(const open_interval& interval) {
		return interval.a <= interval.b;
	}

	inline double midpoint(const open_interval& interval){
		assert(invariant(interval));
		return (interval.a + interval.b)/2;
	}

	inline double length(const open_interval& interval){
		assert(invariant(interval));
		return interval.b - interval.a;
	}
}

#endif
