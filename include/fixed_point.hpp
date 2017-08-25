#ifndef CALCNUM_FIXED_POINT_HPP_0801338869
#define CALCNUM_FIXED_POINT_HPP_0801338869

#include "intervals.hpp"
#include "utils.hpp"

#include <functional>
#include <cassert>
namespace calcnum{

	// aitken
	// pto fisso normale

	struct fixed_point_status{
		double f;
	};
	class fixed_point{
		std::function<double(double)> f;
		fixed_point_status st;
	public:
		explicit fixed_point(const std::function<double(double)>& f_, double x0) : f(f_) {
			st.f = f(x0);
		}

		fixed_point& operator++(){
			st.f = f(st.f);
			return *this;
		}

		const fixed_point_status& operator*() const {
			return st;
		}
		const fixed_point_status* operator->() const {
			return &st;
		}
	};
}

#endif
