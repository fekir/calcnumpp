//          Copyright Federico Kircheis 2017-2018
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CALCNUM_EULER_HPP_1437731362
#define CALCNUM_EULER_HPP_1437731362

#include <functional>

#include "utils.hpp"
#include "ode.hpp"

namespace calcnum{
	struct explicit_euler_status{
		double t;
		double y;
		double f;
	};
	class explicit_euler{
		cauchy_eq f;
		explicit_euler_status st;
		double h;

	public:
		bool invariant() const {
			return h > 0 && !std::isinf(h) && f;
		}

		explicit explicit_euler(const cauchy_prob& cp, double h_) : f(cp.f), h(h_){
			test_invariant<explicit_euler> _(*this);
			st.t = cp.t0;
			st.y = cp.y0;
			st.f = f(cp.t0,cp.y0);
		}

		explicit_euler& operator++(){
			test_invariant<explicit_euler> _(*this);
			st.y = std::fma(h, st.f, st.y);
			st.t += h;
			st.f = f(st.t, st.y);
			return *this;
		}

		const explicit_euler_status& operator*() const {
			test_invariant<explicit_euler> _(*this);
			return st;
		}
		const explicit_euler_status* operator->() const {
			test_invariant<explicit_euler> _(*this);
			return &st;
		}
	};


	struct implicit_euler_status{
		double t;
		double y;
	};
	class implicit_euler{
		cauchy_eq f;
		double h;
		std::function<cauchy_var(cauchy_var, cauchy_time, cauchy_eq, double)> callback;
		implicit_euler_status st;

	public:
		bool invariant() const {
			return h > 0 && !std::isinf(h) && f && callback;
		}

		explicit implicit_euler(const cauchy_prob& cp, double h_, std::function<cauchy_var(cauchy_var, cauchy_time, cauchy_eq, double)> c_) : f(cp.f), h(h_), callback(c_){
			test_invariant<implicit_euler> _(*this);
			st.t = cp.t0;
			st.y = cp.y0;
		}

		implicit_euler& operator++(){
			test_invariant<implicit_euler> _(*this);
			st.t += h;
			st.y = callback(st.y, st.t, f, h);
			return *this;
		}

		implicit_euler_status& operator*() {
			test_invariant<implicit_euler> _(*this);
			return st;
		}
		implicit_euler_status* operator->() {
			test_invariant<implicit_euler> _(*this);
			return &st;
		}
	};

}

#endif
