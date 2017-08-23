#ifndef CALCNUM_NEWTON_HPP_1797832091
#define CALCNUM_NEWTON_HPP_1797832091

#include "intervals.hpp"
#include "utils.hpp"

#include <functional>
#include <cassert>
namespace calcnum{

	// unlike bisection, the definition of status is not that natural, given the Newton formula
	// f(x[k])/f'(x[k]) (equiv to x[k+1]-x[k]) or f(x[k]) can be used as a stopping criteria
	struct newton_status{
		double fx;
		double dfx;
		double x;
	};

	class newton_iter{
		std::function<double(double)> f;
		std::function<double(double)> df;
		newton_status st;

	public:
		explicit newton_iter(const std::function<double(double)>& f_, const std::function<double(double)>& df_, const double& x0) : f(f_), df(df_){
			st.x = x0;
			st.fx = f(st.x);
			st.dfx = df(st.x);
		}
		newton_iter& operator++(){
			auto oldst = st;
			st.x = oldst.x - oldst.fx/oldst.dfx;
			st.fx = f(st.x);
			st.dfx = df(st.x);
			return *this;
		}
		const newton_status& operator*() const {
			return st;
		}
		const newton_status* operator->() const {
			return &st;
		}
	};


	struct secanti_status{
		double xn;
		double xn_1;
		double fxn;
		double fxn_1;
	};

	class secanti_iter{
		std::function<double(double)> f;
		secanti_status st;

	public:
		explicit secanti_iter(const std::function<double(double)>& f_, double x0, double x1) : f(f_){
			st.xn_1 = x0;
			st.xn = x1;
			st.fxn_1 = f(x0);
			st.fxn = f(x1);
		}
		secanti_iter& operator++(){
			const auto x_new = st.xn - ((st.xn - st.xn_1)/(st.fxn - st.fxn_1))*st.fxn;

			// watch ordering!
			st.xn_1 = st.xn;
			st.fxn_1 = st.fxn;

			st.xn = x_new;
			st.fxn = f(x_new);

			return *this;
		}
		const secanti_status& operator*() const {
			return st;
		}
		const secanti_status* operator->() const {
			return &st;
		}
	};
}

#endif
