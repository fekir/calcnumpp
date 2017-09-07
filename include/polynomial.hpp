//          Copyright Federico Kircheis 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CALCNUM_POLYNOMIAL_HPP_9495053367
#define CALCNUM_POLYNOMIAL_HPP_9495053367


#include "utils.hpp"
#include "intervals.hpp"

#include <cmath>

#include <iterator>
#include <numeric>
#include <cassert>
#include <vector>
#include <functional>
#include <algorithm>


namespace calcnum{
	template<class iter>
	double eval_with_hoerners_rule(iter first, iter last, double x0){
		auto rfirst = std::reverse_iterator<iter>(last);
		auto rlast = std::reverse_iterator<iter>(first);

		// result_{n+1}=result_n*x0+coeff_n
		auto acc = [x0](double result, double coeff){
			return std::fma(result, x0, coeff);
		};
		return std::accumulate(rfirst, rlast, 0.0, acc);
	}

	template<class OutputIt, class numeric>
	void assert_poly_0_init(OutputIt out, numeric size){
		(void)out;
		for(numeric i = 0; i != size ; ++i){
			assert(approx_equal(out[i], 0, 0.000001));
		}
	}

	// out needs to be 0-init!
	// FIXME: there are more efficient versions
	// do not like function signature very much
	template<class InputIt1, class InputIt2, class OutputIt>
	void multiply_poly(InputIt1 first1, InputIt1 last1, InputIt2 first2, InputIt2 last2, OutputIt out)
	{
		assert_poly_0_init(out, std::distance(first1, last1 ) + std::distance(first2, last2 ) -1);
		for(auto it = first1; it != last1; ++it){
			const auto i = std::distance(first1, it ) ;
			for(auto jt = first2; jt != last2; ++jt){
				const auto j = std::distance(first2, jt ) ;
				// out[i+j] += it*jt;
				out[i+j] = std::fma(*it, *jt, out[i+j]);
			}
		}
	}

	template<class InputIt1, class InputIt2, class OutputIt>
	void multiply_poly(InputIt1 first1, InputIt1 last1, double c, OutputIt out)
	{
		std::transform(first1, last1, out, [c](double v){return v*c;});
	}

	template<class InputIt1, class InputIt2, class OutputIt>
	void sum_poly(InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt out)
	{
		std::transform(first1, last1, first2, out, [](double c1, double c2){return c1+c2;});
	}


	// could specialize, when size is known at compile time, with std::array<double, N>
	class polynomial{
		std::vector<double> coeff = {0};
	public:

		static void normalize(std::vector<double>& coeff){
			auto it = coeff.rbegin();
			for((void)it; it != coeff.rend(); ++it){
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
				if(*it != 0){
#pragma GCC diagnostic pop
					break;
				}
			}
			if(it == coeff.rend()){
				coeff = {0};
			} else{
				coeff.erase(it.base(), coeff.end());
			}
		}
		bool invariant() const{
			if(coeff.size() > 1){
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
				return coeff.back() != 0;
#pragma GCC diagnostic pop
			}
			return !coeff.empty();
		}
		polynomial() = default;
		explicit polynomial(double c0) : coeff{c0}{
			test_invariant<polynomial> _(*this);
		}
		explicit polynomial(const std::vector<double>& coeff_) : coeff(coeff_){
			normalize(coeff);
			test_invariant<polynomial> _(*this);
		}

		friend polynomial operator+(polynomial lhs, const polynomial& rhs){
			{
				test_invariant<polynomial> _(lhs);
				test_invariant<polynomial> __(rhs);
				lhs.coeff.resize(std::max(lhs.coeff.size(), rhs.coeff.size()));
				const auto minsize = std::min(lhs.coeff.size(), rhs.coeff.size());
				std::transform(std::begin(lhs.coeff), std::begin(lhs.coeff)+minsize, std::begin(rhs.coeff), std::begin(lhs.coeff), std::plus<double>());
				normalize(lhs.coeff);
			}
			return lhs;
		}

		friend polynomial operator+(polynomial lhs, double coeff0){
			{
				test_invariant<polynomial> _(lhs);
				lhs.coeff[0] += coeff0;
			}
			return lhs;
		}
		friend polynomial operator+(double coeff0, polynomial rhs){
			{
				test_invariant<polynomial> _(rhs);
				rhs.coeff[0] += coeff0;
			}
			return rhs;
		}

		friend polynomial operator*(const polynomial& lhs, const polynomial& rhs){
			test_invariant<polynomial> _(lhs);
			test_invariant<polynomial> __(rhs);
			polynomial toret;
			{
				test_invariant<polynomial> ___(toret);
				toret.coeff.resize(lhs.coeff.size()+rhs.coeff.size()-1);
				multiply_poly(std::begin(lhs.coeff), std::end(lhs.coeff), std::begin(rhs.coeff), std::end(rhs.coeff), std::begin(toret.coeff));
				normalize(toret.coeff);
			}
			return toret;

		}
		friend polynomial operator*(polynomial lhs, double coeff0){
			{
				test_invariant<polynomial> _(lhs);
				std::transform(std::begin(lhs.coeff), std::end(lhs.coeff), std::begin(lhs.coeff), [coeff0](double v){return v*coeff0;});
				normalize(lhs.coeff);
			}
			return lhs;
		}
		friend polynomial operator*(double coeff0, polynomial lhs){
			{
				test_invariant<polynomial> _(lhs);
				std::transform(std::begin(lhs.coeff), std::end(lhs.coeff), std::begin(lhs.coeff), [coeff0](double v){return v*coeff0;});
				normalize(lhs.coeff);
			}
			return lhs;
		}
		friend polynomial derive(polynomial p){
			{
				test_invariant<polynomial> _(p);
				if(p.coeffs().size() == 1){
					return polynomial(0);
				}
				for(std::size_t i = 1; i != p.coeff.size(); ++i){
					p.coeff[i] = p.coeff[i]*d(i);
				}
				p.coeff.erase(p.coeff.begin());
			}
			return p;
		}

		friend polynomial integrate(polynomial p){
			{
				test_invariant<polynomial> _(p);
				for(std::size_t i = 0; i != p.coeff.size(); ++i){
					p.coeff[i] = p.coeff[i]/d(i+1);
				}
				p.coeff.insert(p.coeff.begin(), 0);
			}
			return p;
		}

		double operator()(double x) const {
			test_invariant<polynomial> _(*this);
			return eval_with_hoerners_rule(std::begin(coeff), std::end(coeff), x);
		}

		const std::vector<double>& coeffs() const{
			test_invariant<polynomial> _(*this);
			return coeff;
		}
	};


	inline double integrate(polynomial p, const closed_interval interv){
		p = integrate(std::move(p));
		return p(interv.b) - p(interv.a);
	}
	inline double integrate(polynomial p, const open_interval interv){
		p = integrate(std::move(p));
		return p(interv.b) - p(interv.a);
	}


	inline void print_polynomial(std::ostream& os, const polynomial pol){
		os << "poly: ";
		for(const auto coeff: pol.coeffs()){
			os << coeff << ", ";
		}
		os << "\n";
	}

}

#endif
