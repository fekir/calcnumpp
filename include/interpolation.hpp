//          Copyright Federico Kircheis 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CALCNUM_INTERPOLATION_HPP_7418023189
#define CALCNUM_INTERPOLATION_HPP_7418023189


#include "utils.hpp"
#include "intervals.hpp"
#include "polynomial.hpp"
#include "constants.hpp"

#include <cmath>

#include <iterator>
#include <numeric>
#include <cassert>
#include <vector>
#include <functional>
#include <algorithm>


namespace calcnum{

	template<class InputIt>
	polynomial lagr_base_poly(InputIt first, InputIt last, const std::size_t j)
	{
		const double x_j = first[j];
		polynomial prevpoly(1);
		std::size_t k = 0;
		for(auto it = first; it != last; ++it, ++k){
			if(k == j){
				continue;
			}
			const double x_i = *it;
			polynomial curr_poly({- x_i/(x_j-x_i), 1/(x_j-x_i)});
			prevpoly = prevpoly * curr_poly;
		}

		return prevpoly;
	}

	template<class InputIt>
	polynomial lagr_poly(InputIt first, InputIt last, InputIt first2, InputIt last2)
	{
		polynomial toret;
		std::size_t i{};
		for(auto it = first; it != last; ++it, ++first2, ++i){
			assert(first2 != last2);
			(void)last2;
			double y_i = *first2;
			toret = toret + y_i * lagr_base_poly(first, last, i);
		}
		return toret;
	}


	class chebyshev_nodes{
		open_interval interv;
		std::size_t num_nodes;
	public:
		bool invariant() const{
			return num_nodes > 0 && calcnum::invariant(interv);;
		}
		chebyshev_nodes(const open_interval& interval, std::size_t num_nodes_) : interv(interval), num_nodes(num_nodes_){
			test_invariant<chebyshev_nodes> _(*this);
		}
		double operator[](std::size_t i) const {
			test_invariant<chebyshev_nodes> _(*this);
			assert(i<num_nodes);
			auto toret =  midpoint(interv) - (length(interv)/2)*std::cos(pi*d(2*i+1)/d(2*num_nodes));
			return toret;
		}
		std::size_t size() const noexcept {
			test_invariant<chebyshev_nodes> _(*this);
			return num_nodes;
		}

		class iterator{
			std::size_t i;
			std::reference_wrapper<const chebyshev_nodes> nodes;
		public:

			typedef double value_type;
			typedef double& reference;
			typedef double* pointer;
			typedef std::bidirectional_iterator_tag iterator_category;
			typedef std::ptrdiff_t difference_type;

			bool invariant() const{
				return i <= nodes.get().size();
			}

			iterator(std::size_t i_, const chebyshev_nodes& n) :i(i_), nodes(n){
				test_invariant<iterator> _(*this);
			}
			iterator& operator++ () {
				{
					test_invariant<iterator> _(*this);
					++i;
				}
				return *this;
			}
			iterator& operator+= (std::size_t t) {
				{
					test_invariant<iterator> _(*this);
					i+=t;
				}
				return *this;
			}
			iterator& operator-- () {
				assert(i>0);
				{
					test_invariant<iterator> _(*this);
					--i;
				}
				return *this;
			}

			iterator& operator-= (std::size_t t) {
				assert(i>=t);
				{
					test_invariant<iterator> _(*this);
					i-=t;
				}
				return *this;
			}
			double operator*(){
				test_invariant<iterator> _(*this);
				assert(i < nodes.get().size());
				return nodes.get()[i];
			}

			friend bool operator==(const iterator c1, const iterator c2){
				test_invariant<iterator> _1(c1);
				test_invariant<iterator> _2(c2);
				return c1.i == c2.i;
			}

		};
		iterator begin() const {
			test_invariant<chebyshev_nodes> _(*this);
			return iterator(0, *this);
		}
		iterator end() const {
			test_invariant<chebyshev_nodes> _(*this);
			return iterator(num_nodes, *this);
		}
	};

	inline bool operator!=(const chebyshev_nodes::iterator c1, const chebyshev_nodes::iterator c2){
		return !(c1 == c2);
	}
	inline chebyshev_nodes::iterator operator-- (chebyshev_nodes::iterator it, int) {
		--it;
		return it;
	}
	inline chebyshev_nodes::iterator operator++ (chebyshev_nodes::iterator it, int) {
		++it;
		return it;
	}
	inline chebyshev_nodes::iterator operator- (chebyshev_nodes::iterator it, std::size_t t) {
		it-=t;
		return it;
	}
	inline chebyshev_nodes::iterator operator+ (chebyshev_nodes::iterator it, std::size_t t) {
		it+=t;
		return it;
	}

	class equidistant_nodes_open{
		open_interval interv;
		std::size_t n;
		double h;
	public:
		bool invariant() const{
			return n > 0 && calcnum::invariant(interv);
		}
		equidistant_nodes_open(const open_interval& interval, std::size_t num_nodes) : interv(interval), n(num_nodes){
			h = length(interv)/d(n+1);
			test_invariant<equidistant_nodes_open> _(*this);
		}
		double operator[](std::size_t i) const {
			test_invariant<equidistant_nodes_open> _(*this);
			assert(i<n);
			return std::fma(i+1, h, interv.a);
		}
		std::size_t size() const noexcept {
			test_invariant<equidistant_nodes_open> _(*this);
			return n;
		}

		class iterator{
			std::size_t i;
			std::reference_wrapper<const equidistant_nodes_open> nodes;
		public:

			typedef double value_type;
			typedef double& reference;
			typedef double* pointer;
			typedef std::bidirectional_iterator_tag iterator_category;
			typedef std::size_t difference_type;

			bool invariant() const{
				return i <= nodes.get().size();
			}
			iterator(std::size_t i_, const equidistant_nodes_open& n_) :i(i_), nodes(n_){
				test_invariant<iterator> _(*this);
			}
			iterator& operator++ () {
				{
					test_invariant<iterator> _(*this);
					++i;
				}
				return *this;
			}
			iterator& operator+= (std::size_t t) {
				{
					test_invariant<iterator> _(*this);
					i+=t;
				}
				return *this;
			}
			iterator& operator-- () {
				assert(i>0);
				{
					test_invariant<iterator> _(*this);
					--i;
				}
				return *this;
			}

			iterator& operator-= (std::size_t t) {
				assert(i>=t);
				{
					test_invariant<iterator> _(*this);
					i-=t;
				}
				return *this;
			}

			double operator*() const {
				test_invariant<iterator> _(*this);
				return nodes.get()[i];
			}

			friend bool operator==(const iterator c1, const iterator c2){
				test_invariant<iterator> _1(c1);
				test_invariant<iterator> _2(c2);
				return c1.i == c2.i;
			}
		};
		iterator begin() const {
			test_invariant<equidistant_nodes_open> _1(*this);
			return iterator(0, *this);
		}
		iterator end() const {
			test_invariant<equidistant_nodes_open> _1(*this);
			return iterator(n, *this);
		}
	};
	inline bool operator!=(equidistant_nodes_open::iterator c1, equidistant_nodes_open::iterator c2){
		return !(c1 == c2);
	}
	inline equidistant_nodes_open::iterator operator-- (equidistant_nodes_open::iterator it, int) {
		--it;
		return it;
	}
	inline equidistant_nodes_open::iterator operator++ (equidistant_nodes_open::iterator it, int) {
		++it;
		return it;
	}
	inline equidistant_nodes_open::iterator operator- (equidistant_nodes_open::iterator it, std::size_t t) {
		it-=t;
		return it;
	}
	inline equidistant_nodes_open::iterator operator+ (equidistant_nodes_open::iterator it, std::size_t t) {
		it+=t;
		return it;
	}

	class equidistant_nodes_closed{
		closed_interval interv;
		std::size_t n;
		double h;
	public:
		bool invariant() const{
			return n > 1 && calcnum::invariant(interv);
		}
		equidistant_nodes_closed(const closed_interval& interval, std::size_t num_nodes) : interv(interval), n(num_nodes){
			h = length(interv)/d(n-1);
			test_invariant<equidistant_nodes_closed> _1(*this);
		}
		double operator[](std::size_t i) const {
			test_invariant<equidistant_nodes_closed> _1(*this);
			assert(i<n);
			return std::fma(i, h, interv.a);
		}
		std::size_t size() const noexcept {
			test_invariant<equidistant_nodes_closed> _1(*this);
			return n;
		}

		class iterator{
			std::size_t i;
			std::reference_wrapper<const equidistant_nodes_closed> nodes;
		public:

			typedef double value_type;
			typedef double& reference;
			typedef double* pointer;
			typedef std::bidirectional_iterator_tag iterator_category;
			typedef std::ptrdiff_t difference_type;

			bool invariant() const{
				return i <= nodes.get().size();
			}

			iterator(std::size_t i_, const equidistant_nodes_closed& n_) :i(i_), nodes(n_){
				test_invariant<iterator> _1(*this);
			}
			iterator& operator++ () {
				{
					test_invariant<iterator> _(*this);
					++i;
				}
				return *this;
			}
			iterator& operator+= (std::size_t t) {
				{
					test_invariant<iterator> _(*this);
					i+=t;
				}
				return *this;
			}
			iterator& operator-- () {
				assert(i>0);
				{
					test_invariant<iterator> _(*this);
					--i;
				}
				return *this;
			}

			iterator& operator-= (std::size_t t) {
				assert(i>=t);
				{
					test_invariant<iterator> _(*this);
					i-=t;
				}
				return *this;
			}
			double operator*() const {
				test_invariant<iterator> _1(*this);
				return nodes.get()[i];
			}

			friend bool operator==(const iterator c1, const iterator c2){
				test_invariant<iterator> _1(c1);
				test_invariant<iterator> _2(c2);
				return c1.i == c2.i;
			}
		};
		iterator begin() const {
			test_invariant<equidistant_nodes_closed> _1(*this);
			return iterator(0, *this);
		}
		iterator end() const {
			test_invariant<equidistant_nodes_closed> _1(*this);
			return iterator(n, *this);
		}
	};
	inline bool operator!=(equidistant_nodes_closed::iterator c1, equidistant_nodes_closed::iterator c2){
		return !(c1 == c2);
	}
	inline equidistant_nodes_closed::iterator operator-- (equidistant_nodes_closed::iterator it, int) {
		--it;
		return it;
	}
	inline equidistant_nodes_closed::iterator operator++ (equidistant_nodes_closed::iterator it, int) {
		++it;
		return it;
	}
	inline equidistant_nodes_closed::iterator operator- (equidistant_nodes_closed::iterator it, std::size_t t) {
		it-=t;
		return it;
	}
	inline equidistant_nodes_closed::iterator operator+ (equidistant_nodes_closed::iterator it, std::size_t t) {
		it+=t;
		return it;
	}

	inline equidistant_nodes_closed make_node_gen(const closed_interval& interval, std::size_t num_nodes){
		return equidistant_nodes_closed(interval, num_nodes);
	}
	inline equidistant_nodes_open make_node_gen(const open_interval& interval, std::size_t num_nodes){
		return equidistant_nodes_open(interval, num_nodes);
	}



	struct poly_interv{
		polynomial poly;
		double b; // no need to store the whole interval, just right limit. Left limit is given by previous interval or -inf
	};

	inline bool poly_interv_are_adjacent(const poly_interv& p_i1, const poly_interv& p_i2){
		return p_i1.b < p_i2.b;
	}

	// FIXME: if two adjacent splines are equals (ie same coeffs), we can use the same spline on a bigger interval
	class interp_spline_deg1{
		std::vector<poly_interv> p_is;
	public:
		interp_spline_deg1() = default;
		bool invariant() const {
			bool grade_less_eq2 = !p_is.empty() && std::all_of(std::begin(p_is), std::end(p_is), [](const poly_interv& p_i){return p_i.poly.coeffs().size() <= 2;});
			bool interval_are_adjacent_sorted = (std::is_sorted_until(std::begin(p_is), std::end(p_is), poly_interv_are_adjacent) == std::end(p_is));
			bool last_cond = std::isinf(p_is.back().b) && signum(p_is.back().b) == sign::greaterzero;

			bool border_cond = true;
			for(auto it = std::begin(p_is); it != std::end(p_is)-1; ++it){
				auto jt = it+1;
				border_cond = border_cond && approx_equal(it->poly(it->b), jt->poly(it->b),0.0001);
			}
			return grade_less_eq2 && interval_are_adjacent_sorted && last_cond && border_cond;
		}
		explicit interp_spline_deg1(const std::vector<poly_interv>& p_is_) : p_is(p_is_){
			test_invariant<interp_spline_deg1> _(*this);
		}
		double operator()(double x) const {
			test_invariant<interp_spline_deg1> _(*this);
			auto it = std::find_if(std::begin(p_is), std::end(p_is), [x](const poly_interv& v){return x <= v.b;});
			if(it == std::end(p_is)){
				assert(std::isnan(x) && "we should always find an interval");
				return x;
			}
			return it->poly(x);
		}
	};


	template<class InputIt>
	interp_spline_deg1 create_interp_spline_deg1(InputIt first, InputIt last, InputIt first2, InputIt last2){
		std::vector<poly_interv> toreturn;
		preallocate_if_randomaccess(first, last, toreturn);
		auto it = first;
		auto jt = first +1;

		auto it_y = first2;
		auto jt_y = first2 +1;
		for(; jt != last; ++it, ++jt,  ++it_y, ++jt_y){
			assert(jt_y != last2);
			(void)last2;
			polynomial p = lagr_poly(it, jt+1, it_y, jt_y+1);
			assert(p.coeffs().size() <= 2 && "should create only linear splines");
			poly_interv p_i;
			p_i.poly = p;
			p_i.b = *jt;
			assert(toreturn.empty() || (approx_equal(p(*it), toreturn.back().poly(*it), 0.0001) && "border condition not true"));
			toreturn.push_back(p_i);
		}
		assert(!toreturn.empty());
		// extend intervals of spline
		toreturn.back().b = std::numeric_limits<double>::infinity();
		return interp_spline_deg1(toreturn);
	}
}

#endif
