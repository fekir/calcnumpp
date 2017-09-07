//          Copyright Federico Kircheis 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#ifndef CALCNUM_UTILS_HPP_7922187677
#define CALCNUM_UTILS_HPP_7922187677

#include <cassert>
#include <cmath>

#include <vector>
#include <algorithm>
#include <numeric>

namespace calcnum{

	template<class iter, class cont>
	void preallocate_if_randomaccess(const iter begin, const iter end, cont& container, int diff = 0) {
		using iter_type = typename std::iterator_traits<iter>::iterator_category;
		if(std::is_same<iter_type,std::random_access_iterator_tag>::value){
			auto size = std::distance(begin, end);
			if (diff >=0){
				size += diff;
			}else{
				size = std::max(size, static_cast<decltype(size)>(-diff)) + diff;
			}
			container.reserve(static_cast<decltype(container.size())>(size));
		}
	}

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


	inline std::vector<double> calculate_convergency(std::vector<double> err){
		assert(!err.empty() && "makes no sense to analyze empty vector");
		const auto err_size = err.size();
		std::transform(err.begin(), err.end(), err.begin(), [](double d){return std::log(std::fabs(d));});

		auto& diff_lerr = err;
		std::adjacent_difference(err.begin(), err.end(), diff_lerr.begin());
		diff_lerr.erase(diff_lerr.begin());
		assert(diff_lerr.size() == err_size-1);

		const auto diff_lerr_size = diff_lerr.size();

		auto& conv = diff_lerr;
		std::transform(diff_lerr.begin(), diff_lerr.end()-1, diff_lerr.begin()+1, conv.begin(),
		               [](double a ,double b){return b/a;}
		);
		conv.pop_back();
		assert(conv.size() == diff_lerr_size-1);

		return conv;
	}

	inline std::vector<double> calculate_convergency_no_cancellation(std::vector<double> err){
		assert(!err.empty() && "makes no sense to analyze empty vector");
		const auto err_size = err.size();
		std::transform(err.begin(), err.end(), err.begin(), [](double d){return std::log(std::fabs(d));});

		auto& conv = err;
		std::transform(err.begin(), err.end()-1, err.begin()+1, conv.begin(),
		               [](double a ,double b){return b/a;}
		);
		conv.pop_back();
		assert(conv.size() == err_size-1);
		return conv;
	}

	inline std::vector<double> clear_from_inf_nan(std::vector<double> c){
		c.erase(std::remove_if(c.begin(), c.end(), [](double d){ return std::isinf(d) || std::isnan(d);}), c.end());
		return c;
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

	// Add Neumaier alternative algorithm
	class kahan_sum{
		double sum = 0;
		double compensation = 0;
	public:
		kahan_sum() = default;
		explicit kahan_sum(double val) : sum(val){}

		kahan_sum& operator+=(const double val){

			auto diff = val - compensation;
			auto tmp = sum + diff;
			compensation = (tmp - sum) - diff;
			sum = tmp;

			return *this;
		}
		explicit operator double () const {
			return sum;
		}
	};
	inline kahan_sum operator+(kahan_sum lhs, const double rhs){
		return lhs += rhs;
	}
	inline kahan_sum operator+(const double lhs, kahan_sum rhs){
		return rhs += lhs;
	}
	inline kahan_sum& operator-=(kahan_sum& lhs, const double rhs){
		return lhs += (-rhs);
	}
	inline kahan_sum operator-(kahan_sum lhs, const double rhs){
		return lhs += (-rhs);
	}
	inline kahan_sum operator-(const double lhs, kahan_sum rhs){
		return rhs += (-lhs);
	}
	inline kahan_sum& operator++(kahan_sum& lhs){
		return lhs += 1;
	}
	inline kahan_sum operator++(kahan_sum lhs,int){
		return lhs += 1;
	}
	inline kahan_sum& operator--(kahan_sum& lhs){
		return lhs += (-1);
	}
	inline kahan_sum operator--(kahan_sum lhs,int){
		return lhs += (-1);
	}
}

#endif
