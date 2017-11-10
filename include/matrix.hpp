#pragma once

#include "utils.hpp"

#include <cmath>
#include <array>
#include <utility>
#include <numeric>
#include <functional>
#include <ostream>

namespace calcnum{
	// not to be confused with std::vector
	template<std::size_t N>
	class vector{
	public:
		std::array<double, N> val = {}; // 0-init
		vector() = default;
		explicit vector( const double (& val_) [N]) {
			std::copy(std::begin(val_), std::end(val_), std::begin(val));
		}
		double& operator [](std::size_t idx) {
			return val.at(idx);
		}
		const double& operator [](std::size_t idx) const {
			return val.at(idx);
		}

		// FIXME: should provide iterator interface, then all those function do not need to be marked as friend
		friend vector<N> operator+(vector<N> lhs, const vector<N>& rhs){
			std::transform(std::begin(lhs.val), std::end(lhs.val), std::begin(rhs.val), std::begin(lhs.val), std::plus<double>());
			return lhs;
		}
		friend vector<N> operator-(vector<N> lhs, const vector<N>& rhs){
			std::transform(std::begin(lhs.val), std::end(lhs.val), std::begin(rhs.val), std::begin(lhs.val), std::minus<double>());
			return lhs;
		}
		friend double operator*(const calcnum::vector<N>& lhs, const calcnum::vector<N>& rhs){
			// may not be optimal because of possible cancellation
			return std::inner_product(std::begin(lhs.val), std::end(lhs.val), std::begin(rhs.val), 0.0);
		}
		friend calcnum::vector<N> operator*(calcnum::vector<N> lhs, double rhs){
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
			if(rhs!=1){
#pragma GCC diagnostic pop
				std::transform(std::begin(lhs.val), std::end(lhs.val), std::begin(lhs.val), [rhs](double v){ return v*rhs;});
			}
			return lhs;
		}
		friend calcnum::vector<N> operator*(double lhs, calcnum::vector<N> rhs){
			return rhs*lhs;
		}
		friend calcnum::vector<N> operator/(calcnum::vector<N> lhs, double rhs){
			std::transform(std::begin(lhs.val), std::end(lhs.val), std::begin(lhs.val), [rhs](double v){ return v/rhs;});
			return lhs;
		}
	};

	template<std::size_t ROW, std::size_t COL>
	class matrix{
		std::array<std::array<double, COL>, ROW> val = {}; // 0-init
	public:
		matrix() = default;
		explicit matrix(const double (&val_)[ROW][COL]) {
			for(std::size_t row = 0; row != ROW; ++row){
				for(std::size_t col= 0; col != COL; ++col){
					val.at(row).at(col) = val_[row][col];
				}
			}
		}
		struct index{
			std::size_t row;
			std::size_t col;
		};
		double& operator [](index idx) {
			return val.at(idx.row).at(idx.col);
		}
		const double& operator [](index idx) const {
			return val.at(idx.row).at(idx.col);
		}
	};

	template<std::size_t ROW, std::size_t COL>
	calcnum::vector<COL> copy_row(const matrix<ROW,COL>& m, std::size_t row) {
		calcnum::vector<COL> v;
		for(std::size_t col = 0; col != COL; ++col){
			v[col] = m[{row,col}];
		}
		return v;
	}
	template<std::size_t ROW, std::size_t COL>
	calcnum::vector<ROW> copy_col(const matrix<ROW,COL>& m, std::size_t col) {
		calcnum::vector<ROW> v;
		for(std::size_t row = 0; row != ROW; ++row){
			v[row] = m[{row,col}];
		}
		return v;
	}
	template<std::size_t ROW, std::size_t COL>
	void set_row(matrix<ROW,COL>& m, std::size_t row, const calcnum::vector<COL> v) {
		for(std::size_t col = 0; col != COL; ++col){
			m[{row,col}] = v[col];
		}
	}
	template<std::size_t ROW, std::size_t COL>
	void set_col(matrix<ROW,COL>& m, std::size_t col, const calcnum::vector<ROW> v) {
		for(std::size_t row = 0; row != ROW; ++row){
			m[{row,col}] = v[row];
		}
	}

	template<std::size_t ROW, std::size_t COL>
	constexpr typename matrix<ROW,COL>::index msize(const matrix<ROW,COL>&) {
		return {ROW, COL};
	}

	template<std::size_t ROW, std::size_t COL>
	matrix<ROW, COL> identity(){
		matrix<ROW,COL> toreturn;
		for(std::size_t n = 0; n != std::min(COL,ROW); ++n){
			toreturn[{n,n}] = 1;
		}
		return toreturn;
	}
	template<std::size_t ROW, std::size_t COL>
	matrix<ROW, COL> operator+(matrix<ROW, COL> lhs, const matrix<ROW, COL>& rhs){
		for (std::size_t row=0;row != ROW;++row){
			for (std::size_t col = 0;col != COL; ++col){
				lhs[{row,col}] += rhs[{row,col}];
			}
		}
		return lhs;
	}
	template<std::size_t ROW, std::size_t COL>
	matrix<ROW, COL> operator-(matrix<ROW, COL> lhs, const matrix<ROW, COL>& rhs){
		for (std::size_t col = 0;col != COL; ++col){
			for (std::size_t row=0;row != ROW;++row){
				lhs[{row,col}] -= rhs[{row,col}];
			}
		}
		return lhs;
	}
	template<std::size_t ROW, std::size_t COL>
	matrix<ROW, COL> operator*(matrix<ROW, COL> lhs, double rhs){
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
		if(rhs!=1){
#pragma GCC diagnostic pop
			for (std::size_t row=0;row != ROW;++row){
				for (std::size_t col = 0;col != COL; ++col){
					lhs[{row,col}] *= rhs;
				}
			}
		}
		return lhs;
	}
	template<std::size_t ROW, std::size_t COL>
	matrix<ROW, COL> operator*(double lhs, matrix<ROW, COL> rhs){
		return rhs*lhs;
	}

	template<std::size_t N, std::size_t M, std::size_t P>
	matrix<N, P> operator*(const matrix<N, M>& lhs, const matrix<M, P>& rhs){
		// FIXME: use Coppersmith-Winograd or Strassen algorithm to yeld better complecity
		matrix<N, P> result={};
		for (std::size_t row=0; row != N;++row){
			for (std::size_t k = 0; k != M; ++k) {
				for (std::size_t col = 0; col != P; ++col){
					result[{row,col}] += lhs[{row,k}]*rhs[{k,col}];
				}
			}
		}
		return result;
	}

	template<std::size_t N>
	double norm_1(const vector<N>& x){
		double res = 0;
		for(std::size_t i = 0; i != N; ++i){
			res += std::fabs(x[i]);
		}
		return res;
	}
	template<std::size_t N>
	double norm_2(const vector<N>& x){
		return std::sqrt(x*x);
	}
	template<std::size_t N, std::size_t p>
	double norm_p(const vector<N>& x){
		static_assert(p>0, "0 norm does not exist");
		switch(p){
		case 1: return norm_1(x);
		case 2: return norm_2(x);
		default:;// continue in this function
		}
		double res = 0;
		for(std::size_t i = 0; i != N; ++i){
			res += std::pow(std::fabs(x[i]), p);
		}
		return std::pow(res, 1.0/p);
	}
	template<std::size_t N>
	double norm_inf(const calcnum::vector<N>& x){
		// could use max_element with custom comparator if we had an iterator...
		double res = 0;
		for(std::size_t i = 0; i != N; ++i){
			res = std::max(res, std::fabs(x[i]));
		}
		return res;
	}

	template<std::size_t ROW, std::size_t COL>
	calcnum::vector<ROW> operator*(const matrix<ROW, COL>& lhs, const calcnum::vector<COL>& rhs){
		calcnum::vector<ROW> res;
		for (std::size_t col = 0;col != COL; ++col){
			for (std::size_t row=0;row != ROW;++row){
				res[row] += lhs[{row,col}] * rhs[col];
			}
		}
		return res;
	}
	template<std::size_t ROW, std::size_t COL>
	calcnum::vector<COL> operator*(const calcnum::vector<ROW>& lhs, const matrix<ROW, COL>& rhs){
		calcnum::vector<COL> res;
		for (std::size_t row=0;row != ROW;++row){
			for (std::size_t col = 0;col != COL; ++col){
				res[row] += lhs[{row,col}] * rhs[col];
			}
		}
		return res;
	}

	template<std::size_t ROW, std::size_t COL>
	calcnum::vector<COL> operator*(const matrix<1, ROW>& lhs, const matrix<ROW, COL>& rhs){
		calcnum::vector<COL> res;
		for (std::size_t row = 0; row != ROW; ++row) {
			for (std::size_t col = 0; col != COL; ++col){
				res[col] += lhs[{0,row}]*rhs[{row,col}];
			}
		}
		return res;
	}

	template<std::size_t ROW, std::size_t COL>
	bool is_diag(const matrix<ROW, COL>& m, double err, std::size_t banda = 1){
		assert((banda%2==1) && "banda"); // only values like 1,3,5,.. are allowed
		const auto diff = (banda-1)/2;
		for(std::size_t r = 0; r != ROW; ++r){
			for(std::size_t c = 0; c != COL; ++c){
				const auto diff_rc = std::max(r,c)- std::min(r,c); // avoid overflow
				if( (diff_rc > diff) && !approx_equal(m[{r,c}], 0, err)){
					return false;
				}
			}
		}
		return true;
	}

	template<std::size_t ROW, std::size_t COL>
	bool is_triang_low(const matrix<ROW, COL>& m, double err){
		for(std::size_t r = 0; r != ROW; ++r){
			for(std::size_t c = std::min(COL,r+1); c != COL; ++c){
				if( !approx_equal(m[{r,c}], 0, err)){
					return false;
				}
			}
		}
		return true;
	}
	template<std::size_t ROW, std::size_t COL>
	bool is_triang_up(const matrix<ROW, COL>& m, double err){
		for(std::size_t r = 0; r != ROW; ++r){
			for(std::size_t c = 0; c != std::min(COL,r); ++c){
				if( !approx_equal(m[{r,c}], 0, err)){
					return false;
				}
			}
		}
		return true;
	}


	template<std::size_t N>
	double tr(const matrix<N, N>& m){
		double t = 0;
		for(std::size_t i = 0; i != N; ++i){
			t+=m[{i,i}];
		}
		return t;
	}

	inline double det(const matrix<1, 1>& m){
		return m[{0,0}];
	}

	inline double det(const matrix<2, 2>& m){
		return m[{0,0}]*m[{1,1}] - m[{0,1}]*m[{1,0}];
}

template<std::size_t ROW, std::size_t COL>
std::ostream& operator<<(std::ostream& os, const matrix<ROW, COL>& M){
	os << "{";
	for(size_t col = 0; col != COL; ++col){
		os << "{";
		for(size_t row = 0; row != ROW-1; ++row){
			os << M[{col,row}] << ", ";
		}
		os << M[{col,ROW-1}] << "},";
	}
	os << "}";
	return os;
}

// return B minus column and row pos
template<std::size_t N>
matrix<N-1, N-1> submat(const matrix<N, N>& B, const typename matrix<N, N>::index& pos){
	matrix<N-1, N-1> C;
	for(std::size_t col = 0, col2 = 0; col != N; ++col){
		if(col == pos.col){
			continue;
		}
		for(std::size_t row = 0, row2 = 0; row != N; ++row){
			if(row == pos.row){
				continue;
			}
			C[{row2,col2}] = B[{row,col}];
			++row2;
		}
		++col2;
	}
	return C;
}
template<std::size_t N>
double det_lap(const matrix<N, N>& B){
	// recursive, use Laplace expansion, prefer using the LU algorithm
	double determinant = 0;
	std::size_t row = 0; // fixed, can be chosen at runtime to stabilize alg
	for(std::size_t col = 0; col != N; ++col){
		matrix<N-1, N-1> C = submat(B,{row,col});
		auto detc = det(C);
		int sign = ((col+row) % 2 == 0) ? 1 : -1;
		determinant += sign * B[{row,col}]* detc;
	}
	return determinant;
}

template<std::size_t N>
matrix<N,N> pivotize_by_row(matrix<N, N>& A){
	auto A_ = A;
	auto P = identity<N,N>();
	for (std::size_t col = 0; col != N; ++col) {
		double a_max = A[{col,col}];
		std::size_t row = col;
		for(std::size_t j = col; j != N; ++j){ // column fixed, search biggest value on row
			if(A[{j,col}]>a_max){
				a_max = A[{j,col}];
				row = j;
			}
		}
		if(row != col){
			for(std::size_t j = 0; j != N; ++j){ // swap rows
				std::swap(A_[{col,j}], A_[{row,j}]);
				std::swap(P[{col,j}], P[{row,j}]);
			}
		}
	}
	A = A_;
	return P;
}

template<std::size_t N>
struct lu_fact_res{
	matrix<N, N> P;
	matrix<N, N> L;
	matrix<N, N> U;
};
template<std::size_t N>
lu_fact_res<N> LU(matrix<N, N>& A){
	const auto P = pivotize_by_row(A);
	// now apply LU
	matrix<N, N> L;
	matrix<N, N> U;

	for (std::size_t j = 0; j != N; ++j) {
		L[{j,j}] = 1;
		for (std::size_t i = 0; i != j + 1; ++i) {
			double s1 = 0;
			for (std::size_t k = 0; k != i; ++k){
				s1 += U[{k,j}] * L[{i,k}];
			}
			U[{i,j}] = A[{i,j}] - s1;
		}
		for (std::size_t i = j; i != N; ++i) {
			double s2 = 0;
			for (std::size_t k = 0; k != j; k++){
				s2 += U[{k,j}] * L[{i,k}];
			}
			L[{i,j}] = (A[{i,j}] - s2) / U[{j,j}];
		}
	}
	return {P, L, U};
}

// special impl for pivoting matrix
template<std::size_t N>
int det_pivot(const matrix<N, N>& B){
	// pivot is special, it should be doable without recursion
	auto detp = det_lap(B);
	return static_cast<int>(std::ceil(detp));
}

template<std::size_t N>
double det(const matrix<N, N>& B_){
	auto B = B_;
	auto lu_fact = LU(B);
	double toret = 1;
	for(std::size_t n = 0; n != N; ++n){
		toret *= lu_fact.U[{n,n}];
	}
	// P matrix has 1 or -1 as determinant
	return toret*det_lap(lu_fact.P);
}

template<std::size_t N>
calcnum::vector<N> solve_triang_inf(const matrix<N, N>& A, calcnum::vector<N> b){
	assert(is_triang_low(A, 0.001));
	auto& x = b;
	x[0] = b[0]/A[{0,0}];

	for(std::size_t i = 1; i != N; ++i){
		double sumprod = 0;
		for (std::size_t col = 0;col != i; ++col){
			sumprod += A[{i,col}] * x[col];
		}
		x[i] = (1/A[{i,i}]) * (b[i] - sumprod);
	}
	return x;
}

template<std::size_t N>
calcnum::vector<N> solve_triang_sup(const matrix<N, N>& A, calcnum::vector<N> b){
	assert(is_triang_up(A, 0.001));
	auto& x = b;
	x[N-1] = b[N-1]/A[{N-1,N-1}];

	for(std::size_t i_ = 1; i_ != N; ++i_){
		auto i = N-i_-1;
		double sumprod = 0;
		for (std::size_t col = i+1;col != N; ++col){
			sumprod += A[{i,col}] * x[col];
		}
		x[i] = (1/A[{i,i}]) * (b[i] - sumprod);
	}
	return x;
}


template<std::size_t N>
struct jacobi_status{
	calcnum::vector<N> x;
};

template<std::size_t N>
class jacobi_iter{
	matrix<N, N> BJ;
	calcnum::vector<N> db;
	jacobi_status<N> st;

public:
	explicit jacobi_iter(const matrix<N, N>& A_, const calcnum::vector<N>& b_, const calcnum::vector<N>& x0) : BJ(A_), db(b_){
		st.x = x0;
		// create iteration matrix
		for(std::size_t i = 0; i != N; ++i){
			for(std::size_t j = 0; j != N; ++j){
				if(i != j){
					BJ[{i,j}] = -BJ[{i,j}]/BJ[{i,i}];
				}
			}
			db[i] = db[i]/BJ[{i,i}];
			BJ[{i,i}] = 0;
		}
	}
	jacobi_iter& operator++(){
		st.x = BJ*st.x + db;
		return *this;
	}
	const jacobi_status<N>* operator->() const {
		return &st;
	}
};

template<std::size_t N>
struct gauss_seidel_status{
	calcnum::vector<N> x;
};

template<std::size_t N>
class gauss_seidel_iter{
	matrix<N, N> U;
	matrix<N, N> DL;
	calcnum::vector<N> db;
	gauss_seidel_status<N> st;

public:
	explicit gauss_seidel_iter(const matrix<N, N>& A, calcnum::vector<N> b_, calcnum::vector<N> x0){
		st.x = x0;
		// create iteration matrixes
		for(std::size_t i = 0; i != N; ++i){
			for(std::size_t j = 0; j != N; ++j){
				if(i < j){
					U[{i,j}] = A[{i,j}];
				}else{
					DL[{i,j}] = A[{i,j}];
				}
			}
		}
		db = solve_triang_inf(DL,b_);
	}
	gauss_seidel_iter& operator++(){
		auto y = U*st.x;
		st.x = db - solve_triang_inf(DL,y);
		return *this;
	}
	gauss_seidel_status<N>& operator*() {
		return st;
	}
	gauss_seidel_status<N>* operator->() {
		return &st;
	}
};

// for omega = 1, we have Gauss-Seidel
template<std::size_t N>
class SOR_iter{
	gauss_seidel_iter<N> it;
	double omega;

public:
	explicit SOR_iter(const matrix<N, N>& A, calcnum::vector<N> b_, calcnum::vector<N> x0, double omega_) :it(A,b_,x0), omega(omega_){
	}
	SOR_iter& operator++(){
		auto xk = it->x;
		++it;
		auto x = it->x;
		x = omega*x + (1-omega)*xk;
		it->x = x;
		return *this;
	}
	gauss_seidel_status<N>& operator*() {
		return it->st;
	}
	gauss_seidel_status<N>* operator->() {
		return &(*it);
	}
};


template<std::size_t N>
struct richardson_status{
	calcnum::vector<N> x;
};

// with constant gradient
template<std::size_t N>
class richardson_iter{
	matrix<N, N> BR;
	calcnum::vector<N> db;
	jacobi_status<N> st;

public:
	explicit richardson_iter(const matrix<N, N>& A, const calcnum::vector<N>& b_, const calcnum::vector<N>& x0, double alpha = 1) : db(b_){
		st.x = x0;
		// create iteration matrix
		BR = calcnum::identity<N,N> - alpha * A;
		db = db*alpha;
	}
	richardson_iter& operator++(){
		st.x = BR*st.x + db;
		return *this;
	}
	const jacobi_status<N>* operator->() const {
		return &st;
	}
};


template<std::size_t N>
struct power_status{
	calcnum::vector<N> z;
	double lambda;
};

template<std::size_t N>
class power_iter{
	matrix<N, N> A;
	power_status<N> st;

public:
	explicit power_iter(const matrix<N, N>& A_, calcnum::vector<N> z0) :A(A_){
		st.z = z0;
	}
	power_iter& operator++(){
		auto y = st.z * (1.0/norm_2(st.z));
		st.z = A*y;
		st.lambda = st.z*y;
		return *this;
	}
	power_status<N>& operator*() {
		return st;
	}
	power_status<N>* operator->() {
		return &st;
	}
};

}

