//          Copyright Federico Kircheis 2017
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE.txt or copy at
//          http://www.boost.org/LICENSE_1_0.txt)

#include <catch.hpp>

#include "matrix.hpp"
#include "utils.hpp"

#include <random>

namespace{

	calcnum::matrix<4,4> create_rand_matrix (unsigned int seed){
		std::mt19937 gen(seed);
		std::uniform_real_distribution<> dis(-10, 10);
		auto drand= [&dis, &gen](){
			return dis(gen);
		};
		calcnum::matrix<4,4> m({
		                           {drand(), drand(), drand(), drand()},
		                           {drand(), drand(), drand(), drand()},
		                           {drand(), drand(), drand(), drand()},
		                           {drand(), drand(), drand(), drand()}
		                       });
		return m;
	}
}

TEST_CASE("matrix operators"){
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
	SECTION("assignment"){
		calcnum::matrix<2, 2> I;
		I[{0,0}] = 1;
		REQUIRE(1 == Approx(I[{0,0}]));
	}
	
	SECTION("constructor"){
		const calcnum::matrix<2, 3> M({{1,2,3},{4,5,6}});
		REQUIRE((M[{0,0}]) == 1);
		REQUIRE((M[{0,1}]) == 2);
		REQUIRE((M[{0,2}]) == 3);
		REQUIRE((M[{1,0}]) == 4);
		REQUIRE((M[{1,1}]) == 5);
		REQUIRE((M[{1,2}]) == 6);
	}

	SECTION("identity"){
		auto I = calcnum::identity<2,3>();
		for(std::size_t n = 0; n != 2; ++n){
			REQUIRE((I[{n,n}]) == 1);
		}
	}

	calcnum::matrix<3, 2> M;
	M[{0,0}] = 1;
	M[{0,1}] = 2;
	M[{1,0}] = 3;
	M[{1,1}] = 4;
	M[{2,0}] = 5;
	M[{2,1}] = 6;
	SECTION("sum"){
		auto M2 = M + M;

		REQUIRE((M2[{0,0}]) == (M[{0,0}]+M[{0,0}]));
		REQUIRE((M2[{0,1}]) == (M[{0,1}]+M[{0,1}]));
		REQUIRE((M2[{1,0}]) == (M[{1,0}]+M[{1,0}]));
		REQUIRE((M2[{1,1}]) == (M[{1,1}]+M[{1,1}]));
		REQUIRE((M2[{2,0}]) == (M[{2,0}]+M[{2,0}]));
		REQUIRE((M2[{2,1}]) == (M[{2,1}]+M[{2,1}]));
	}


	SECTION("product-scalar"){
		auto M2 = 2*M;

		REQUIRE((M2[{0,0}]) == 2*(M[{0,0}]));
		REQUIRE((M2[{0,1}]) == 2*(M[{0,1}]));
		REQUIRE((M2[{1,0}]) == 2*(M[{1,0}]));
		REQUIRE((M2[{1,1}]) == 2*(M[{1,1}]));
		REQUIRE((M2[{2,0}]) == 2*(M[{2,0}]));
		REQUIRE((M2[{2,1}]) == 2*(M[{2,1}]));
	}
#pragma GCC diagnostic pop
	SECTION("product"){
		calcnum::matrix<2, 4> M2;
		M2[{0,0}] = 1;
		M2[{0,1}] = 2;
		M2[{0,2}] = 3;
		M2[{0,3}] = 4;
		M2[{1,0}] = 5;
		M2[{1,1}] = 6;
		M2[{1,2}] = 7;
		M2[{1,3}] = 8;

		auto M3 = M*M2;

		REQUIRE(Approx(M3[{0,0}]) == 11);
		REQUIRE(Approx(M3[{0,1}]) == 14);
		REQUIRE(Approx(M3[{0,2}]) == 17);
		REQUIRE(Approx(M3[{0,3}]) == 20);
		REQUIRE(Approx(M3[{1,0}]) == 23);
		REQUIRE(Approx(M3[{1,1}]) == 30);
		REQUIRE(Approx(M3[{1,2}]) == 37);
		REQUIRE(Approx(M3[{1,3}]) == 44);
		REQUIRE(Approx(M3[{2,0}]) == 35);
		REQUIRE(Approx(M3[{2,1}]) == 46);
		REQUIRE(Approx(M3[{2,2}]) == 57);
		REQUIRE(Approx(M3[{2,3}]) == 68);
	}
}

TEST_CASE("matrix properties"){

	SECTION("diagn"){
		SECTION("diag"){
			SECTION("good"){
				calcnum::matrix<3, 2> m1;
				m1[{0,0}] = 1;
				m1[{1,1}] = 4;
				REQUIRE(calcnum::is_diag(m1, 0.0001));
			}
			SECTION("bad"){
				calcnum::matrix<3, 2> m1;
				m1[{0,0}] = 1;
				m1[{1,1}] = 4;
				m1[{2,1}] = 4;
				REQUIRE_FALSE(calcnum::is_diag(m1, 0.0001));
			}
		}
		SECTION("tridiag"){
			SECTION("good"){
				calcnum::matrix<5, 4> m1;
				for(std::size_t n = 0; n != 4; ++n){
					m1[{n,n}] = 1;
					m1[{n+1,n}] = 1;
				}
				REQUIRE(calcnum::is_diag(m1, 0.0001,3));
			}
			SECTION("bad"){
				calcnum::matrix<5, 4> m1;
				for(std::size_t n = 0; n != 4; ++n){
					m1[{n,n}] = 1;
					m1[{n+1,n}] = 1;
				}
				m1[{0,2}] = 1;
				REQUIRE_FALSE(calcnum::is_diag(m1, 0.0001, 1));
			}
		}
	}

	SECTION("triang"){
		SECTION("inf"){
			SECTION("good"){
				calcnum::matrix<3, 2> m1;
				m1[{0,0}] = 1;
				m1[{1,1}] = 4;
				m1[{2,1}] = 4;
				REQUIRE(calcnum::is_triang_low(m1, 0.0001));
			}
			SECTION("bad"){
				calcnum::matrix<3, 2> m1;
				m1[{0,0}] = 1;
				m1[{1,1}] = 4;
				m1[{0,1}] = 4;
				REQUIRE_FALSE(calcnum::is_triang_low(m1, 0.0001));
			}
		}
		SECTION("sup"){
			SECTION("good"){
				calcnum::matrix<3, 2> m1;
				m1[{0,0}] = 1;
				m1[{0,1}] = 1;
				m1[{1,1}] = 1;
				REQUIRE(calcnum::is_triang_up(m1, 0.0001));
			}
			SECTION("bad"){
				calcnum::matrix<3, 2> m1;
				m1[{0,0}] = 1;
				m1[{1,1}] = 4;
				m1[{2,1}] = 4;
				REQUIRE_FALSE(calcnum::is_triang_up(m1, 0.0001));
			}
		}

	}
}

TEST_CASE("LU"){
	SECTION("pivotize"){

		calcnum::matrix<2, 2> m1;
		m1[{0,0}] = 1;
		m1[{0,1}] = 4;
		m1[{1,0}] = 2;
		m1[{1,1}] = 3;

		auto p = calcnum::pivotize_by_row(m1);
		REQUIRE(Approx(p[{0,0}]) == 0);
		REQUIRE(Approx(p[{0,1}]) == 1);
		REQUIRE(Approx(p[{1,0}]) == 1);
		REQUIRE(Approx(p[{1,1}]) == 0);
	}
	SECTION("LU (random)"){
		std::random_device rd;
		auto seed = rd();
		auto m1 = create_rand_matrix(seed);
		INFO("seed: " << seed);
		const auto A = m1;
		auto res = calcnum::LU(m1);
		{
			INFO("L: " << res.L);
			INFO("U: " << res.U);
			REQUIRE(calcnum::is_triang_low(res.L, 0.0001));
			REQUIRE(calcnum::is_triang_up(res.U, 0.0001));
		}
		INFO("A: " <<A);
		INFO("PA (res): " << m1);
		INFO("P: " << res.P);
		INFO("PA: " << res.P*A);
		INFO("LU: " << res.L*res.U);

		for(size_t i = 0; i != msize(A).col; ++i){
			for(size_t j = 0; j != msize(A).row; ++j){
				REQUIRE(Approx((res.P*A)[{i,j}]) == (m1[{i,j}]));
				REQUIRE(Approx((res.L*res.U)[{i,j}]) == ((res.P*A)[{i,j}]));
			}
		}
	}
}

TEST_CASE("determinant"){
	SECTION("determinant"){
		SECTION("1x1"){
			calcnum::matrix<1, 1> p;
			p[{0,0}]=1;
			auto det = calcnum::det(p);
			REQUIRE(Approx(det) == 1);
		}
		SECTION("2x2"){
			calcnum::matrix<2, 2> p;
			p[{0,1}]=1;
			p[{1,0}]=1;

			auto det = calcnum::det(p);
			REQUIRE(Approx(det) == -1);

		}
		SECTION("pivot"){
			calcnum::matrix<3, 3> p;
			SECTION("t1"){
				p[{0,1}]=1;
				p[{1,0}]=1;
				p[{2,2}] = 1;

				auto det = calcnum::det_pivot(p);
				REQUIRE(det==-1);
			}
			SECTION("t2"){
				p[{0,0}] = 1;
				p[{1,1}] = 1;
				p[{2,2}] = 1;

				auto det = calcnum::det_pivot(p);
				REQUIRE(det==1);
			}
		}
		SECTION(">2x2"){
			calcnum::matrix<3, 3> m1;
			m1[{0,0}] = 1;
			m1[{0,1}] = 3;
			m1[{0,2}] = 5;
			m1[{1,0}] = 2;
			m1[{1,1}] = 4;
			m1[{1,2}] = 7;
			m1[{2,0}] = 1;
			m1[{2,1}] = 1;
			m1[{2,2}] = 0;
			auto det = calcnum::det(m1);
			REQUIRE(Approx(det) == 4);

			auto det2 = calcnum::det_lap(m1);
			REQUIRE(Approx(det2) == det);
		}
		SECTION(">2x2, singular"){
			calcnum::matrix<3, 3> m1;
			m1[{0,0}] = 1;
			m1[{0,1}] = 2;
			m1[{0,2}] = 3;
			m1[{1,0}] = 4;
			m1[{1,1}] = 5;
			m1[{1,2}] = 6;
			m1[{2,0}] = 7;
			m1[{2,1}] = 8;
			m1[{2,2}] = 9;
			auto det = calcnum::det(m1);
			REQUIRE(Approx(det) == 0);

			auto det2 = calcnum::det_lap(m1);
			REQUIRE(Approx(det2) == det);

		}
	}
}

TEST_CASE("solve system"){
	SECTION("solve system exactly"){
		SECTION("triang sup"){
			constexpr std::size_t N = 3;
			calcnum::matrix<N,N> A;
			A[{0,0}] =  5;
			A[{0,1}] =  4;
			A[{0,2}] = -1;

			A[{1,1}] = 10;
			A[{1,2}] = -3;

			A[{2,2}] = 1;

			calcnum::vector<N> b({0,11,3});

			auto x = calcnum::solve_triang_sup(A,b);

			REQUIRE(Approx(x[0]) ==  -1);
			REQUIRE(Approx(x[1]) ==  2);
			REQUIRE(Approx(x[2]) ==  3);
		}
		SECTION("triang inf"){
			constexpr std::size_t N = 3;
			calcnum::matrix<N,N> A;
			A[{0,0}] = 1;

			A[{1,0}] = -3;
			A[{1,1}] = 10;

			A[{2,0}] = -1;
			A[{2,1}] =  4;
			A[{2,2}] =  5;

			calcnum::vector<N> b({3,11,0});

			auto x = calcnum::solve_triang_inf(A,b);

			REQUIRE(Approx(x[0]) ==  3);
			REQUIRE(Approx(x[1]) ==  2);
			REQUIRE(Approx(x[2]) == -1);
		}
		SECTION("jacobi"){
			constexpr std::size_t N = 3;
			calcnum::matrix<N,N> A;
			A[{0,0}] = 3;
			A[{0,1}] = 1;
			A[{0,2}] = 1;

			A[{1,0}] = 2;
			A[{1,1}] = 3;
			A[{1,2}] = 0;

			A[{2,0}] =-1;
			A[{2,1}] = 2;
			A[{2,2}] = 4;

			calcnum::vector<N> b({5,5,5});
			calcnum::vector<N> sol({1,1,1});
			calcnum::jacobi_iter<N> it(A, b, {});
			auto err = std::numeric_limits<double>::infinity();
			while(err > 0.0001){
				++it;
				auto newerr = calcnum::norm_inf(sol - it->x);
				REQUIRE(newerr < err);
				err= newerr;
			}
		}
		SECTION("gauss-seider"){
			constexpr std::size_t N = 3;
			calcnum::matrix<N,N> A;
			A[{0,0}] = 3;
			A[{0,1}] = 0;
			A[{0,2}] = 4;

			A[{1,0}] = 7;
			A[{1,1}] = 4;
			A[{1,2}] = 2;

			A[{2,0}] = 1;
			A[{2,1}] = 1;
			A[{2,2}] = 2;

			calcnum::vector<N> b({7,13,4});
			calcnum::vector<N> sol({1,1,1});
			calcnum::gauss_seidel_iter<N> it(A, b, {});
			auto err = std::numeric_limits<double>::infinity();
			while(err > 0.0001){
				++it;
				auto newerr = calcnum::norm_inf(sol - it->x);
				REQUIRE(newerr < err);
				err= newerr;
			}
		}
		SECTION("SOR"){
			constexpr std::size_t N = 3;
			calcnum::matrix<N,N> A;
			A[{0,0}] = 3;
			A[{0,1}] = 0;
			A[{0,2}] = 4;

			A[{1,0}] = 7;
			A[{1,1}] = 4;
			A[{1,2}] = 2;

			A[{2,0}] = 1;
			A[{2,1}] = 1;
			A[{2,2}] = 2;

			calcnum::vector<N> b({7,13,4});
			calcnum::vector<N> sol({1,1,1});
			calcnum::SOR_iter<N> it(A, b, {}, 1); // with omega = 1 it's the same as GS
			auto err = std::numeric_limits<double>::infinity();
			while(err > 0.0001){
				++it;
				auto newerr = calcnum::norm_inf(sol - it->x);
				REQUIRE(newerr < err);
				err= newerr;
			}
		}
		SECTION("eigenvalues"){
			constexpr std::size_t N = 2;
			calcnum::matrix<N,N> A;
			A[{0,0}] =  2;
			A[{0,1}] = -4;

			A[{1,0}] = -1;
			A[{1,1}] = -1;

			calcnum::vector<N> y0({0,1});

			calcnum::power_iter<N> it(A, y0);
			const double val = 3;
			double err = std::numeric_limits<double>::infinity();
			double max_err = 0.00001;
			for(int num_it = 0; (num_it != 100) && (err > max_err); ++num_it){
				++it;
				err = std::fabs(val - it->lambda);
			}
			REQUIRE(err < max_err);
			auto eigv = it->z * (1/calcnum::norm_2(it->z));
			REQUIRE(Approx(std::fabs(eigv[0])) == 4.0/std::sqrt(17));
			REQUIRE(Approx(std::fabs(eigv[1])) == 1.0/std::sqrt(17));
		}
	}

}

