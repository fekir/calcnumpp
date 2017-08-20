#include <catch.hpp>

#include "intervals.hpp"

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"

TEST_CASE("closed_interval"){
    SECTION("case1"){
        const calcnum::closed_interval interval{0,1};
        REQUIRE(calcnum::length(interval)==1);
        REQUIRE(calcnum::midpoint(interval)==1./2);
    }
    SECTION("case2"){
        const calcnum::closed_interval interval{-1,1};
        REQUIRE(calcnum::length(interval)==(1.- (-1.)));
        REQUIRE(calcnum::midpoint(interval)==(1. + -1.)/2);
    }
}

TEST_CASE("open_interval"){
    SECTION("case1"){
        const calcnum::open_interval interval{0,1};
        REQUIRE(calcnum::length(interval)==1);
        REQUIRE(calcnum::midpoint(interval)==1./2);
    }
    SECTION("case2"){
        const calcnum::open_interval interval{-1,1};
        REQUIRE(calcnum::length(interval)==(1.- (-1.)));
        REQUIRE(calcnum::midpoint(interval)==(1. + -1.)/2);
    }
}
#pragma GCC diagnostic pop