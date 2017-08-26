#include <catch.hpp>

#include "utils.hpp"

TEST_CASE("approx_equal", "[utils][approx_equal]"){
    SECTION("good case 1"){
        REQUIRE(calcnum::approx_equal(1.1, 1.2, 0.2));
    }
    SECTION("good case 2"){
        REQUIRE(calcnum::approx_equal(1.2, 1.1, 0.2));
    }
    SECTION("good case 3"){
        REQUIRE(calcnum::approx_equal(-1.2, -1.1, 0.2));
    }
    SECTION("bad case"){
        REQUIRE_FALSE(calcnum::approx_equal(1.1, 1.2, 0.01));
    }
}

TEST_CASE("signum", "[utils][signum]"){
    REQUIRE(calcnum::signum( 1.) > 0);
    REQUIRE(calcnum::signum( 0.) == 0);
    REQUIRE(calcnum::signum(-1.) < 0);
}
