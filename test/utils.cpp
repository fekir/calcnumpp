#include <catch.hpp>

#include "utils.hpp"

TEST_CASE("approx_equal"){
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
