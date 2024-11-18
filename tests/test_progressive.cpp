#include "test_commons.h"

#include <catch2/catch_test_macros.hpp>

TEST_CASE("jfif (grayscale)") {
    CheckImage("progressive.jpg");
}

TEST_CASE("jfif (4:4:4)") {
    CheckImage("progressive_small.jpg");
}

TEST_CASE("jfif/exif (4:2:2)") {
    CheckImage("progressive_2.jpg", "such decoder");
}
