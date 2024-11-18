#include "test_commons.h"

#include <ranges>

#include <catch2/catch_test_macros.hpp>

TEST_CASE("small jfif (4:2:0)") {
    CheckImage("small.jpg", ":)");
}

TEST_CASE("jfif (4:4:4)") {
    CheckImage("lenna.jpg");
}

TEST_CASE("jfif (4:2:2)") {
    CheckImage("bad_quality.jpg", "so quality");
}

TEST_CASE("tiny jfif (4:4:4)") {
    CheckImage("tiny.jpg");
}

TEST_CASE("exif (4:2:2)") {
    CheckImage("chroma_halfed.jpg");
}

TEST_CASE("exif (grayscale)") {
    CheckImage("grayscale.jpg");
}

TEST_CASE("jfif/exif (4:2:0)") {
    CheckImage("test.jpg");
}

TEST_CASE("exif (4:4:4)") {
    CheckImage("colors.jpg");
}

TEST_CASE("photoshop (4:4:4)") {
    CheckImage("save_for_web.jpg");
}

TEST_CASE("Error handling") {
    for (auto i : std::views::iota(1, 25)) {
        ExpectFail("bad" + std::to_string(i) + ".jpg");
    }
}
