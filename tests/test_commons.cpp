#include "png.h"
#include "jpg.h"
#include "util.h"
#include "../decoder.h"

#include <cmath>

#include <catch2/catch_test_macros.hpp>

namespace {

double Distance(const RGB& lhs, const RGB& rhs) {
    auto square = [](int x) { return x * x; };
    auto r_sqr = square(lhs.r - rhs.r);
    auto g_sqr = square(lhs.g - rhs.g);
    auto b_sqr = square(lhs.b - rhs.b);
    return std::sqrt(r_sqr + g_sqr + b_sqr);
}

void Compare(const Image& actual, const Image& expected) {
    REQUIRE(actual.Height() == expected.Height());
    REQUIRE(actual.Width() == expected.Width());

    auto sum = 0.;
    for (auto y : std::views::iota(0u, actual.Height())) {
        for (auto x : std::views::iota(0u, actual.Width())) {
            auto actual_pixel = actual.GetPixel(y, x);
            auto expected_pixel = expected.GetPixel(y, x);
            sum += Distance(actual_pixel, expected_pixel);
        }
    }
    auto mean = sum / (actual.Width() * actual.Height());
    REQUIRE(mean < 5);
}

const auto kImagesPath = GetFileDir(__FILE__) / "images";

}  // namespace

void CheckImage(std::string_view filename, std::string_view comment = "") {
    auto png_filename = std::filesystem::path{filename}.replace_extension("png");
    auto image = Decode(kImagesPath / filename);
    CHECK(image.GetComment() == comment);
    WritePng(kImagesPath / png_filename, image);
    auto ok_image = ReadJpg(kImagesPath / filename);
    Compare(image, ok_image);
}

void ExpectFail(std::string_view filename) {
    INFO(filename);
    auto path = kImagesPath / "bad" / filename;
    if (!std::filesystem::exists(path)) {
        FAIL(path << " does not exist");
    }
    CHECK_THROWS(Decode(path));
}
