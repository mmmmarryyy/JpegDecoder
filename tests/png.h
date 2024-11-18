#pragma once

#include "../image.h"

#include <png.h>
#include <ranges>
#include <filesystem>

inline void WritePng(const std::filesystem::path& path, const Image& image) {
    auto* fp = std::fopen(path.c_str(), "wb");
    if (!fp) {
        throw std::runtime_error{"Can't open file for writing " + path.string()};
    }
    auto png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    auto info = png_create_info_struct(png);
    if (setjmp(png_jmpbuf(png))) {
        throw std::runtime_error{"Shit happens"};
    }

    png_init_io(png, fp);
    png_set_IHDR(png, info, image.Width(), image.Height(), 8, PNG_COLOR_TYPE_RGBA,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);
    png_write_info(png, info);

    std::vector<png_bytep> bytes(image.Height());
    for (size_t y : std::views::iota(0u, image.Height())) {
        bytes[y] = static_cast<png_bytep>(malloc(png_get_rowbytes(png, info)));
        for (size_t x : std::views::iota(0u, image.Width())) {
            auto pixel = image.GetPixel(y, x);
            bytes[y][x * 4] = pixel.r;
            bytes[y][x * 4 + 1] = pixel.g;
            bytes[y][x * 4 + 2] = pixel.b;
            bytes[y][x * 4 + 3] = 255;
        }
    }
    png_write_image(png, bytes.data());
    png_write_end(png, nullptr);

    std::fclose(fp);
    png_destroy_write_struct(&png, &info);
    for (auto i : std::views::iota(0u, image.Height())) {
        free(bytes[i]);
    }
}
