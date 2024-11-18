#pragma once

#include "../image.h"

#include <jpeglib.h>
#include <filesystem>
#include <ranges>

inline Image ReadJpg(const std::filesystem::path& path) {
    jpeg_decompress_struct cinfo;
    jpeg_error_mgr err;
    auto* infile = std::fopen(path.c_str(), "rb");
    if (!infile) {
        throw std::runtime_error{"can't open " + path.string()};
    }

    cinfo.err = jpeg_std_error(&err);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);

    jpeg_read_header(&cinfo, true);
    jpeg_start_decompress(&cinfo);

    auto row_stride = cinfo.output_width * cinfo.output_components;
    auto buffer = (*cinfo.mem->alloc_sarray)(reinterpret_cast<j_common_ptr>(&cinfo), JPOOL_IMAGE,
                                             row_stride, 1);

    Image result{cinfo.output_width, cinfo.output_height};
    size_t y = 0;
    while (cinfo.output_scanline < cinfo.output_height) {
        jpeg_read_scanlines(&cinfo, buffer, 1);
        for (auto x : std::views::iota(0u, result.Width())) {
            RGB pixel;
            if (cinfo.output_components == 3) {
                pixel.r = buffer[0][x * 3];
                pixel.g = buffer[0][x * 3 + 1];
                pixel.b = buffer[0][x * 3 + 2];
            } else {
                pixel.r = pixel.g = pixel.b = buffer[0][x];
            }
            result.SetPixel(y, x, pixel);
        }
        ++y;
    }

    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    std::fclose(infile);
    return result;
}
