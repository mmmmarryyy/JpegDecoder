#pragma once

#include <vector>
#include <string>

struct RGB {
    int r, g, b;
};

class Image {
public:
    Image() = default;
    Image(size_t width, size_t height) {
        SetSize(width, height);
    }

    void SetSize(size_t width, size_t height) {
        data_.assign(height, std::vector<RGB>(width));
    }

    size_t Width() const {
        return data_.empty() ? 0 : data_[0].size();
    }

    size_t Height() const {
        return data_.size();
    }

    void SetPixel(size_t y, size_t x, const RGB& pixel) {
        data_[y][x] = pixel;
    }

    RGB GetPixel(size_t y, size_t x) const {
        return data_[y][x];
    }

    RGB& GetPixel(size_t y, size_t x) {
        return data_[y][x];
    }

    void SetComment(std::string comment) {
        comment_ = std::move(comment);
    }

    const std::string& GetComment() const {
        return comment_;
    }

private:
    std::vector<std::vector<RGB>> data_;
    std::string comment_;
};
