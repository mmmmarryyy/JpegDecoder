#pragma once

#include "image.h"

#include <filesystem>

Image Decode(const std::filesystem::path& path);
