#pragma once

#include <string_view>

void CheckImage(std::string_view filename, std::string_view comment = "");
void ExpectFail(std::string_view filename);
