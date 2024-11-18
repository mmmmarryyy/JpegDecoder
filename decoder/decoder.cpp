#include "../decoder.h"
#include "../image.h"
#include "fftw3.h"

#include "huffman.cpp"
#include "bits_worker.cpp"
#include "consts_and_structs.cpp"

#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <vector>
#include <memory>
#include <numeric>

class Decoder {
private:
    bool soi_flag_ = false;
    bool eoi_flag_ = false;
    bool sof0_flag_ = false;

    uint8_t number_of_components_ = 0;

    size_t image_width_, image_height_;

    BitsWorker bits_worker_;

    std::vector<SoF0Component> sof0_;
    std::vector<SosComponent> sos_;
    std::map<std::pair<uint8_t, uint8_t>, std::vector<Node>> dht_;
    std::map<uint8_t, std::vector<std::vector<uint16_t>>> dqt_;

    std::string comment_;

    void ParseSoF0(std::istream& i_stream) {
        uint8_t content_length[2];
        i_stream.read(reinterpret_cast<char*>(&content_length), 2);

        int length =
            (static_cast<int>(content_length[0]) * 256 + static_cast<int>(content_length[1]));
        length -= 2;

        i_stream.seekg(1, std::ios::cur);  // precision
        length--;

        uint8_t height[2];
        uint8_t width[2];

        i_stream.read(reinterpret_cast<char*>(&height), 2);
        i_stream.read(reinterpret_cast<char*>(&width), 2);

        length -= 4;

        if ((width[0] == 0 && width[1] == 0) || (height[0] == 0 && height[1] == 0)) {
            throw "Wrong image";
        }

        image_width_ = (static_cast<size_t>(width[0]) * 256 + static_cast<size_t>(width[1]));
        image_height_ = (static_cast<size_t>(height[0]) * 256 + static_cast<size_t>(height[1]));

        uint8_t inner_number_of_components;
        i_stream.read(reinterpret_cast<char*>(&inner_number_of_components), 1);
        length--;

        if (inner_number_of_components > 3 || inner_number_of_components == 0) {
            throw "wrong number of components";
        }

        if (number_of_components_ == 0) {
            number_of_components_ = inner_number_of_components;
        } else {
            if (number_of_components_ != inner_number_of_components) {
                throw "wrong number of components";
            }
        }

        for (int i = 0; i < number_of_components_; ++i) {
            SoF0Component sof0_component;

            i_stream.read(reinterpret_cast<char*>(&sof0_component.component_id), 1);
            length--;

            if (sof0_component.component_id > number_of_components_ ||
                sof0_component.component_id != i + 1) {
                throw "wrong component id";
            }

            uint8_t samplings_reader;
            i_stream.read(reinterpret_cast<char*>(&samplings_reader), 1);
            length--;

            sof0_component.vartical_sampl = bits_worker_.BitParser(samplings_reader, 0, 3);
            sof0_component.horizontal_sampl = bits_worker_.BitParser(samplings_reader, 4, 7);

            i_stream.read(reinterpret_cast<char*>(&sof0_component.qtable_id), 1);
            length--;

            sof0_.push_back(sof0_component);
        }

        if (length != 0) {
            throw "read wrong number of bytes";
        } else {
            return;
        }
    }

    void ParseDht(std::istream& i_stream) {
        uint8_t content_length[2];
        i_stream.read(reinterpret_cast<char*>(&content_length), 2);

        int length =
            (static_cast<int>(content_length[0]) * 256 + static_cast<int>(content_length[1]));
        length -= 2;

        while (length > 0) {
            std::vector<Node> huffman_tree;
            uint8_t number_of_ht, type_of_ht, ht_imformation;
            std::vector<uint8_t> number_of_symbols, symbols;

            i_stream.read(reinterpret_cast<char*>(&ht_imformation), 1);
            length--;

            number_of_ht = bits_worker_.BitParser(ht_imformation, 0, 3);
            type_of_ht = bits_worker_.BitParser(ht_imformation, 4, 4);

            if (bits_worker_.BitParser(ht_imformation, 5, 7) != 0) {
                throw "wrong value for not_used parameter";
            }

            int sum = 0;
            for (int i = 0; i < 16; ++i) {
                uint8_t counter;
                i_stream.read(reinterpret_cast<char*>(&counter), 1);
                length--;

                number_of_symbols.push_back(counter);
                sum += counter;
            }

            while (number_of_symbols.back() == 0 && !number_of_symbols.empty()) {
                number_of_symbols.pop_back();
            }

            for (int i = 0; i < sum; ++i) {
                uint8_t symbol;
                i_stream.read(reinterpret_cast<char*>(&symbol), 1);
                length--;

                symbols.push_back(symbol);
            }

            std::pair<uint8_t, uint8_t> key = std::make_pair(number_of_ht, type_of_ht);
            if (dht_.find(key) != dht_.end()) {
                throw "try to add huffman table with same parameters";
            }

            HuffmanTree tree_builder;
            huffman_tree = tree_builder.Build(number_of_symbols, symbols);
            dht_[key] = huffman_tree;
        }

        if (length != 0) {
            throw "read wrong number of bytes";
        } else {
            return;
        }
    }

    void ParseDqt(std::istream& i_stream) {
        uint8_t content_length[2];
        i_stream.read(reinterpret_cast<char*>(&content_length), 2);

        int length =
            (static_cast<int>(content_length[0]) * 256 + static_cast<int>(content_length[1]));
        length -= 2;

        while (length > 0) {
            uint8_t qt_information_byte;
            i_stream.read(reinterpret_cast<char*>(&qt_information_byte), 1);
            length--;

            uint8_t qtable_id;
            qtable_id = bits_worker_.BitParser(qt_information_byte, 0, 3);

            if (qtable_id > 3) {
                throw "bad qtable id";
            }

            uint8_t precision = bits_worker_.BitParser(qt_information_byte, 4, 7) == 0 ? 1 : 2;
            int data_length = 64 * precision;

            std::vector<uint8_t> qt_table_data(data_length);
            std::vector<std::vector<uint16_t>> qt_values(8, std::vector<uint16_t>(8));

            i_stream.read(reinterpret_cast<char*>(&qt_table_data[0]), data_length);
            length -= data_length;

            for (int i = 0; i < data_length / precision; ++i) {
                auto zig_zag_index = zig_zag_order[i];

                if (precision == 1) {
                    qt_values[zig_zag_index / 8][zig_zag_index % 8] =
                        static_cast<uint16_t>(qt_table_data[i]);
                } else {
                    qt_values[zig_zag_index / 8][zig_zag_index % 8] =
                        static_cast<uint16_t>(qt_table_data[2 * i]) * 256 +
                        static_cast<uint16_t>(qt_table_data[2 * i + 1]);
                }
            }

            dqt_[qtable_id] = qt_values;
        }

        if (length != 0) {
            throw "read wrong number of bytes";
        } else {
            return;
        }
    }

    void ParseSos(std::istream& i_stream) {
        uint8_t content_length[2];
        i_stream.read(reinterpret_cast<char*>(&content_length), 2);

        int length =
            (static_cast<int>(content_length[0]) * 256 + static_cast<int>(content_length[1]));
        length -= 2;

        uint8_t inner_number_of_components;
        i_stream.read(reinterpret_cast<char*>(&inner_number_of_components), 1);
        length--;

        if (inner_number_of_components > 3 || inner_number_of_components == 0) {
            throw "wrong number of components";
        }

        if (number_of_components_ == 0) {
            number_of_components_ = inner_number_of_components;
        } else {
            if (number_of_components_ != inner_number_of_components) {
                throw "wrong number of components";
            }
        }

        for (uint8_t i = 0; i < inner_number_of_components; ++i) {
            SosComponent component;
            i_stream.read(reinterpret_cast<char*>(&component.component_id), 1);
            length--;

            if (component.component_id > inner_number_of_components) {
                throw "bad component id";
            }

            uint8_t haffman_tree_id;
            i_stream.read(reinterpret_cast<char*>(&haffman_tree_id), 1);
            length--;

            component.haffman_table_AC = bits_worker_.BitParser(haffman_tree_id, 0, 3);
            component.haffman_table_DC = bits_worker_.BitParser(haffman_tree_id, 4, 7);
            sos_.push_back(component);
        }

        uint8_t first_skip_byte, second_skip_byte, third_skip_byte;

        i_stream.read(reinterpret_cast<char*>(&first_skip_byte), 1);
        length--;

        i_stream.read(reinterpret_cast<char*>(&second_skip_byte), 1);
        length--;

        i_stream.read(reinterpret_cast<char*>(&third_skip_byte), 1);
        length--;

        if (first_skip_byte != 0x00 || second_skip_byte != 0x3F || third_skip_byte != 0x00) {
            throw "wronh skip bytes for baseline";
        }

        if (length != 0) {
            throw "read wrong number of bytes";
        } else {
            return;
        }
    }

    void ParseCom(std::istream& i_stream) {
        uint8_t content_length[2];
        i_stream.read(reinterpret_cast<char*>(&content_length), 2);
        int comment_part_length =
            (static_cast<int>(content_length[0]) * 256 + static_cast<int>(content_length[1])) - 2;

        std::string comment_part;
        comment_part.resize(comment_part_length);
        i_stream.read(reinterpret_cast<char*>(&comment_part[0]), comment_part_length);
        comment_ += comment_part;
    }

    uint8_t GetMaxSampling(uint8_t type) {
        uint8_t max_sampling = 0;

        for (auto component : sof0_) {
            if (type == 0) {
                max_sampling = std::max(component.horizontal_sampl, max_sampling);
            } else {
                max_sampling = std::max(component.vartical_sampl, max_sampling);
            }
        }

        return max_sampling;
    }

    uint8_t GetHaffmanTreeValue(std::istream& i_stream, std::vector<Node> huffman_tree) {
        auto current_node = huffman_tree[0];
        int counter = 0;
        while (counter++ < 16) {
            if (bits_worker_.GetNewBit(i_stream) == 1) {
                if (current_node.right_node_ == 0) {
                    throw "wrong huffman code";
                }

                current_node = huffman_tree[current_node.right_node_];
            } else {
                if (current_node.left_node_ == 0) {
                    throw "wrong huffman code";
                }

                current_node = huffman_tree[current_node.left_node_];
            }

            if (current_node.left_node_ == 0 && current_node.right_node_ == 0) {
                return current_node.value_;
            }
        }

        throw "wrong huffman code: length > 16";
    }

    int16_t GetCoefByHuffmanTreeValue(std::istream& i_stream, uint8_t n) {
        if (n == 0) {
            return 0;
        } else {
            int16_t result = bits_worker_.GetNewBit(i_stream);
            bool riest_bit_value = (result == 0);

            for (int i = 1; i < n; ++i) {
                result *= 2;
                result += bits_worker_.GetNewBit(i_stream);
            }

            if (riest_bit_value) {
                result = result - (1 << n) + 1;
            }

            return result;
        }
    }

    std::vector<std::vector<int16_t>> ParseTable(std::istream& i_stream, uint8_t table_id_dc,
                                                 uint8_t table_id_ac) {
        std::vector<int16_t> table;

        if (dht_.find({table_id_dc, 0}) == dht_.end() ||
            dht_.find({table_id_ac, 1}) == dht_.end()) {
            throw "wrong id of huffman table";
        }

        uint8_t huffman_tree_value = GetHaffmanTreeValue(i_stream, dht_[{table_id_dc, 0}]);
        if (huffman_tree_value == 0) {
            table.push_back(0);
        } else {
            table.push_back(GetCoefByHuffmanTreeValue(i_stream, huffman_tree_value));
        }

        while (table.size() < 64) {
            uint8_t huffman_tree_value = GetHaffmanTreeValue(i_stream, dht_[{table_id_ac, 1}]);

            if (huffman_tree_value == 0) {
                while (table.size() < 64) {
                    table.push_back(0);
                }
                break;
            } else {
                int number_of_zeros = bits_worker_.BitParser(huffman_tree_value, 4, 7);

                if (table.size() + number_of_zeros >= 64) {
                    throw "wrong number of zeroes: table won't be 64 in size";
                }

                for (int _ = 0; _ < number_of_zeros; ++_) {
                    table.push_back(0);
                }

                table.push_back(GetCoefByHuffmanTreeValue(
                    i_stream, bits_worker_.BitParser(huffman_tree_value, 0, 3)));
            }
        }

        std::vector<std::vector<int16_t>> result(8, std::vector<int16_t>(8));

        for (size_t i = 0; i < table.size(); ++i) {
            result[zig_zag_order[i] / 8][zig_zag_order[i] % 8] = table[i];
        }

        return result;
    }

    std::vector<std::vector<std::vector<std::vector<std::vector<int16_t>>>>> ParseTables(
        std::istream& i_stream) {

        int16_t max_h = GetMaxSampling(0);
        int16_t max_v = GetMaxSampling(1);
        std::vector<int16_t> previous_dc_coefs(3, 0);

        std::vector<std::vector<std::vector<std::vector<std::vector<int16_t>>>>> tables(3);

        uint32_t number_of_mcu_blocks_h = (image_width_ + 8 * max_h - 1) / (8 * max_h);
        uint32_t number_of_mcu_blocks_w = (image_height_ + 8 * max_v - 1) / (8 * max_v);

        for (uint16_t mcu_blocks_v_id = 0; mcu_blocks_v_id < number_of_mcu_blocks_w;
             ++mcu_blocks_v_id) {
            for (uint16_t mcu_blocks_h_id = 0; mcu_blocks_h_id < number_of_mcu_blocks_h;
                 ++mcu_blocks_h_id) {
                for (SoF0Component component : sof0_) {
                    tables[component.component_id - 1].resize(
                        number_of_mcu_blocks_w * component.vartical_sampl,
                        std::vector<std::vector<std::vector<int16_t>>>(number_of_mcu_blocks_h *
                                                                       component.horizontal_sampl));
                    for (size_t table_in_mcu_block_v_id = 0;
                         table_in_mcu_block_v_id < component.vartical_sampl;
                         table_in_mcu_block_v_id++) {
                        for (size_t table_in_mcu_block_h_id = 0;
                             table_in_mcu_block_h_id < component.horizontal_sampl;
                             table_in_mcu_block_h_id++) {
                            for (SosComponent sos_component : sos_) {
                                if (component.component_id == sos_component.component_id) {
                                    std::vector<std::vector<int16_t>> table =
                                        ParseTable(i_stream, sos_component.haffman_table_DC,
                                                   sos_component.haffman_table_AC);
                                    previous_dc_coefs[component.component_id - 1] += table[0][0];
                                    table[0][0] = previous_dc_coefs[component.component_id - 1];

                                    if (dqt_.find(component.qtable_id) == dqt_.end()) {
                                        throw "wrong dqt_ id";
                                    }

                                    std::vector<std::vector<uint16_t>> qt =
                                        dqt_[component.qtable_id];

                                    for (size_t i = 0; i < 8; ++i) {
                                        for (size_t j = 0; j < 8; ++j) {
                                            table[i][j] *= qt[i][j];
                                        }
                                    }

                                    tables[component.component_id - 1]
                                          [mcu_blocks_v_id * component.vartical_sampl +
                                           table_in_mcu_block_v_id]
                                          [mcu_blocks_h_id * component.horizontal_sampl +
                                           table_in_mcu_block_h_id] = table;
                                }
                            }
                        }
                    }
                }
            }
        }

        return tables;
    };

    RGB YCbCrToRgb(double y, double cb, double cr) {
        int r = static_cast<int>(round(1.402 * (cr - 128.0) + y));
        int g = static_cast<int>(round(y - 0.344136 * (cb - 128.0) - 0.714136 * (cr - 128.0)));
        int b = static_cast<int>(round(1.772 * (cb - 128.0) + y));

        r = std::min(std::max(0, r), 255);
        g = std::min(std::max(0, g), 255);
        b = std::min(std::max(0, b), 255);

        return RGB{r, g, b};
    }

    Image Parse(std::istream& i_stream) {
        if (!sof0_flag_) {
            throw "don't have necessary data";
        }

        std::vector<std::vector<std::vector<std::vector<std::vector<int16_t>>>>> tables =
            ParseTables(i_stream);

        int16_t max_h = GetMaxSampling(0);
        int16_t max_v = GetMaxSampling(1);
        uint32_t number_of_mcu_blocks_h = (image_width_ + 8 * max_h - 1) / (8 * max_h);
        uint32_t number_of_mcu_blocks_w = (image_height_ + 8 * max_v - 1) / (8 * max_v);

        Image image(image_width_, image_height_);
        image.SetComment(comment_);

        double data[64];
        double result[64];
        for (uint16_t mcu_blocks_v_id = 0; mcu_blocks_v_id < number_of_mcu_blocks_w;
             ++mcu_blocks_v_id) {
            for (uint16_t mcu_blocks_h_id = 0; mcu_blocks_h_id < number_of_mcu_blocks_h;
                 ++mcu_blocks_h_id) {
                std::vector<std::vector<int16_t>> local_tables[3];
                for (SoF0Component component : sof0_) {
                    local_tables[component.component_id - 1].resize(
                        8 * component.vartical_sampl,
                        std::vector<int16_t>(8 * component.horizontal_sampl));

                    for (size_t table_in_mcu_block_v_id = 0;
                         table_in_mcu_block_v_id < component.vartical_sampl;
                         table_in_mcu_block_v_id++) {
                        for (size_t table_in_mcu_block_h_id = 0;
                             table_in_mcu_block_h_id < component.horizontal_sampl;
                             table_in_mcu_block_h_id++) {

                            for (size_t i = 0; i < 8; ++i) {
                                for (size_t j = 0; j < 8; ++j) {
                                    data[i * 8 + j] =
                                        tables[component.component_id - 1]
                                              [mcu_blocks_v_id * component.vartical_sampl +
                                               table_in_mcu_block_v_id]
                                              [mcu_blocks_h_id * component.horizontal_sampl +
                                               table_in_mcu_block_h_id][i][j];
                                }
                            }

                            for (size_t i = 0; i < 8; ++i) {
                                data[i] *= std::sqrt(2);
                                data[i * 8] *= std::sqrt(2);
                            }

                            for (size_t i = 0; i < 64; ++i) {
                                data[i] /= 16.0;
                            }

                            auto* execution_plan = fftw_plan_r2r_2d(
                                8, 8, data, result, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE);
                            fftw_execute(execution_plan);
                            fftw_destroy_plan(execution_plan);

                            for (size_t i = 0; i < 8; ++i) {
                                for (size_t j = 0; j < 8; ++j) {
                                    local_tables[component.component_id -
                                                 1][8 * table_in_mcu_block_v_id +
                                                    i][8 * table_in_mcu_block_h_id + j] =
                                        std::min<int16_t>(
                                            255, std::max<int16_t>(
                                                     0, 128 + std::round(result[i * 8 + j])));
                                }
                            }
                        }
                    }
                }

                for (int16_t i = 0; i < 8 * max_v; ++i) {
                    for (int16_t j = 0; j < 8 * max_h; ++j) {
                        int16_t y_cb_cr[3] = {0, 0, 0};

                        for (auto component : sof0_) {
                            y_cb_cr[component.component_id - 1] =
                                local_tables[component.component_id - 1]
                                            [i * component.vartical_sampl / max_v]
                                            [j * component.horizontal_sampl / max_h];
                        }
                        const auto rgb = number_of_components_ == 1
                                             ? RGB{y_cb_cr[0], y_cb_cr[0], y_cb_cr[0]}
                                             : YCbCrToRgb(y_cb_cr[0], y_cb_cr[1], y_cb_cr[2]);
                        size_t y_index = mcu_blocks_v_id * 8 * max_v + i;
                        size_t x_index = mcu_blocks_h_id * 8 * max_h + j;
                        if (y_index < image_height_ && x_index < image_width_) {
                            image.SetPixel(y_index, x_index, rgb);
                        }
                    }
                }
            }
        }
        return image;
    }

public:
    Image DecodeJpefFile(std::istream& i_stream) {  // TODO: convert to std::istream&
        Image result_image;

        uint8_t marker[2];
        i_stream.read(reinterpret_cast<char*>(&marker[0]), 2);
        if (marker[0] != 0xFF || marker[1] != 0xD8) {
            throw "don't have soi at the beginning of i_stream";
        }

        soi_flag_ = true;

        while (!i_stream.eof()) {
            uint8_t first_byte;
            uint8_t second_byte;
            i_stream.read(reinterpret_cast<char*>(&first_byte), 1);
            if (first_byte != 0xFF) {
                throw "first_byte not 0xFF";
            }

            i_stream.read(reinterpret_cast<char*>(&second_byte), 1);
            if (second_byte == 0xD8) {
                if (soi_flag_) {
                    throw "already have soi";
                }

                soi_flag_ = true;
            } else if (second_byte == 0xC0) {
                if (sof0_flag_) {
                    throw "already have sof0_";
                }

                sof0_flag_ = true;
                ParseSoF0(i_stream);
            } else if (second_byte == 0xC4) {
                ParseDht(i_stream);
            } else if (second_byte == 0xDB) {
                ParseDqt(i_stream);
            } else if ((second_byte & 0xF0) == 0xE0) {
                uint8_t content_length[2];
                i_stream.read(reinterpret_cast<char*>(&content_length), 2);
                int data_length = (static_cast<int>(content_length[0]) * 256 +
                                   static_cast<int>(content_length[1])) -
                                  2;

                std::vector<uint8_t> data(data_length);
                i_stream.read(reinterpret_cast<char*>(&data[0]), data_length);
            } else if (second_byte == 0xFE) {
                ParseCom(i_stream);
            } else if (second_byte == 0xDA) {
                ParseSos(i_stream);
                result_image = Parse(i_stream);

                uint8_t checker;
                i_stream.read(reinterpret_cast<char*>(&checker), 1);
                if (checker != 0xFF) {
                    throw "don't have eoi after sos_ part";
                }
                i_stream.read(reinterpret_cast<char*>(&checker), 1);
                if (checker != 0xD9) {
                    throw "don't have eoi after sos_ part";
                }

                eoi_flag_ = true;
                break;
            } else if (second_byte == 0xD9) {
                eoi_flag_ = true;
                break;
            } else {
                throw "wrong header";
            }
        }
        if (!eoi_flag_) {
            throw "don't have eoi part at the end of i_stream";
        }

        return result_image;
    }
};

Image Decode(const std::filesystem::path& path) {
    Decoder decoder;
    std::ifstream file(path, std::ios::binary);
    if (!file.is_open()) {
        throw "Error opening file";
    }
    std::istream& i_stream = file;
    auto result_image = decoder.DecodeJpefFile(i_stream);
    file.close();
    return result_image;
};
