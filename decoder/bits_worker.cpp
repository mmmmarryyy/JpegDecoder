#include <fstream>
#include <vector>

class BitsWorker {
public:
    uint8_t BitParser(uint8_t byte, uint8_t start_bit, uint8_t end_bit) {
        uint8_t mask = ((1 << (end_bit - start_bit + 1)) - 1) << start_bit;
        return (byte & mask) >> start_bit;
    }

    void GetNewByte(std::istream& stream) {
        stream.read(reinterpret_cast<char*>(&current_byte_), 1);
        if (current_byte_ == 0xFF) {
            uint8_t test_byte;
            stream.read(reinterpret_cast<char*>(&test_byte), 1);

            if (test_byte != 0x00 && test_byte != 0xFF) {
                throw "unexpected byte";
            }
        }
    }

    uint8_t GetNewBit(std::istream& stream) {
        if (current_position_in_byte_ == 7) {
            GetNewByte(stream);
        }

        uint8_t bit =
            (current_byte_ & (1 << current_position_in_byte_)) >> current_position_in_byte_;

        if (current_position_in_byte_ == 0) {
            current_position_in_byte_ = 8;
        }

        current_position_in_byte_--;
        return bit;
    }

private:
    uint8_t current_byte_ = 0;
    uint8_t current_position_in_byte_ = 7;
};