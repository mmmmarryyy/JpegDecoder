#include <fstream>
#include <vector>
#include <numeric>

struct Node {
    size_t left_node_;
    size_t right_node_;
    uint8_t value_;
};

class HuffmanTree {
private:
    std::vector<Node> tree_;

    void BuildTree(size_t index, size_t depth, std::vector<uint8_t>& number_of_symbols,
                   std::vector<uint8_t> symbols, size_t& index_of_symbol) {
        if (depth > 0 && number_of_symbols[depth - 1] > 0) {
            tree_[index].value_ = symbols[index_of_symbol++];
            number_of_symbols[depth - 1]--;
            return;
        }

        if (index_of_symbol < symbols.size() && depth < number_of_symbols.size()) {
            tree_.emplace_back(0, 0, 0);
            tree_[index].left_node_ = tree_.size() - 1;
            BuildTree(tree_[index].left_node_, depth + 1, number_of_symbols, symbols,
                      index_of_symbol);
            if (index_of_symbol < symbols.size() && depth < number_of_symbols.size()) {
                tree_.emplace_back(0, 0, 0);
                tree_[index].right_node_ = tree_.size() - 1;
                BuildTree(tree_[index].right_node_, depth + 1, number_of_symbols, symbols,
                          index_of_symbol);
            }
        }
    }

public:
    std::vector<Node> Build(const std::vector<uint8_t>& number_of_symbols,
                            const std::vector<uint8_t>& symbols) {
        if (std::accumulate(number_of_symbols.begin(), number_of_symbols.end(), 0u) !=
            symbols.size()) {
            throw "Invalid huffman intro information";
        }

        tree_.emplace_back(0, 0, 0);

        size_t index_of_symbol = 0;
        std::vector<uint8_t> number_of_symbols_copy = number_of_symbols;

        BuildTree(0, 0, number_of_symbols_copy, symbols, index_of_symbol);

        for (auto i : number_of_symbols_copy) {
            if (i != 0) {
                throw "Bad huffman tree";
            }
        }

        return tree_;
    }
};