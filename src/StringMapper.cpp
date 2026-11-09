#include "StringMapper.h"

std::vector<std::string> StringMapper::extract_ngrams(const std::string& s, int n) {
    std::vector<std::string> ngrams;
    
    // 如果字串是空的，直接回傳
    if (s.empty()) return ngrams;

    // 頭尾補上 (n-1) 個 '#' 作為邊界符號，這是為了讓字首和字尾特徵也被捕捉
    std::string padded_s = std::string(n - 1, '#') + s + std::string(n - 1, '#');

    // 切出長度為 n 的滑動視窗 (Sliding Window)
    for (size_t i = 0; i <= padded_s.length() - n; ++i) {
        ngrams.push_back(padded_s.substr(i, n));
    }

    return ngrams;
}