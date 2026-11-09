#include <iostream>
#include "StringMapper.h"

int main() {
    std::string test_str = "Jim Gray";
    
    std::cout << "Original String: \"" << test_str << "\"\n";
    std::cout << "Extracting 2-grams...\n";

    // 呼叫我們剛寫好的模組
    auto ngrams = StringMapper::extract_ngrams(test_str, 2);

    // 印出結果看看
    std::cout << "Result: [ ";
    for (const auto& gram : ngrams) {
        std::cout << "'" << gram << "' ";
    }
    std::cout << "]\n";
    std::cout << "Total n-grams: " << ngrams.size() << "\n";

    return 0;
}