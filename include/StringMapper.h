#ifndef STRING_MAPPER_H
#define STRING_MAPPER_H

#include <string>
#include <vector>

class StringMapper {
public:
    // 將字串切成 n-gram，預設 n=2 (2-gram)
    static std::vector<std::string> extract_ngrams(const std::string& s, int n = 2);
};

#endif // STRING_MAPPER_H