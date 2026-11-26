// 字串 -> n-gram -> bucket 計數向量
#ifndef STRING_MAPPER_H
#define STRING_MAPPER_H

#include <cstdint>
#include <string>
#include <vector>

class StringMapper {
public:
    // 將字串切成 n-gram，頭尾各補 (n-1) 個 '#'，共 |s|+n-1 個 gram
    static std::vector<std::string> extract_ngrams(const std::string& s, int n = 2);

    // 每個 n-gram 的 64-bit hash (固定演算法，跨平台結果一致)
    static std::vector<uint64_t> gram_hashes(const std::string& s, int n = 2);

    // 把 n-gram 雜湊到 L 個 bucket 並計數 (論文 Figure 1 的 gram counting vector)
    static std::vector<uint32_t> bucket_vector(const std::string& s, int n, int L);

    // 一次掃描同時算出「低維排序向量 (L 維)」與「高維 signature (S 維)」。
    // L 或 S 為 0 時略過該向量。
    static void two_resolution(const std::string& s, int n, int L, int S,
                               std::vector<uint32_t>& vecL, std::vector<uint32_t>& sig);

    static uint64_t hash_gram(const char* p, int n);
};

#endif // STRING_MAPPER_H
