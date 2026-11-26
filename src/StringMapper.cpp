#include "StringMapper.h"

static inline uint64_t mix64(uint64_t x) {
    x ^= x >> 33;
    x *= 0xff51afd7ed558ccdULL;
    x ^= x >> 33;
    x *= 0xc4ceb9fe1a85ec53ULL;
    x ^= x >> 33;
    return x;
}

uint64_t StringMapper::hash_gram(const char* p, int n) {
    uint64_t h = 1469598103934665603ULL;  // FNV-1a
    for (int i = 0; i < n; ++i) {
        h ^= (unsigned char)p[i];
        h *= 1099511628211ULL;
    }
    return mix64(h);
}

static std::string pad(const std::string& s, int n) {
    return std::string(n - 1, '#') + s + std::string(n - 1, '#');
}

std::vector<std::string> StringMapper::extract_ngrams(const std::string& s, int n) {
    std::vector<std::string> ngrams;
    if (s.empty() || n < 1) return ngrams;
    std::string p = pad(s, n);
    for (size_t i = 0; i + n <= p.size(); ++i) ngrams.push_back(p.substr(i, n));
    return ngrams;
}

std::vector<uint64_t> StringMapper::gram_hashes(const std::string& s, int n) {
    std::vector<uint64_t> hs;
    if (s.empty() || n < 1) return hs;
    std::string p = pad(s, n);
    for (size_t i = 0; i + n <= p.size(); ++i) hs.push_back(hash_gram(p.data() + i, n));
    return hs;
}

std::vector<uint32_t> StringMapper::bucket_vector(const std::string& s, int n, int L) {
    std::vector<uint32_t> v(L, 0), unused;
    two_resolution(s, n, L, 0, v, unused);
    return v;
}

void StringMapper::two_resolution(const std::string& s, int n, int L, int S,
                                  std::vector<uint32_t>& vecL, std::vector<uint32_t>& sig) {
    vecL.assign(L > 0 ? L : 0, 0);
    sig.assign(S > 0 ? S : 0, 0);
    if (s.empty()) return;
    std::string p = pad(s, n);
    for (size_t i = 0; i + n <= p.size(); ++i) {
        uint64_t h = hash_gram(p.data() + i, n);
        if (L > 0) vecL[h % (uint64_t)L]++;
        if (S > 0) sig[h % (uint64_t)S]++;
    }
}
