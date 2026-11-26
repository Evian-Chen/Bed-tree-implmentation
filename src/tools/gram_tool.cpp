// gram_tool：單獨測試 n-gram 切分、bucket 向量、z-order key、signature
//   gram_tool "Jim Gray"                         預設 n=2, L=4
//   gram_tool "Jim Gray" --n 2 --L 4 --sig 64    同時顯示 64 維 signature
#include <iostream>

#include "Cli.h"
#include "StringMapper.h"
#include "ZOrder.h"

using namespace bed;

int main(int argc, char** argv) {
    Args a(argc, argv);
    if (a.pos.empty()) {
        std::cerr << "usage: gram_tool <string> [--n 2] [--L 4] [--sig 64] [--B bits]\n";
        return 2;
    }
    const std::string s = a.pos[0];
    int n = a.geti("n", 2), L = a.geti("L", 4), S = a.geti("sig", 0);
    int B = a.geti("B", ZOrder::default_bits(L));

    std::vector<std::string> grams = StringMapper::extract_ngrams(s, n);
    std::vector<uint64_t> hs = StringMapper::gram_hashes(s, n);
    std::cout << "string = \"" << s << "\"   n = " << n << "\n";
    std::cout << "n-grams (" << grams.size() << " = |s|+n-1 = " << (s.size() + n - 1) << "):\n";
    for (size_t i = 0; i < grams.size(); ++i) {
        std::cout << "  \"" << grams[i] << "\"  hash%" << L << " = " << (hs[i] % (uint64_t)L);
        if (S > 0) std::cout << "   hash%" << S << " = " << (hs[i] % (uint64_t)S);
        std::cout << "\n";
    }

    std::vector<uint32_t> v, sig;
    StringMapper::two_resolution(s, n, L, S, v, sig);
    std::cout << "\nL=" << L << " gram counting vector = <";
    for (int d = 0; d < L; ++d) std::cout << v[d] << (d + 1 < L ? ", " : "");
    std::cout << ">\n";

    ZOrder z(L, B);
    ZKey k = z.encode(v);
    std::cout << "z-order key (" << z.T << " bits, " << B << " bits per dim, saturates at " << z.maxv
              << ") = " << k.to_bit_string(z.T) << "\n";

    if (S > 0) {
        std::cout << "\nsignature (S=" << S << " buckets), non-zero entries (bucket:count):\n  ";
        for (int d = 0; d < S; ++d) if (sig[d]) std::cout << d << ":" << sig[d] << " ";
        uint32_t tot = 0;
        for (int d = 0; d < S; ++d) tot += sig[d];
        std::cout << "\n  sum = " << tot << " (should equal number of n-grams)\n";
    }
    return 0;
}
