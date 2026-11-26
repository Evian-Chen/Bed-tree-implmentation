// zorder_tool：單獨測試 z-order 編碼、解碼、區間上下界、編輯距離下界
//   zorder_tool encode --L 4 --B 2 3 2 1 3            -> 11011011   (論文 Figure 2 範例)
//   zorder_tool decode --L 4 --B 2 11011011           -> 3 2 1 3
//   zorder_tool bounds --L 4 --B 2 11011011 11011110  -> 每維 [lb,ub]  (論文 Figure 3 範例)
//   zorder_tool lb --q "Jam Gray" --s1 "Jim Gray" --s2 "Jim Grey" --L 4 --n 2
//        以 s1,s2 的 key 為區間端點，算 query 到該區間的編輯距離下界，並列出區間內字串的真實距離
#include <algorithm>
#include <iostream>

#include "Cli.h"
#include "EditDistance.h"
#include "StringMapper.h"
#include "ZOrder.h"

using namespace bed;

int main(int argc, char** argv) {
    Args a(argc, argv);
    if (a.pos.empty()) {
        std::cerr << "usage: zorder_tool encode|decode|bounds|lb ... (see header comment)\n";
        return 2;
    }
    const std::string cmd = a.pos[0];
    int L = a.geti("L", 4);
    int B = a.geti("B", ZOrder::default_bits(L));

    if (cmd == "encode") {
        std::vector<uint32_t> v;
        for (size_t i = 1; i < a.pos.size(); ++i) v.push_back((uint32_t)std::atoi(a.pos[i].c_str()));
        L = (int)v.size();
        if (!a.has("B")) B = ZOrder::default_bits(L);
        ZOrder z(L, B);
        std::cout << z.encode(v).to_bit_string(z.T) << "\n";
    } else if (cmd == "decode" && a.pos.size() >= 2) {
        int T = (int)a.pos[1].size();
        if (!a.has("B")) B = T / L;
        ZOrder z(L, B);
        std::vector<uint32_t> v(L);
        z.decode(ZKey::from_bit_string(a.pos[1]), v.data());
        for (int d = 0; d < L; ++d) std::cout << v[d] << (d + 1 < L ? " " : "\n");
    } else if (cmd == "bounds" && a.pos.size() >= 3) {
        int T = (int)a.pos[1].size();
        if (!a.has("B")) B = T / L;
        ZOrder z(L, B);
        ZKey lo = ZKey::from_bit_string(a.pos[1]), hi = ZKey::from_bit_string(a.pos[2]);
        std::vector<uint32_t> lb(L), ub(L);
        z.bounds(lo, hi, lb.data(), ub.data());
        std::cout << "common prefix length = " << ZOrder::common_prefix(lo, hi, z.T) << " of " << z.T << " bits\n";
        for (int d = 0; d < L; ++d)
            std::cout << "  bucket " << (d + 1) << ": [" << lb[d] << ", " << ub[d] << "]"
                      << (ub[d] == z.maxv ? "  (ub saturated = unbounded)" : "") << "\n";
    } else if (cmd == "lb") {
        int n = a.geti("n", 2);
        std::string q = a.get("q"), s1 = a.get("s1"), s2 = a.get("s2");
        ZOrder z(L, B);
        ZKey k1 = z.encode(StringMapper::bucket_vector(s1, n, L));
        ZKey k2 = z.encode(StringMapper::bucket_vector(s2, n, L));
        if (k2 < k1) { std::swap(k1, k2); std::swap(s1, s2); }
        std::vector<uint32_t> lb(L), ub(L);
        z.bounds(k1, k2, lb.data(), ub.data());
        std::vector<uint32_t> vq = StringMapper::bucket_vector(q, n, L);
        int lbv = z.lower_bound_ed(vq.data(), lb.data(), ub.data(), n);
        std::cout << "interval: [" << k1.to_bit_string(z.T) << " (" << s1 << "), "
                  << k2.to_bit_string(z.T) << " (" << s2 << ")]\n";
        std::cout << "per-bucket bounds:";
        for (int d = 0; d < L; ++d) std::cout << " [" << lb[d] << "," << ub[d] << "]";
        std::cout << "\nquery vector     :";
        for (int d = 0; d < L; ++d) std::cout << " " << vq[d];
        std::cout << "\nlower bound LB(q, interval) = " << lbv << "\n";
        std::cout << "true ED(q,s1) = " << edit_distance(q, s1) << ", ED(q,s2) = " << edit_distance(q, s2)
                  << "   (LB must be <= every true ED of strings inside the interval)\n";
    } else {
        std::cerr << "bad command or arguments\n";
        return 2;
    }
    return 0;
}
