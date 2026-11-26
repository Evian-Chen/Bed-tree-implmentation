// Z-order (Morton) 編碼：把 L 維計數向量交錯成一個一維 key
// 以及：由 key 區間推出每一維計數的上下界 (論文 Figure 3)，並計算編輯距離下界 (Eq. 4)
#ifndef BED_ZORDER_H
#define BED_ZORDER_H

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace bed {

const int kKeyWords = 4;  // key 最多 256 bit

struct ZKey {
    uint64_t w[kKeyWords];  // w[0] 是最高位字組
    ZKey() { for (int i = 0; i < kKeyWords; ++i) w[i] = 0; }
    bool operator<(const ZKey& o) const {
        for (int i = 0; i < kKeyWords; ++i)
            if (w[i] != o.w[i]) return w[i] < o.w[i];
        return false;
    }
    bool operator==(const ZKey& o) const {
        for (int i = 0; i < kKeyWords; ++i) if (w[i] != o.w[i]) return false;
        return true;
    }
    bool operator!=(const ZKey& o) const { return !(*this == o); }
    // g: 由最高位算起的 bit 編號 (0-based)
    bool bit(int g) const { return (w[g >> 6] >> (63 - (g & 63))) & 1ULL; }
    void set_bit(int g, bool v) {
        uint64_t m = 1ULL << (63 - (g & 63));
        if (v) w[g >> 6] |= m; else w[g >> 6] &= ~m;
    }
    std::string to_bit_string(int T) const {
        std::string s;
        for (int g = 0; g < T; ++g) s += bit(g) ? '1' : '0';
        return s;
    }
    static ZKey from_bit_string(const std::string& s) {
        ZKey k;
        for (size_t g = 0; g < s.size() && g < 256; ++g) k.set_bit((int)g, s[g] == '1');
        return k;
    }
};

struct ZOrder {
    int L;          // 維度數 (bucket 數)
    int B;          // 每維 bit 數；計數超過 maxv 時飽和
    int T;          // 總 bit 數 = L*B
    uint32_t maxv;  // 2^B - 1，等於 maxv 代表「>= maxv」(上界未知)

    ZOrder(int L_ = 4, int B_ = 8) : L(L_), B(B_), T(L_ * B_) {
        if (L < 1 || B < 1 || B > 16 || T > 64 * kKeyWords)
            throw std::invalid_argument("ZOrder: need 1<=B<=16 and L*B<=256");
        maxv = (B >= 32) ? 0xffffffffu : ((1u << B) - 1u);
    }

    // 預設每維 bit 數：最多 8，且 L*B 不超過 256
    static int default_bits(int L) {
        int b = (64 * kKeyWords) / L;
        return b > 8 ? 8 : (b < 1 ? 1 : b);
    }

    // 由最高位開始，依序輪流取各維的下一個 bit
    ZKey encode(const uint32_t* v) const {
        ZKey k;
        for (int j = 0; j < B; ++j)
            for (int d = 0; d < L; ++d) {
                uint32_t x = v[d] > maxv ? maxv : v[d];  // 飽和
                k.set_bit(j * L + d, (x >> (B - 1 - j)) & 1u);
            }
        return k;
    }
    ZKey encode(const std::vector<uint32_t>& v) const { return encode(v.data()); }

    void decode(const ZKey& k, uint32_t* out) const {
        for (int d = 0; d < L; ++d) out[d] = 0;
        for (int j = 0; j < B; ++j)
            for (int d = 0; d < L; ++d)
                out[d] = (out[d] << 1) | (k.bit(j * L + d) ? 1u : 0u);
    }

    static int common_prefix(const ZKey& a, const ZKey& b, int T) {
        int g = 0;
        while (g < T && a.bit(g) == b.bit(g)) ++g;
        return g;
    }

    // key 區間 [lo, hi] 內任何 key 都共用前 p 個 bit，其餘 bit 未知：
    // 未知 bit 全 0 得到下界 lb，全 1 得到上界 ub。
    void bounds(const ZKey& lo, const ZKey& hi, uint32_t* lb, uint32_t* ub) const {
        int p = common_prefix(lo, hi, T);
        for (int d = 0; d < L; ++d) { lb[d] = 0; ub[d] = 0; }
        for (int j = 0; j < B; ++j)
            for (int d = 0; d < L; ++d) {
                int g = j * L + d;
                uint32_t bl, bu;
                if (g < p) { bl = bu = lo.bit(g) ? 1u : 0u; }
                else { bl = 0; bu = 1; }
                lb[d] = (lb[d] << 1) | bl;
                ub[d] = (ub[d] << 1) | bu;
            }
    }

    // 論文 Eq.(4)：編輯距離 >= ceil( max( Σ_{q>ub}(q-ub), Σ_{lb>q}(lb-q) ) / n )
    // 其中 ub == maxv 視為「無上界」，該維在第一項貢獻 0，確保飽和時仍是合法下界。
    int lower_bound_ed(const uint32_t* vq, const uint32_t* lb, const uint32_t* ub, int n) const {
        uint32_t s1 = 0, s2 = 0;
        for (int d = 0; d < L; ++d) {
            if (ub[d] != maxv && vq[d] > ub[d]) s1 += vq[d] - ub[d];
            if (lb[d] > vq[d]) s2 += lb[d] - vq[d];
        }
        uint32_t m = s1 > s2 ? s1 : s2;
        return (int)((m + n - 1) / n);
    }
};

}  // namespace bed

#endif
