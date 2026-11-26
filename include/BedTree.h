// Bed-tree (gram counting order) 與 Bed+-tree (two-resolution gram signatures)
//
//   cfg.sig == 0  : 原版 Bed-tree。key = L 維 gram 計數向量的 z-order。
//   cfg.sig  > 0  : Bed+-tree。排序仍用 L 維 (低維)，但每個 leaf entry 額外存 sig 維
//                   (高維) 的計數 signature，進 DP 之前先用它算 Eq.(4) 下界過濾。
// 兩者共用同一份程式，差別只在是否存/使用 signature，比較才公平。
#ifndef BED_BEDTREE_H
#define BED_BEDTREE_H

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#include "BPlusTree.h"
#include "EditDistance.h"
#include "StringMapper.h"
#include "ZOrder.h"

namespace bed {

const int kMaxSig = 64;

struct BedConfig {
    int n = 2;                  // gram 長度
    int L = 4;                  // 排序用的 bucket 數 (低維)
    int B = 0;                  // z-order 每維 bit 數；0 = 自動
    int sig = 0;                // signature bucket 數；0 = 不使用 (原版 Bed-tree)
    size_t page_bytes = 4096;   // 模擬的 page 大小，用來決定節點容量
    size_t string_bytes = 24;   // 估計 leaf 內每個字串佔的 bytes (只影響節點容量)
    double fill = 0.9;          // bulk load 填充率
};

struct Entry {
    ZKey key;
    uint32_t id;
    uint8_t has_sig;            // 0 表示該字串 gram 數 > 255，signature 飽和不可用
    uint8_t sig[kMaxSig];
    Entry() : id(0), has_sig(0) { std::memset(sig, 0, sizeof(sig)); }
};

struct QueryStats {
    uint64_t internal_visited = 0;  // 讀到的內部節點數 (≈ page reads)
    uint64_t leaf_visited = 0;      // 讀到的 leaf 數 (≈ page reads)
    uint64_t children_pruned = 0;   // 被內部節點下界剪掉的子樹數
    uint64_t entries_scanned = 0;   // 在 leaf 內被檢查的 entry 數
    uint64_t pruned_by_key = 0;     // 被「L 維 key 向量」下界濾掉
    uint64_t pruned_by_sig = 0;     // 被「signature」下界濾掉 (Bed+ 才有)
    uint64_t len_skipped = 0;       // 長度差 > tau 直接略過
    uint64_t dp_calls = 0;          // 真正跑 DP 驗證的次數
    uint64_t results = 0;
    uint64_t nodes() const { return internal_visited + leaf_visited; }
    void add(const QueryStats& o) {
        internal_visited += o.internal_visited; leaf_visited += o.leaf_visited;
        children_pruned += o.children_pruned; entries_scanned += o.entries_scanned;
        pruned_by_key += o.pruned_by_key; pruned_by_sig += o.pruned_by_sig;
        len_skipped += o.len_skipped; dp_calls += o.dp_calls; results += o.results;
    }
};

struct StorageInfo {
    size_t leaves, internals, height, entries;
    size_t leaf_entry_bytes, leaf_capacity, internal_capacity;
    size_t total_pages() const { return leaves + internals; }
};

// signature 下界：Eq.(4) 套用在 sig 維向量，除以 n 後無條件進位
inline int sig_lower_bound(const uint32_t* sq, const uint8_t* sv, int S, int n) {
    uint32_t s1 = 0, s2 = 0;
    for (int d = 0; d < S; ++d) {
        uint32_t a = sq[d], b = sv[d];
        if (a > b) s1 += a - b; else s2 += b - a;
    }
    uint32_t m = s1 > s2 ? s1 : s2;
    return (int)((m + n - 1) / n);
}

class BedTree {
public:
    explicit BedTree(const BedConfig& cfg)
        : cfg_(fix(cfg)), z_(cfg_.L, cfg_.B), tree_(0, 0) {
        if (cfg_.sig < 0 || cfg_.sig > kMaxSig) throw std::invalid_argument("sig must be 0..64");
        size_t key_bytes = (size_t)(z_.T + 7) / 8;
        leaf_entry_bytes_ = key_bytes + 4 + cfg_.string_bytes + (size_t)cfg_.sig;
        size_t internal_entry_bytes = 2 * key_bytes + 8;
        tree_ = BPlusTree<Entry>(std::max<size_t>(2, cfg_.page_bytes / leaf_entry_bytes_),
                                 std::max<size_t>(3, cfg_.page_bytes / internal_entry_bytes));
    }

    const BedConfig& config() const { return cfg_; }
    const ZOrder& zorder() const { return z_; }
    size_t size() const { return tree_.size(); }
    const BPlusTree<Entry>& tree() const { return tree_; }

    // 批次建樹 (排序後 bulk load)
    void build(const std::vector<std::string>& data) {
        strs_ = data;
        alive_.assign(data.size(), 1);
        std::vector<Entry> es;
        es.reserve(data.size());
        for (size_t i = 0; i < data.size(); ++i) es.push_back(make_entry((uint32_t)i, data[i]));
        std::sort(es.begin(), es.end(), [](const Entry& a, const Entry& b) {
            return BPlusTree<Entry>::kid(a) < BPlusTree<Entry>::kid(b);
        });
        tree_.bulk_load(es, cfg_.fill);
    }

    uint32_t insert(const std::string& s) {
        uint32_t id = (uint32_t)strs_.size();
        strs_.push_back(s);
        alive_.push_back(1);
        tree_.insert(make_entry(id, s));
        return id;
    }

    bool erase(uint32_t id) {
        if (id >= strs_.size() || !alive_[id]) return false;
        Entry e = make_entry(id, strs_[id]);
        if (!tree_.erase(KeyId(e.key, id))) return false;
        alive_[id] = 0;
        return true;
    }

    const std::string& str(uint32_t id) const { return strs_[id]; }
    bool alive(uint32_t id) const { return id < alive_.size() && alive_[id]; }
    size_t total_ids() const { return strs_.size(); }

    // 範圍查詢：回傳所有 ED(q, s) <= tau 的字串 id (遞增排序)
    std::vector<uint32_t> range_query(const std::string& q, int tau, QueryStats* st = nullptr) const {
        QueryStats local;
        QueryStats& S = st ? *st : local;
        Ctx c;
        c.q = &q;
        c.tau = tau;
        StringMapper::two_resolution(q, cfg_.n, cfg_.L, cfg_.sig, c.vq, c.sq);
        c.lb.resize(cfg_.L);
        c.ub.resize(cfg_.L);
        c.pt.resize(cfg_.L);
        std::vector<uint32_t> out;
        search(tree_.root(), c, out, S);
        std::sort(out.begin(), out.end());
        S.results += out.size();
        return out;
    }

    // 對照組：不用索引，逐一驗證 (length filter + DP)
    std::vector<uint32_t> scan_query(const std::string& q, int tau, QueryStats* st = nullptr) const {
        QueryStats local;
        QueryStats& S = st ? *st : local;
        std::vector<uint32_t> out;
        for (size_t i = 0; i < strs_.size(); ++i) {
            if (!alive_[i]) continue;
            ++S.entries_scanned;
            if (std::abs((int)strs_[i].size() - (int)q.size()) > tau) { ++S.len_skipped; continue; }
            ++S.dp_calls;
            if (verify_ed(q, strs_[i], tau) <= tau) out.push_back((uint32_t)i);
        }
        S.results += out.size();
        return out;
    }

    StorageInfo storage() const {
        StorageInfo s;
        tree_.count_nodes(&s.leaves, &s.internals);
        s.height = (size_t)tree_.height();
        s.entries = tree_.size();
        s.leaf_entry_bytes = leaf_entry_bytes_;
        s.leaf_capacity = tree_.leaf_capacity();
        s.internal_capacity = tree_.internal_capacity();
        return s;
    }

    Entry make_entry(uint32_t id, const std::string& s) const {
        Entry e;
        e.id = id;
        std::vector<uint32_t> v, sg;
        StringMapper::two_resolution(s, cfg_.n, cfg_.L, cfg_.sig, v, sg);
        e.key = z_.encode(v);
        if (cfg_.sig > 0) {
            e.has_sig = 1;
            for (int d = 0; d < cfg_.sig; ++d) {
                if (sg[d] > 255) { e.has_sig = 0; break; }
                e.sig[d] = (uint8_t)sg[d];
            }
        }
        return e;
    }

private:
    struct Ctx {
        const std::string* q;
        int tau;
        std::vector<uint32_t> vq, sq, lb, ub, pt;
    };

    static BedConfig fix(BedConfig c) {
        if (c.B == 0) c.B = ZOrder::default_bits(c.L);
        return c;
    }

    void search(const BPlusTree<Entry>::Node* nd, Ctx& c, std::vector<uint32_t>& out, QueryStats& S) const {
        const std::string& q = *c.q;
        if (nd->leaf) {
            ++S.leaf_visited;
            for (size_t i = 0; i < nd->entries.size(); ++i) {
                const Entry& e = nd->entries[i];
                ++S.entries_scanned;
                // 第一層：由 key 還原 L 維向量 (精確值，飽和維除外) 算下界
                z_.decode(e.key, c.pt.data());
                if (z_.lower_bound_ed(c.vq.data(), c.pt.data(), c.pt.data(), cfg_.n) > c.tau) {
                    ++S.pruned_by_key;
                    continue;
                }
                // 第二層 (Bed+ 專有)：高維 signature 下界
                if (cfg_.sig > 0 && e.has_sig &&
                    sig_lower_bound(c.sq.data(), e.sig, cfg_.sig, cfg_.n) > c.tau) {
                    ++S.pruned_by_sig;
                    continue;
                }
                const std::string& s = strs_[e.id];
                if (std::abs((int)s.size() - (int)q.size()) > c.tau) { ++S.len_skipped; continue; }
                ++S.dp_calls;
                if (verify_ed(q, s, c.tau) <= c.tau) out.push_back(e.id);
            }
            return;
        }
        ++S.internal_visited;
        for (size_t i = 0; i < nd->kids.size(); ++i) {
            const auto& k = nd->kids[i];
            z_.bounds(k.lo.key, k.hi.key, c.lb.data(), c.ub.data());
            if (z_.lower_bound_ed(c.vq.data(), c.lb.data(), c.ub.data(), cfg_.n) > c.tau) {
                ++S.children_pruned;
                continue;
            }
            search(k.node.get(), c, out, S);
        }
    }

    BedConfig cfg_;
    ZOrder z_;
    BPlusTree<Entry> tree_;
    size_t leaf_entry_bytes_;
    std::vector<std::string> strs_;
    std::vector<char> alive_;
};

}  // namespace bed

#endif
