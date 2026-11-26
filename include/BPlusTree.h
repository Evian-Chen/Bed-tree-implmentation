// 自製的記憶體內 B+-tree (不使用任何資料庫系統)。
// 與一般 B+-tree 的差異：每個內部節點的 child 都記錄子樹實際的 [lo, hi] key，
// 並開放節點讓 Bed-tree 自訂「以 lower bound 剪枝」的走訪。
// 刪除採 lazy 策略：只移除空節點，不做 merge/redistribute (期末專案的簡化，見 README)。
#ifndef BED_BPLUSTREE_H
#define BED_BPLUSTREE_H

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "ZOrder.h"

namespace bed {

// 排序鍵：(z-order key, 字串 id)；id 用來讓相同 key 的字串有全序
struct KeyId {
    ZKey key;
    uint32_t id;
    KeyId() : id(0) {}
    KeyId(const ZKey& k, uint32_t i) : key(k), id(i) {}
    bool operator<(const KeyId& o) const { return key < o.key || (key == o.key && id < o.id); }
    bool operator==(const KeyId& o) const { return key == o.key && id == o.id; }
};

// Entry 需有成員 `ZKey key; uint32_t id;`
template <class Entry>
class BPlusTree {
public:
    struct Node;
    struct Child {
        KeyId lo, hi;               // 此子樹內最小/最大 (key,id)
        std::unique_ptr<Node> node;
    };
    struct Node {
        bool leaf;
        std::vector<Entry> entries;  // leaf 用
        std::vector<Child> kids;     // internal 用
        explicit Node(bool l) : leaf(l) {}
    };

    BPlusTree(size_t leaf_cap = 64, size_t internal_cap = 64)
        : leaf_cap_(leaf_cap < 2 ? 2 : leaf_cap),
          internal_cap_(internal_cap < 3 ? 3 : internal_cap),
          size_(0) {
        root_.reset(new Node(true));
    }

    size_t size() const { return size_; }
    const Node* root() const { return root_.get(); }
    size_t leaf_capacity() const { return leaf_cap_; }
    size_t internal_capacity() const { return internal_cap_; }

    static KeyId kid(const Entry& e) { return KeyId(e.key, e.id); }

    // 由已排序 (依 KeyId) 的 entries 一次建好整棵樹；fill 為每個節點的填充率
    void bulk_load(const std::vector<Entry>& sorted, double fill = 0.9) {
        root_.reset(new Node(true));
        size_ = 0;
        if (sorted.empty()) return;
        size_t lc = std::max<size_t>(1, (size_t)(leaf_cap_ * fill));
        size_t ic = std::max<size_t>(2, (size_t)(internal_cap_ * fill));
        std::vector<Child> level;
        for (size_t i = 0; i < sorted.size(); i += lc) {
            std::unique_ptr<Node> nd(new Node(true));
            size_t e = std::min(sorted.size(), i + lc);
            nd->entries.assign(sorted.begin() + i, sorted.begin() + e);
            Child c;
            c.node = std::move(nd);
            refresh(c);
            level.push_back(std::move(c));
        }
        while (level.size() > 1) {
            std::vector<Child> next;
            for (size_t i = 0; i < level.size(); i += ic) {
                std::unique_ptr<Node> nd(new Node(false));
                size_t e = std::min(level.size(), i + ic);
                for (size_t j = i; j < e; ++j) nd->kids.push_back(std::move(level[j]));
                Child c;
                c.node = std::move(nd);
                refresh(c);
                next.push_back(std::move(c));
            }
            level.swap(next);
        }
        root_ = std::move(level[0].node);
        size_ = sorted.size();
    }

    void insert(const Entry& e) {
        std::unique_ptr<Node> sib = insert_rec(root_.get(), e);
        if (sib) {  // 根分裂，樹長高
            std::unique_ptr<Node> nr(new Node(false));
            Child a, b;
            a.node = std::move(root_);
            b.node = std::move(sib);
            refresh(a);
            refresh(b);
            nr->kids.push_back(std::move(a));
            nr->kids.push_back(std::move(b));
            root_ = std::move(nr);
        }
        ++size_;
    }

    // 刪除 (key,id)；找不到回傳 false
    bool erase(const KeyId& k) {
        bool found = erase_rec(root_.get(), k);
        if (!found) return false;
        --size_;
        while (!root_->leaf && root_->kids.size() == 1) {
            std::unique_ptr<Node> only = std::move(root_->kids[0].node);
            root_ = std::move(only);
        }
        if (!root_->leaf && root_->kids.empty()) root_.reset(new Node(true));
        return true;
    }

    int height() const {
        int h = 1;
        const Node* n = root_.get();
        while (!n->leaf) { n = n->kids[0].node.get(); ++h; }
        return h;
    }

    void count_nodes(size_t* leaves, size_t* internals) const {
        *leaves = 0; *internals = 0;
        count_rec(root_.get(), leaves, internals);
    }

    // 依序走訪所有 entry
    template <class F>
    void for_each(F f) const { each_rec(root_.get(), f); }

    // 檢查結構不變量；失敗時 err 裡有原因
    bool check_invariants(std::string* err) const {
        std::ostringstream os;
        int leaf_depth = -1;
        size_t cnt = 0;
        bool ok = check_rec(root_.get(), 1, &leaf_depth, &cnt, true, os);
        if (ok && cnt != size_) { os << "size mismatch: counted " << cnt << " vs size_ " << size_; ok = false; }
        if (!ok && err) *err = os.str();
        return ok;
    }

private:
    size_t leaf_cap_, internal_cap_, size_;
    std::unique_ptr<Node> root_;

    static KeyId min_of(const Node* n) { return n->leaf ? kid(n->entries.front()) : n->kids.front().lo; }
    static KeyId max_of(const Node* n) { return n->leaf ? kid(n->entries.back()) : n->kids.back().hi; }
    static void refresh(Child& c) { c.lo = min_of(c.node.get()); c.hi = max_of(c.node.get()); }

    // 找到該往哪個 child：最後一個 lo <= k 的 child (都比 k 大則取 0)
    static size_t child_index(const Node* nd, const KeyId& k) {
        size_t lo = 0, hi = nd->kids.size();  // 找第一個 lo > k
        while (lo < hi) {
            size_t mid = (lo + hi) / 2;
            if (k < nd->kids[mid].lo) hi = mid; else lo = mid + 1;
        }
        return lo == 0 ? 0 : lo - 1;
    }

    std::unique_ptr<Node> insert_rec(Node* nd, const Entry& e) {
        if (nd->leaf) {
            KeyId k = kid(e);
            size_t lo = 0, hi = nd->entries.size();
            while (lo < hi) {
                size_t mid = (lo + hi) / 2;
                if (k < kid(nd->entries[mid])) hi = mid; else lo = mid + 1;
            }
            nd->entries.insert(nd->entries.begin() + lo, e);
            if (nd->entries.size() <= leaf_cap_) return nullptr;
            size_t mid = nd->entries.size() / 2;
            std::unique_ptr<Node> sib(new Node(true));
            sib->entries.assign(nd->entries.begin() + mid, nd->entries.end());
            nd->entries.erase(nd->entries.begin() + mid, nd->entries.end());
            return sib;
        }
        size_t idx = child_index(nd, kid(e));
        std::unique_ptr<Node> sib = insert_rec(nd->kids[idx].node.get(), e);
        refresh(nd->kids[idx]);
        if (sib) {
            Child c;
            c.node = std::move(sib);
            refresh(c);
            nd->kids.insert(nd->kids.begin() + idx + 1, std::move(c));
        }
        if (nd->kids.size() <= internal_cap_) return nullptr;
        size_t mid = nd->kids.size() / 2;
        std::unique_ptr<Node> right(new Node(false));
        for (size_t j = mid; j < nd->kids.size(); ++j) right->kids.push_back(std::move(nd->kids[j]));
        nd->kids.erase(nd->kids.begin() + mid, nd->kids.end());
        return right;
    }

    bool erase_rec(Node* nd, const KeyId& k) {
        if (nd->leaf) {
            for (size_t i = 0; i < nd->entries.size(); ++i)
                if (kid(nd->entries[i]) == k) {
                    nd->entries.erase(nd->entries.begin() + i);
                    return true;
                }
            return false;
        }
        size_t idx = child_index(nd, k);
        Node* ch = nd->kids[idx].node.get();
        if (!erase_rec(ch, k)) return false;
        bool empty = ch->leaf ? ch->entries.empty() : ch->kids.empty();
        if (empty) nd->kids.erase(nd->kids.begin() + idx);
        else refresh(nd->kids[idx]);
        return true;
    }

    void count_rec(const Node* n, size_t* l, size_t* in) const {
        if (n->leaf) { ++*l; return; }
        ++*in;
        for (size_t i = 0; i < n->kids.size(); ++i) count_rec(n->kids[i].node.get(), l, in);
    }

    template <class F>
    void each_rec(const Node* n, F& f) const {
        if (n->leaf) { for (size_t i = 0; i < n->entries.size(); ++i) f(n->entries[i]); return; }
        for (size_t i = 0; i < n->kids.size(); ++i) each_rec(n->kids[i].node.get(), f);
    }

    bool check_rec(const Node* n, int depth, int* leaf_depth, size_t* cnt, bool is_root,
                   std::ostringstream& os) const {
        if (n->leaf) {
            if (*leaf_depth == -1) *leaf_depth = depth;
            else if (*leaf_depth != depth) { os << "leaves at different depths"; return false; }
            if (n->entries.empty() && !is_root) { os << "empty leaf"; return false; }
            if (n->entries.size() > leaf_cap_) { os << "leaf over capacity"; return false; }
            for (size_t i = 1; i < n->entries.size(); ++i)
                if (!(kid(n->entries[i - 1]) < kid(n->entries[i]))) { os << "leaf not strictly sorted"; return false; }
            *cnt += n->entries.size();
            return true;
        }
        if (n->kids.empty()) { os << "empty internal node"; return false; }
        if (n->kids.size() > internal_cap_) { os << "internal over capacity"; return false; }
        for (size_t i = 0; i < n->kids.size(); ++i) {
            const Child& c = n->kids[i];
            if (!check_rec(c.node.get(), depth + 1, leaf_depth, cnt, false, os)) return false;
            if (!(c.lo == min_of(c.node.get())) || !(c.hi == max_of(c.node.get()))) {
                os << "stale child lo/hi"; return false;
            }
            if (i > 0 && !(n->kids[i - 1].hi < c.lo)) { os << "children ranges overlap/unsorted"; return false; }
        }
        return true;
    }
};

}  // namespace bed

#endif
