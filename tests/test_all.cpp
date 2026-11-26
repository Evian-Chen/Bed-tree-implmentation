// test_all：單元測試與「對照暴力搜尋」的正確性測試
//   test_all            全部執行
//   test_all ed|gram|zorder|bptree|bounds|bedtree   只跑某一組
// 最重要的是 bedtree 組：Bed-tree / Bed+-tree 的結果必須和暴力搜尋完全相同 (不漏、不多)。
#include <iostream>
#include <random>
#include <set>

#include "BedTree.h"
#include "BPlusTree.h"
#include "Cli.h"
#include "EditDistance.h"
#include "StringMapper.h"
#include "ZOrder.h"

using namespace bed;

static int g_fail = 0, g_checks = 0;
#define CHECK(cond)                                                                    \
    do {                                                                               \
        ++g_checks;                                                                    \
        if (!(cond)) { ++g_fail; std::cout << "  FAIL " << __FILE__ << ":" << __LINE__ << "  " #cond "\n"; } \
    } while (0)

static std::string rand_str(std::mt19937_64& g, int maxlen, int alpha) {
    int len = (int)(g() % (uint64_t)(maxlen + 1));
    std::string s;
    for (int i = 0; i < len; ++i) s += (char)('a' + g() % (uint64_t)alpha);
    return s;
}

static void test_ed() {
    std::cout << "[ed] edit distance\n";
    CHECK(edit_distance("kitten", "sitting") == 3);
    CHECK(edit_distance("", "abc") == 3);
    CHECK(edit_distance("abc", "") == 3);
    CHECK(edit_distance("same", "same") == 0);
    CHECK(edit_distance("flaw", "lawn") == 2);
    CHECK(edit_distance("Jim Gray", "Jam Gray") == 1);
    std::mt19937_64 g(1);
    for (int it = 0; it < 20000; ++it) {
        std::string a = rand_str(g, 14, 3), b = rand_str(g, 14, 3);
        int d = edit_distance(a, b);
        CHECK(edit_distance(b, a) == d);  // 對稱
        for (int tau = 0; tau <= 6; ++tau) {
            int v = verify_ed(a, b, tau);
            if (v != (d <= tau ? d : tau + 1)) {
                CHECK(false);
                std::cout << "   a=" << a << " b=" << b << " tau=" << tau << " d=" << d << " v=" << v << "\n";
            } else CHECK(true);
        }
    }
}

static void test_gram() {
    std::cout << "[gram] n-gram & bucket vectors\n";
    std::vector<std::string> gs = StringMapper::extract_ngrams("Jim Gray", 2);
    CHECK(gs.size() == 9);  // |s|+n-1
    CHECK(gs.front() == "#J" && gs.back() == "y#");
    CHECK(StringMapper::extract_ngrams("ab", 3).size() == 4);
    CHECK(StringMapper::extract_ngrams("", 2).empty());
    std::mt19937_64 g(2);
    for (int it = 0; it < 500; ++it) {
        std::string s = rand_str(g, 30, 26);
        if (s.empty()) continue;
        int n = 2 + (int)(g() % 3);
        for (int L = 1; L <= 64; L *= 2) {
            std::vector<uint32_t> v = StringMapper::bucket_vector(s, n, L);
            uint32_t sum = 0;
            for (size_t i = 0; i < v.size(); ++i) sum += v[i];
            CHECK(sum == s.size() + n - 1);
        }
        // L 整除 S 時，L 維向量就是 S 維向量的摺疊 (Bed+ 低維/高維一致)
        std::vector<uint32_t> v4, s64, v64u, dummy;
        StringMapper::two_resolution(s, n, 4, 64, v4, s64);
        for (int b = 0; b < 4; ++b) {
            uint32_t acc = 0;
            for (int d = b; d < 64; d += 4) acc += s64[d];
            CHECK(acc == v4[b]);
        }
    }
}

static void test_zorder() {
    std::cout << "[zorder] encode/decode/bounds\n";
    ZOrder z(4, 2);
    uint32_t v[4] = {3, 2, 1, 3};
    CHECK(z.encode(v).to_bit_string(z.T) == "11011011");  // 論文 Figure 2
    uint32_t out[4];
    z.decode(z.encode(v), out);
    CHECK(out[0] == 3 && out[1] == 2 && out[2] == 1 && out[3] == 3);
    // 論文 Figure 3：區間 [11011011, 11011110] 的第 1 個 bucket 必為 3
    uint32_t lb[4], ub[4];
    z.bounds(ZKey::from_bit_string("11011011"), ZKey::from_bit_string("11011110"), lb, ub);
    CHECK(lb[0] == 3 && ub[0] == 3);
    CHECK(lb[1] == 2 && ub[1] == 3);   // {2,3}
    CHECK(lb[2] == 0 && ub[2] == 1);   // {0,1}
    CHECK(lb[3] == 2 && ub[3] == 3);   // {2,3}
    // 飽和：超過 maxv 的計數被截斷為 maxv
    uint32_t big[4] = {100, 0, 0, 0};
    z.decode(z.encode(big), out);
    CHECK(out[0] == 3);

    std::mt19937_64 g(3);
    for (int it = 0; it < 3000; ++it) {
        int L = 1 + (int)(g() % 16);
        int B = 1 + (int)(g() % (uint64_t)ZOrder::default_bits(L));
        ZOrder zz(L, B);
        std::vector<uint32_t> a(L), b(L), c(L), d(L);
        for (int i = 0; i < L; ++i) { a[i] = g() % (zz.maxv + 3); b[i] = g() % (zz.maxv + 3); c[i] = g() % (zz.maxv + 3); }
        ZKey ka = zz.encode(a), kb = zz.encode(b), kc = zz.encode(c);
        zz.decode(ka, d.data());
        for (int i = 0; i < L; ++i) CHECK(d[i] == std::min(a[i], zz.maxv));
        ZKey lo = ka < kb ? ka : kb, hi = ka < kb ? kb : ka;
        if (!(kc < lo) && !(hi < kc)) {  // kc 落在區間內 -> 其解碼值必須落在 [lb,ub]
            std::vector<uint32_t> lbv(L), ubv(L);
            zz.bounds(lo, hi, lbv.data(), ubv.data());
            zz.decode(kc, d.data());
            for (int i = 0; i < L; ++i) CHECK(lbv[i] <= d[i] && d[i] <= ubv[i]);
        }
    }
}

static void test_bptree() {
    std::cout << "[bptree] insert/erase/bulk vs std::set, invariants\n";
    struct E { ZKey key; uint32_t id; E() : id(0) {} };
    typedef BPlusTree<E> T;
    std::mt19937_64 g(4);
    for (int cap = 2; cap <= 9; cap += 7) {
        T t(cap, cap < 3 ? 3 : cap);
        std::set<uint64_t> ref;
        for (int i = 0; i < 6000; ++i) {
            uint64_t v = g() % 3000;
            E e; e.key.w[3] = v; e.id = (uint32_t)v;
            if (ref.count(v)) { CHECK(t.erase(KeyId(e.key, e.id))); ref.erase(v); }
            else { t.insert(e); ref.insert(v); }
            if (i % 250 == 0) { std::string err; CHECK(t.check_invariants(&err)); if (!err.empty()) std::cout << "   " << err << "\n"; }
        }
        std::vector<uint64_t> got;
        t.for_each([&](const E& e) { got.push_back(e.key.w[3]); });
        CHECK(got.size() == ref.size() && std::equal(got.begin(), got.end(), ref.begin()));
        CHECK(!t.erase(KeyId(ZKey(), 123456)));  // 刪除不存在的
    }
    std::vector<E> es;
    for (int i = 0; i < 1000; ++i) { E e; e.key.w[3] = (uint64_t)i; e.id = (uint32_t)i; es.push_back(e); }
    T t(8, 5);
    t.bulk_load(es, 0.9);
    std::string err;
    CHECK(t.check_invariants(&err) && t.size() == 1000);
    for (int i = 0; i < 1000; ++i) { E e; e.key.w[3] = 1000 + (uint64_t)i * 2; e.id = 1000 + (uint32_t)i; t.insert(e); }
    CHECK(t.check_invariants(&err) && t.size() == 2000);
}

// 核心性質：下界永遠 <= 真實編輯距離 (否則會漏答案)
static void test_bounds() {
    std::cout << "[bounds] lower bounds never exceed true edit distance\n";
    std::mt19937_64 g(5);
    for (int it = 0; it < 4000; ++it) {
        int n = 2 + (int)(g() % 2);
        int L = 1 << (g() % 5);                       // 1,2,4,8,16
        int S = 1 << (3 + g() % 4);                   // 8..64
        int alpha = 2 + (int)(g() % 5);
        int maxlen = (it % 10 == 0) ? 90 : 12;        // 偶爾用長字串讓計數飽和
        std::string q = rand_str(g, maxlen, alpha), s1 = rand_str(g, maxlen, alpha), s2 = rand_str(g, maxlen, alpha),
                    s = rand_str(g, maxlen, alpha);
        ZOrder z(L, ZOrder::default_bits(L) > 3 ? 3 : ZOrder::default_bits(L));  // 故意用小 B 讓飽和更常發生
        std::vector<uint32_t> vq = StringMapper::bucket_vector(q, n, L);
        // (a) 點下界 (leaf 用 key 還原的向量)
        std::vector<uint32_t> vs = StringMapper::bucket_vector(s, n, L), pt(L);
        z.decode(z.encode(vs), pt.data());
        CHECK(z.lower_bound_ed(vq.data(), pt.data(), pt.data(), n) <= edit_distance(q, s));
        // (b) 區間下界：s 的 key 在 [k1,k2] 內時，LB(q,[k1,k2]) <= ED(q,s)
        ZKey k1 = z.encode(StringMapper::bucket_vector(s1, n, L)), k2 = z.encode(StringMapper::bucket_vector(s2, n, L));
        if (k2 < k1) std::swap(k1, k2);
        ZKey ks = z.encode(vs);
        if (!(ks < k1) && !(k2 < ks)) {
            std::vector<uint32_t> lb(L), ub(L);
            z.bounds(k1, k2, lb.data(), ub.data());
            CHECK(z.lower_bound_ed(vq.data(), lb.data(), ub.data(), n) <= edit_distance(q, s));
        }
        // (c) signature 下界 (只在 has_sig 時使用)
        std::vector<uint32_t> sq, ss, d1;
        StringMapper::two_resolution(q, n, 0, S, d1, sq);
        StringMapper::two_resolution(s, n, 0, S, d1, ss);
        bool ok = true;
        uint8_t sv[kMaxSig] = {0};
        for (int d = 0; d < S; ++d) { if (ss[d] > 255) ok = false; else sv[d] = (uint8_t)ss[d]; }
        if (ok) CHECK(sig_lower_bound(sq.data(), sv, S, n) <= edit_distance(q, s));
    }
}

static void check_tree(BedTree& t, const std::vector<std::string>& queries, int maxtau, const char* tag) {
    for (size_t i = 0; i < queries.size(); ++i)
        for (int tau = 0; tau <= maxtau; ++tau) {
            std::vector<uint32_t> a = t.range_query(queries[i], tau), b = t.scan_query(queries[i], tau);
            if (a != b) {
                ++g_fail;
                std::cout << "  FAIL(" << tag << ") q=\"" << queries[i] << "\" tau=" << tau << " tree=" << a.size()
                          << " scan=" << b.size() << "\n";
            }
            ++g_checks;
        }
}

static void test_bedtree() {
    std::cout << "[bedtree] Bed-tree / Bed+-tree == brute force (including after inserts/deletes)\n";
    std::mt19937_64 g(6);
    struct Case { int n, L, sig, alpha, maxlen, page; };
    Case cases[] = {{2, 4, 0, 4, 10, 256},  {2, 4, 64, 4, 10, 256}, {3, 4, 64, 5, 12, 512}, {2, 2, 32, 26, 14, 256},
                    {2, 8, 16, 3, 9, 128},  {2, 64, 0, 4, 10, 1024}, {2, 4, 64, 4, 120, 4096}, {2, 1, 64, 6, 8, 128}};
    for (size_t ci = 0; ci < sizeof(cases) / sizeof(cases[0]); ++ci) {
        const Case& c = cases[ci];
        std::vector<std::string> data;
        std::set<std::string> seen;
        while (data.size() < 2500) { std::string s = rand_str(g, c.maxlen, c.alpha); if (!s.empty() && seen.insert(s).second) data.push_back(s); }
        std::vector<std::string> queries;
        for (int i = 0; i < 40; ++i) queries.push_back(rand_str(g, c.maxlen, c.alpha));
        for (int i = 0; i < 40; ++i) queries.push_back(data[g() % data.size()]);
        BedConfig cfg; cfg.n = c.n; cfg.L = c.L; cfg.sig = c.sig; cfg.page_bytes = (size_t)c.page; cfg.string_bytes = 8;
        BedTree t(cfg);
        std::vector<std::string> half(data.begin(), data.begin() + 1500);
        t.build(half);
        std::string err;
        CHECK(t.tree().check_invariants(&err));
        check_tree(t, queries, 3, "after build");
        for (size_t i = 1500; i < data.size(); ++i) t.insert(data[i]);
        CHECK(t.tree().check_invariants(&err));
        check_tree(t, queries, 3, "after insert");
        for (int i = 0; i < 800; ++i) t.erase((uint32_t)(g() % data.size()));
        CHECK(t.tree().check_invariants(&err));
        CHECK(t.tree().size() == [&] { size_t a = 0; for (size_t i = 0; i < t.total_ids(); ++i) a += t.alive((uint32_t)i); return a; }());
        check_tree(t, queries, 3, "after erase");
        std::cout << "  case " << ci << " (n=" << c.n << " L=" << c.L << " sig=" << c.sig << ") done\n";
    }
}

int main(int argc, char** argv) {
    std::string which = argc > 1 ? argv[1] : "all";
    bool all = which == "all";
    if (all || which == "ed") test_ed();
    if (all || which == "gram") test_gram();
    if (all || which == "zorder") test_zorder();
    if (all || which == "bptree") test_bptree();
    if (all || which == "bounds") test_bounds();
    if (all || which == "bedtree") test_bedtree();
    std::cout << "\n" << (g_checks - g_fail) << "/" << g_checks << " checks passed, " << g_fail << " failed\n";
    return g_fail ? 1 : 0;
}
