// bench：實驗程式。每個實驗都會「對照暴力搜尋」檢查結果集合完全相同，並輸出表格與 CSV (results/)。
//   bench e1 --data F     維度兩難：原版 Bed-tree 在不同 L 下的候選數 vs 節點數
//   bench e2 --data F     主要比較：scan vs Bed-tree vs Bed+-tree (不同 tau)
//   bench e3 --data F     參數格：排序維度 L × signature 維度 S
//   bench e4 --data F     更新：逐筆插入 / 刪除吞吐量，以及更新後查詢仍正確
//   bench e5 --data F     gram 長度 n 的影響
//   bench e6 --data F     資料量 (scalability)
//   bench all --data F    全部
// 共同選項： --qcount 200 (每個 tau 的查詢數)  --taus 1,2,3  --seed 7  --limit N (只用前 N 筆資料)
//            --page 4096  --reps 3 (取最短時間)  --out results (CSV 目錄)
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <random>

#include "BedTree.h"
#include "Cli.h"

using namespace bed;

struct Row {
    std::string exp, label;
    int tau = 0;
    double us = 0, build_s = 0;
    double dp = 0, nodes = 0, internal = 0, leaf = 0, scanned = 0, p_key = 0, p_sig = 0, results = 0;
    size_t pages = 0, height = 0, leaf_cap = 0;
    int mismatch = 0;
};

static std::string g_out = "results";
static int g_total_mismatch = 0;

static std::string mutate(const std::string& s, int edits, std::mt19937_64& g) {
    std::string t = s;
    for (int e = 0; e < edits; ++e) {
        int op = (int)(g() % 3);
        char c = (char)('a' + g() % 26);
        if (t.empty()) op = 0;
        if (op == 0) t.insert(t.begin() + (long)(g() % (t.size() + 1)), c);
        else if (op == 1) t.erase(t.begin() + (long)(g() % t.size()));
        else t[g() % t.size()] = c;
    }
    return t;
}

static std::vector<std::string> make_queries(const std::vector<std::string>& data, int count, int edits, uint64_t seed) {
    std::mt19937_64 g(seed + (uint64_t)edits * 1000003ULL);
    std::vector<std::string> q;
    for (int i = 0; i < count; ++i) q.push_back(mutate(data[g() % data.size()], edits, g));
    return q;
}

typedef std::vector<std::vector<uint32_t> > Results;

static Results reference(const BedTree& t, const std::vector<std::string>& qs, int tau) {
    Results r;
    for (size_t i = 0; i < qs.size(); ++i) r.push_back(t.scan_query(qs[i], tau));
    return r;
}

static Row eval(const std::string& exp, const std::string& label, const BedTree& t, bool scan,
                const std::vector<std::string>& qs, int tau, const Results& ref, double build_s, int reps) {
    Row r;
    r.exp = exp; r.label = label; r.tau = tau; r.build_s = build_s;
    // 暖機
    for (size_t i = 0; i < qs.size() && i < 5; ++i) scan ? t.scan_query(qs[i], tau) : t.range_query(qs[i], tau);
    double best = 1e30;
    QueryStats st;
    Results got;
    for (int rep = 0; rep < reps; ++rep) {
        QueryStats cur;
        Results res;
        res.reserve(qs.size());
        Timer tm;
        for (size_t i = 0; i < qs.size(); ++i) res.push_back(scan ? t.scan_query(qs[i], tau, &cur) : t.range_query(qs[i], tau, &cur));
        double s = tm.seconds();
        if (s < best) best = s;
        if (rep == 0) { st = cur; got = res; }
    }
    double nq = (double)qs.size();
    r.us = best / nq * 1e6;
    r.dp = st.dp_calls / nq; r.nodes = st.nodes() / nq; r.internal = st.internal_visited / nq;
    r.leaf = st.leaf_visited / nq; r.scanned = st.entries_scanned / nq;
    r.p_key = st.pruned_by_key / nq; r.p_sig = st.pruned_by_sig / nq; r.results = st.results / nq;
    for (size_t i = 0; i < qs.size(); ++i) if (got[i] != ref[i]) ++r.mismatch;
    g_total_mismatch += r.mismatch;
    if (!scan) {
        StorageInfo si = t.storage();
        r.pages = si.total_pages(); r.height = si.height; r.leaf_cap = si.leaf_capacity;
    }
    return r;
}

static void print_header() {
    std::cout << std::left << std::setw(5) << "exp" << std::setw(18) << "config" << std::right << std::setw(4) << "tau"
              << std::setw(10) << "us/query" << std::setw(10) << "DP/query" << std::setw(10) << "nodes/q"
              << std::setw(9) << "internal" << std::setw(8) << "leaf" << std::setw(10) << "scanned"
              << std::setw(10) << "pr_key" << std::setw(10) << "pr_sig" << std::setw(9) << "results"
              << std::setw(8) << "pages" << std::setw(7) << "h" << std::setw(8) << "build_s" << std::setw(5) << "bad"
              << "\n";
}

static void print_row(const Row& r) {
    std::cout << std::left << std::setw(5) << r.exp << std::setw(18) << r.label << std::right << std::setw(4) << r.tau
              << std::fixed << std::setprecision(1) << std::setw(10) << r.us << std::setw(10) << r.dp
              << std::setw(10) << r.nodes << std::setw(9) << r.internal << std::setw(8) << r.leaf
              << std::setw(10) << r.scanned << std::setw(10) << r.p_key << std::setw(10) << r.p_sig
              << std::setw(9) << r.results << std::setw(8) << r.pages << std::setw(7) << r.height
              << std::setprecision(2) << std::setw(8) << r.build_s << std::setw(5) << r.mismatch << "\n";
}

static void write_csv(const std::string& name, const std::vector<Row>& rows) {
    std::string path = g_out + "/" + name + ".csv";
    std::ofstream f(path.c_str());
    if (!f) { std::cerr << "cannot write " << path << " (create the directory first)\n"; return; }
    f << "exp,config,tau,us_per_query,dp_per_query,nodes_per_query,internal_per_query,leaf_per_query,entries_scanned,pruned_by_key,pruned_by_sig,results_per_query,pages,height,leaf_cap,build_s,mismatch_queries\n";
    for (size_t i = 0; i < rows.size(); ++i) {
        const Row& r = rows[i];
        f << r.exp << "," << r.label << "," << r.tau << "," << r.us << "," << r.dp << "," << r.nodes << ","
          << r.internal << "," << r.leaf << "," << r.scanned << "," << r.p_key << "," << r.p_sig << "," << r.results
          << "," << r.pages << "," << r.height << "," << r.leaf_cap << "," << r.build_s << "," << r.mismatch << "\n";
    }
    std::cout << "  -> " << path << "\n";
}

struct Ctx {
    std::vector<std::string> data;
    std::vector<int> taus;
    int qcount, reps;
    uint64_t seed;
    size_t page;
    int n;
};

static BedTree* build_tree(const Ctx& c, int L, int S, int n, double* secs, const std::vector<std::string>& data) {
    BedConfig cfg;
    cfg.n = n; cfg.L = L; cfg.sig = S; cfg.page_bytes = c.page;
    double avg = 0;
    for (size_t i = 0; i < data.size(); ++i) avg += (double)data[i].size();
    cfg.string_bytes = (size_t)(avg / (double)data.size() + 0.5);  // leaf 內字串佔用，依資料平均長度
    BedTree* t = new BedTree(cfg);
    Timer tm;
    t->build(data);
    *secs = tm.seconds();
    return t;
}

static std::string cfg_label(int L, int S) {
    return "L=" + itos(L) + (S ? ",S=" + itos(S) : ",S=-");
}

// 每個 tau 用「該 tau 次編輯」的查詢，並用暴力搜尋建立參考答案
struct TauSet { int tau; std::vector<std::string> qs; Results ref; };

static std::vector<TauSet> make_sets(const Ctx& c, const BedTree& ref_tree, const std::vector<std::string>& data) {
    std::vector<TauSet> v;
    for (size_t i = 0; i < c.taus.size(); ++i) {
        TauSet s;
        s.tau = c.taus[i];
        s.qs = make_queries(data, c.qcount, s.tau, c.seed);
        s.ref = reference(ref_tree, s.qs, s.tau);
        v.push_back(s);
    }
    return v;
}

static void e1(const Ctx& c) {
    std::cout << "\n=== E1: dimensionality dilemma (original Bed-tree, sweep L) ===\n";
    double sec;
    std::unique_ptr<BedTree> ref(build_tree(c, 4, 0, c.n, &sec, c.data));
    std::vector<TauSet> sets = make_sets(c, *ref, c.data);
    std::vector<Row> rows;
    print_header();
    int Ls[] = {1, 2, 4, 8, 16, 32, 64};
    for (size_t li = 0; li < sizeof(Ls) / sizeof(Ls[0]); ++li)
        for (size_t k = 0; k < sets.size(); ++k) {
            std::unique_ptr<BedTree> t(build_tree(c, Ls[li], 0, c.n, &sec, c.data));
            Row r = eval("E1", "Bed L=" + itos(Ls[li]), *t, false, sets[k].qs, sets[k].tau, sets[k].ref, sec, c.reps);
            rows.push_back(r); print_row(r);
        }
    write_csv("e1_dilemma", rows);
}

static void e2(const Ctx& c) {
    std::cout << "\n=== E2: scan vs Bed-tree vs Bed+-tree ===\n";
    double sec;
    std::unique_ptr<BedTree> scan_t(build_tree(c, 4, 0, c.n, &sec, c.data));
    std::vector<TauSet> sets = make_sets(c, *scan_t, c.data);
    std::vector<Row> rows;
    print_header();
    for (size_t k = 0; k < sets.size(); ++k) {
        Row r = eval("E2", "scan", *scan_t, true, sets[k].qs, sets[k].tau, sets[k].ref, 0, c.reps);
        rows.push_back(r); print_row(r);
        std::unique_ptr<BedTree> bed(build_tree(c, 4, 0, c.n, &sec, c.data));
        r = eval("E2", "Bed (L=4)", *bed, false, sets[k].qs, sets[k].tau, sets[k].ref, sec, c.reps);
        rows.push_back(r); print_row(r);
        std::unique_ptr<BedTree> bp(build_tree(c, 4, 64, c.n, &sec, c.data));
        r = eval("E2", "Bed+ (L=4,S=64)", *bp, false, sets[k].qs, sets[k].tau, sets[k].ref, sec, c.reps);
        rows.push_back(r); print_row(r);
    }
    write_csv("e2_main", rows);
}

static void e3(const Ctx& c) {
    std::cout << "\n=== E3: ordering dimension L x signature dimension S ===\n";
    double sec;
    std::unique_ptr<BedTree> ref(build_tree(c, 4, 0, c.n, &sec, c.data));
    int tau = c.taus.size() > 1 ? c.taus[1] : c.taus[0];
    std::vector<std::string> qs = make_queries(c.data, c.qcount, tau, c.seed);
    Results rf = reference(*ref, qs, tau);
    std::vector<Row> rows;
    print_header();
    int Ls[] = {2, 4, 8, 16};
    int Ss[] = {0, 16, 32, 64};
    for (size_t li = 0; li < 4; ++li)
        for (size_t si = 0; si < 4; ++si) {
            std::unique_ptr<BedTree> t(build_tree(c, Ls[li], Ss[si], c.n, &sec, c.data));
            Row r = eval("E3", cfg_label(Ls[li], Ss[si]), *t, false, qs, tau, rf, sec, c.reps);
            rows.push_back(r); print_row(r);
        }
    write_csv("e3_grid", rows);
}

static void e4(const Ctx& c) {
    std::cout << "\n=== E4: updates (insert / delete throughput, correctness afterwards) ===\n";
    std::vector<Row> rows;
    std::vector<std::string> shuffled = c.data;
    std::mt19937_64 g(c.seed);
    std::shuffle(shuffled.begin(), shuffled.end(), g);
    size_t half = shuffled.size() / 2;
    std::vector<std::string> first(shuffled.begin(), shuffled.begin() + (long)half);
    int S_list[] = {0, 64};
    std::cout << std::left << std::setw(18) << "config" << std::right << std::setw(14) << "bulk_load_s"
              << std::setw(16) << "insert_us/op" << std::setw(16) << "insert_ops/s" << std::setw(16) << "delete_us/op"
              << std::setw(16) << "delete_ops/s" << std::setw(10) << "pages" << std::setw(12) << "query_bad" << "\n";
    for (size_t si = 0; si < 2; ++si) {
        BedConfig cfg;
        cfg.n = c.n; cfg.L = 4; cfg.sig = S_list[si]; cfg.page_bytes = c.page;
        double avg = 0;
        for (size_t i = 0; i < c.data.size(); ++i) avg += (double)c.data[i].size();
        cfg.string_bytes = (size_t)(avg / (double)c.data.size() + 0.5);
        BedTree t(cfg);
        Timer tb;
        t.build(first);
        double bulk = tb.seconds();
        Timer ti;
        for (size_t i = half; i < shuffled.size(); ++i) t.insert(shuffled[i]);
        double ins = ti.seconds();
        size_t ninsert = shuffled.size() - half;
        // 隨機刪除 20% 的字串
        std::vector<uint32_t> ids;
        for (uint32_t i = 0; i < (uint32_t)t.total_ids(); ++i) ids.push_back(i);
        std::shuffle(ids.begin(), ids.end(), g);
        size_t ndel = ids.size() / 5;
        Timer td;
        for (size_t i = 0; i < ndel; ++i) t.erase(ids[i]);
        double del = td.seconds();
        std::string err;
        bool inv = t.tree().check_invariants(&err);
        // 更新後的查詢仍須等於暴力搜尋
        int bad = 0;
        for (size_t k = 0; k < c.taus.size(); ++k) {
            std::vector<std::string> qs = make_queries(c.data, c.qcount / 2 > 0 ? c.qcount / 2 : 1, c.taus[k], c.seed + 99);
            for (size_t i = 0; i < qs.size(); ++i)
                if (t.range_query(qs[i], c.taus[k]) != t.scan_query(qs[i], c.taus[k])) ++bad;
        }
        if (!inv) { std::cout << "INVARIANT FAILED: " << err << "\n"; ++bad; }
        g_total_mismatch += bad;
        std::string label = S_list[si] ? "Bed+ (L=4,S=64)" : "Bed (L=4)";
        std::cout << std::left << std::setw(18) << label << std::right << std::fixed << std::setprecision(3)
                  << std::setw(14) << bulk << std::setprecision(2) << std::setw(16) << ins / (double)ninsert * 1e6
                  << std::setprecision(0) << std::setw(16) << (double)ninsert / ins << std::setprecision(2)
                  << std::setw(16) << del / (double)ndel * 1e6 << std::setprecision(0) << std::setw(16)
                  << (double)ndel / del << std::setw(10) << t.storage().total_pages() << std::setw(12) << bad << "\n";
        Row r; r.exp = "E4"; r.label = label; r.us = ins / (double)ninsert * 1e6; r.build_s = bulk;
        r.dp = del / (double)ndel * 1e6;  // 借用欄位：delete us/op
        r.pages = t.storage().total_pages(); r.mismatch = bad;
        rows.push_back(r);
    }
    std::string path = g_out + "/e4_updates.csv";
    std::ofstream f(path.c_str());
    f << "config,bulk_load_s,insert_us_per_op,delete_us_per_op,pages,query_mismatches\n";
    for (size_t i = 0; i < rows.size(); ++i)
        f << rows[i].label << "," << rows[i].build_s << "," << rows[i].us << "," << rows[i].dp << "," << rows[i].pages << "," << rows[i].mismatch << "\n";
    std::cout << "  -> " << path << "\n";
}

static void e5(const Ctx& c) {
    std::cout << "\n=== E5: gram length n ===\n";
    double sec;
    int tau = c.taus.size() > 1 ? c.taus[1] : c.taus[0];
    std::unique_ptr<BedTree> ref(build_tree(c, 4, 0, 2, &sec, c.data));
    std::vector<std::string> qs = make_queries(c.data, c.qcount, tau, c.seed);
    Results rf = reference(*ref, qs, tau);
    std::vector<Row> rows;
    print_header();
    for (int n = 2; n <= 4; ++n)
        for (int S = 0; S <= 64; S += 64) {
            std::unique_ptr<BedTree> t(build_tree(c, 4, S, n, &sec, c.data));
            Row r = eval("E5", "n=" + itos(n) + (S ? ",Bed+" : ",Bed"), *t, false, qs, tau, rf, sec, c.reps);
            rows.push_back(r); print_row(r);
        }
    write_csv("e5_gram_n", rows);
}

static void e6(const Ctx& c, const std::vector<int>& pct) {
    std::cout << "\n=== E6: data size scalability ===\n";
    std::vector<Row> rows;
    print_header();
    int tau = c.taus.size() > 1 ? c.taus[1] : c.taus[0];
    for (size_t p = 0; p < pct.size(); ++p) {
        std::vector<std::string> sub(c.data.begin(), c.data.begin() + (long)(c.data.size() * (size_t)pct[p] / 100));
        double sec;
        std::unique_ptr<BedTree> st(build_tree(c, 4, 0, c.n, &sec, sub));
        std::vector<std::string> qs = make_queries(sub, c.qcount, tau, c.seed);
        Results rf = reference(*st, qs, tau);
        std::string tag = itos(pct[p]) + "%";
        Row r = eval("E6", "scan " + tag, *st, true, qs, tau, rf, 0, c.reps);
        rows.push_back(r); print_row(r);
        std::unique_ptr<BedTree> b(build_tree(c, 4, 0, c.n, &sec, sub));
        r = eval("E6", "Bed " + tag, *b, false, qs, tau, rf, sec, c.reps);
        rows.push_back(r); print_row(r);
        std::unique_ptr<BedTree> bp(build_tree(c, 4, 64, c.n, &sec, sub));
        r = eval("E6", "Bed+ " + tag, *bp, false, qs, tau, rf, sec, c.reps);
        rows.push_back(r); print_row(r);
    }
    write_csv("e6_scale", rows);
}

int main(int argc, char** argv) {
    Args a(argc, argv);
    if (a.pos.empty() || !a.has("data")) {
        std::cerr << "usage: bench e1|e2|e3|e4|e5|e6|all --data FILE [--qcount 200] [--taus 1,2,3] [--limit N]\n";
        return 2;
    }
    Ctx c;
    c.data = load_lines(a.get("data"));
    if (a.has("limit") && (size_t)a.geti("limit", 0) < c.data.size()) c.data.resize((size_t)a.geti("limit", 0));
    c.taus = parse_int_list(a.get("taus", "1,2,3"));
    c.qcount = a.geti("qcount", 200);
    c.reps = a.geti("reps", 3);
    c.seed = (uint64_t)a.geti("seed", 7);
    c.page = (size_t)a.geti("page", 4096);
    c.n = a.geti("n", 2);
    g_out = a.get("out", "results");
    double avg = 0;
    size_t mx = 0;
    for (size_t i = 0; i < c.data.size(); ++i) { avg += (double)c.data[i].size(); mx = std::max(mx, c.data[i].size()); }
    std::cout << "dataset: " << a.get("data") << "  strings=" << c.data.size() << "  avg_len=" << std::fixed
              << std::setprecision(1) << avg / (double)c.data.size() << "  max_len=" << mx << "\n";
    std::cout << "queries per tau=" << c.qcount << "  page=" << c.page << "B  n=" << c.n << "  reps=" << c.reps << "\n";
    std::string w = a.pos[0];
    bool all = w == "all";
    if (all || w == "e1") e1(c);
    if (all || w == "e2") e2(c);
    if (all || w == "e3") e3(c);
    if (all || w == "e4") e4(c);
    if (all || w == "e5") e5(c);
    if (all || w == "e6") { std::vector<int> p; p.push_back(25); p.push_back(50); p.push_back(100); e6(c, p); }
    std::cout << "\nTOTAL queries whose result differed from brute force: " << g_total_mismatch << "\n";
    return g_total_mismatch ? 1 : 0;
}
