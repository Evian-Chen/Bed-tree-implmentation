// bedtree_cli：建索引並執行範圍查詢 (ED <= tau)
//   bedtree_cli --data data/sample_words.txt --q "recieve" --tau 2
//   bedtree_cli --data F --q "recieve" --tau 2 --mode bedplus --L 4 --sig 64
//   bedtree_cli --data F --queries data/sample_queries.txt --tau 1 --compare
// 選項：
//   --mode scan|bed|bedplus   (預設 --compare 同時跑三種；單獨指定時只跑一種)
//   --n 2   gram 長度         --L 4   排序用 bucket 數 (低維)
//   --sig 64  Bed+ 的 signature 維度 (--mode bedplus 預設 64)
//   --page 4096  模擬 page 大小      --show 10  最多列出幾筆結果
#include <iomanip>
#include <iostream>

#include "BedTree.h"
#include "Cli.h"

using namespace bed;

struct Run {
    std::string name;
    BedConfig cfg;
    bool scan;
};

static void print_stats(const std::string& name, const QueryStats& s, size_t nq, double sec) {
    std::cout << "  " << std::left << std::setw(8) << name << std::right
              << " time/query=" << std::fixed << std::setprecision(1) << (sec / nq * 1e6) << "us"
              << "  nodes=" << (double)s.nodes() / nq << " (internal " << (double)s.internal_visited / nq
              << ", leaf " << (double)s.leaf_visited / nq << ")"
              << "  entries_scanned=" << (double)s.entries_scanned / nq
              << "  pruned_key=" << (double)s.pruned_by_key / nq
              << "  pruned_sig=" << (double)s.pruned_by_sig / nq
              << "  len_skip=" << (double)s.len_skipped / nq
              << "  DP=" << (double)s.dp_calls / nq
              << "  results=" << (double)s.results / nq << "\n";
}

int main(int argc, char** argv) {
    Args a(argc, argv);
    if (!a.has("data") || (!a.has("q") && !a.has("queries"))) {
        std::cerr << "usage: bedtree_cli --data F (--q STR | --queries F) --tau T [--mode scan|bed|bedplus] [--compare]\n";
        return 2;
    }
    std::vector<std::string> data = load_lines(a.get("data"));
    std::vector<std::string> queries;
    if (a.has("q")) queries.push_back(a.get("q"));
    else queries = load_lines(a.get("queries"));
    int tau = a.geti("tau", 1);
    int show = a.geti("show", 10);

    BedConfig base;
    base.n = a.geti("n", 2);
    base.L = a.geti("L", 4);
    base.page_bytes = (size_t)a.geti("page", 4096);

    std::vector<Run> runs;
    std::string mode = a.get("mode", a.has("compare") ? "all" : "bedplus");
    BedConfig c0 = base; c0.sig = 0;
    BedConfig c1 = base; c1.sig = a.geti("sig", 64);
    if (mode == "scan" || mode == "all") runs.push_back({"scan", c0, true});
    if (mode == "bed" || mode == "all") runs.push_back({"bed", c0, false});
    if (mode == "bedplus" || mode == "all") runs.push_back({"bed+", c1, false});
    if (runs.empty()) { std::cerr << "bad --mode\n"; return 2; }

    std::cout << "data=" << data.size() << " strings, queries=" << queries.size() << ", tau=" << tau
              << ", n=" << base.n << ", L=" << base.L << ", page=" << base.page_bytes << "B\n";

    std::vector<std::vector<uint32_t> > reference;  // 第一個 run 的結果，用來互相比對
    bool mismatch = false;
    for (size_t r = 0; r < runs.size(); ++r) {
        BedTree t(runs[r].cfg);
        Timer tb;
        t.build(data);
        double build_s = tb.seconds();
        QueryStats total;
        Timer tq;
        std::vector<std::vector<uint32_t> > res;
        for (size_t i = 0; i < queries.size(); ++i)
            res.push_back(runs[r].scan ? t.scan_query(queries[i], tau, &total)
                                       : t.range_query(queries[i], tau, &total));
        double q_s = tq.seconds();
        StorageInfo si = t.storage();
        std::cout << "\n[" << runs[r].name << "] build=" << std::fixed << std::setprecision(3) << build_s << "s";
        if (!runs[r].scan)
            std::cout << "  pages=" << si.total_pages() << " (leaf " << si.leaves << ", internal " << si.internals
                      << ", height " << si.height << ", leaf_cap " << si.leaf_capacity << ")";
        std::cout << "\n";
        print_stats(runs[r].name, total, queries.size(), q_s);
        if (r == 0) reference = res;
        else if (res != reference) { mismatch = true; std::cout << "  !! results differ from first run\n"; }
        if (queries.size() == 1 && r + 1 == runs.size()) {
            std::cout << "\nmatches for \"" << queries[0] << "\" (" << res[0].size() << " total):\n";
            for (size_t i = 0; i < res[0].size() && (int)i < show; ++i)
                std::cout << "  " << t.str(res[0][i]) << "   (ed=" << edit_distance(queries[0], t.str(res[0][i])) << ")\n";
        }
    }
    if (runs.size() > 1) std::cout << "\nresult sets identical across modes: " << (mismatch ? "NO" : "yes") << "\n";
    return mismatch ? 1 : 0;
}
