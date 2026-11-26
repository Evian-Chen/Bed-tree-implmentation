// bptree_tool：單獨測試自製 B+-tree (key 為整數，方便觀察)
//   bptree_tool insert 5 3 9 1 7 8 2 --leaf-cap 3 --internal-cap 3   依序插入並印出樹
//   bptree_tool erase  5 3 9 1 7 --leaf-cap 3 --internal-cap 3 --del 3 9
//   bptree_tool bulk 20 --leaf-cap 4 --internal-cap 3                 bulk load 1..20
//   bptree_tool random --count 20000 --leaf-cap 8 --internal-cap 6    隨機插入/刪除並檢查不變量
#include <iostream>
#include <random>
#include <set>

#include "BPlusTree.h"
#include "Cli.h"

using namespace bed;

struct IntEntry {
    ZKey key;
    uint32_t id;
    IntEntry() : id(0) {}
};

static IntEntry mk(uint64_t v) {
    IntEntry e;
    e.key.w[kKeyWords - 1] = v;
    e.id = (uint32_t)v;
    return e;
}

typedef BPlusTree<IntEntry> Tree;

static void dump(const Tree::Node* n, int depth) {
    std::string ind(depth * 4, ' ');
    if (n->leaf) {
        std::cout << ind << "leaf:";
        for (size_t i = 0; i < n->entries.size(); ++i) std::cout << " " << n->entries[i].key.w[kKeyWords - 1];
        std::cout << "\n";
        return;
    }
    std::cout << ind << "internal:\n";
    for (size_t i = 0; i < n->kids.size(); ++i) {
        std::cout << ind << " [" << n->kids[i].lo.key.w[kKeyWords - 1] << ".." << n->kids[i].hi.key.w[kKeyWords - 1] << "]\n";
        dump(n->kids[i].node.get(), depth + 1);
    }
}

static void report(const Tree& t) {
    dump(t.root(), 0);
    std::string err;
    bool ok = t.check_invariants(&err);
    size_t lv, in;
    t.count_nodes(&lv, &in);
    std::cout << "size=" << t.size() << " height=" << t.height() << " leaves=" << lv << " internals=" << in
              << "  invariants: " << (ok ? "OK" : ("FAILED: " + err)) << "\n";
}

int main(int argc, char** argv) {
    Args a(argc, argv);
    if (a.pos.empty()) {
        std::cerr << "usage: bptree_tool insert|erase|bulk|random ...\n";
        return 2;
    }
    Tree t((size_t)a.geti("leaf-cap", 4), (size_t)a.geti("internal-cap", 4));
    const std::string cmd = a.pos[0];
    if (cmd == "insert" || cmd == "erase") {
        for (size_t i = 1; i < a.pos.size(); ++i) {
            t.insert(mk((uint64_t)std::atoll(a.pos[i].c_str())));
            if (a.has("verbose")) { std::cout << "-- after insert " << a.pos[i] << "\n"; dump(t.root(), 0); }
        }
        if (cmd == "erase") {
            std::vector<int> del = parse_int_list(a.get("del"));
            for (size_t i = 0; i < del.size(); ++i) {
                bool ok = t.erase(KeyId(mk((uint64_t)del[i]).key, (uint32_t)del[i]));
                std::cout << "erase " << del[i] << (ok ? ": removed\n" : ": not found\n");
            }
        }
        report(t);
    } else if (cmd == "bulk" && a.pos.size() >= 2) {
        int n = std::atoi(a.pos[1].c_str());
        std::vector<IntEntry> es;
        for (int i = 1; i <= n; ++i) es.push_back(mk((uint64_t)i));
        t.bulk_load(es, a.getd("fill", 0.9));
        report(t);
    } else if (cmd == "random") {
        int count = a.geti("count", 10000);
        std::mt19937_64 rng((uint64_t)a.geti("seed", 1));
        std::set<uint64_t> ref;
        int bad = 0;
        for (int i = 0; i < count; ++i) {
            uint64_t v = rng() % (uint64_t)(count * 2);
            if (ref.count(v) || rng() % 3) {
                if (!ref.count(v)) { t.insert(mk(v)); ref.insert(v); }
                else { t.erase(KeyId(mk(v).key, (uint32_t)v)); ref.erase(v); }
            } else {
                t.insert(mk(v)); ref.insert(v);
            }
            if (i % 1000 == 0) {
                std::string err;
                if (!t.check_invariants(&err)) { std::cout << "invariant failed at op " << i << ": " << err << "\n"; ++bad; }
            }
        }
        std::vector<uint64_t> got;
        struct C { std::vector<uint64_t>* g; void operator()(const IntEntry& e) { g->push_back(e.key.w[kKeyWords - 1]); } } c = {&got};
        t.for_each(c);
        bool same = got.size() == ref.size() && std::equal(got.begin(), got.end(), ref.begin());
        std::cout << "content matches std::set: " << (same ? "yes" : "NO") << "\n";
        std::string err;
        bool ok = t.check_invariants(&err);
        size_t lv, in;
        t.count_nodes(&lv, &in);
        std::cout << "final size=" << t.size() << " height=" << t.height() << " leaves=" << lv << " internals=" << in
                  << " invariants=" << (ok ? "OK" : err) << "\n";
        return (same && ok && !bad) ? 0 : 1;
    } else {
        std::cerr << "bad command\n";
        return 2;
    }
    return 0;
}
