// gen_data：產生「合成」測試資料 (沒有網路資料集時使用；固定 seed 可重現)
//   gen_data words  --count 370105 --seed 1 --out data/words_syn.txt     英文風格單字，平均約 9 字元
//   gen_data names  --count 100000 --seed 1 --out data/names_syn.txt     "First Last" 人名風格
//   gen_data titles --count 50000  --seed 1 --out data/titles_syn.txt    多單字標題，平均約 40 字元
//   gen_data queries --data F --count 200 --edits 2 --seed 7 --out Q     由資料集中抽字串，隨機編輯 k 次
// 注意：這是以「音節」隨機組成的假資料，不是真實詞庫；分佈與真實資料不同 (見 README)。
#include <fstream>
#include <iostream>
#include <random>
#include <set>

#include "Cli.h"

using namespace bed;

static const char* kOnset[] = {"b","c","d","f","g","h","j","k","l","m","n","p","r","s","t","v","w","z",
                               "br","cr","dr","fr","gr","pr","tr","bl","cl","fl","gl","pl","sl","st","sp","sh","ch","th","qu",""};
static const char* kNuc[] = {"a","e","i","o","u","ai","ea","ou","io","ee","oo","ia"};
static const char* kCoda[] = {"","","","n","r","s","t","l","m","d","ng","st","nd","rt","ck"};

template <class T, size_t N>
static const char* pick(std::mt19937_64& g, T (&arr)[N]) { return arr[g() % N]; }

static std::string make_word(std::mt19937_64& g, int min_syl, int max_syl) {
    int syl = min_syl + (int)(g() % (uint64_t)(max_syl - min_syl + 1));
    std::string w;
    for (int i = 0; i < syl; ++i) {
        w += pick(g, kOnset);
        w += pick(g, kNuc);
        w += pick(g, kCoda);
    }
    return w;
}

static std::string cap(std::string s) {
    if (!s.empty() && s[0] >= 'a' && s[0] <= 'z') s[0] = (char)(s[0] - 'a' + 'A');
    return s;
}

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

int main(int argc, char** argv) {
    Args a(argc, argv);
    if (a.pos.empty() || !a.has("out")) {
        std::cerr << "usage: gen_data words|names|titles|queries --out FILE [--count N] [--seed S] ...\n";
        return 2;
    }
    const std::string kind = a.pos[0];
    size_t count = (size_t)a.geti("count", 1000);
    std::mt19937_64 g((uint64_t)a.geti("seed", 1));
    std::vector<std::string> out;

    if (kind == "queries") {
        std::vector<std::string> data = load_lines(a.get("data"));
        int edits = a.geti("edits", 1);
        for (size_t i = 0; i < count; ++i) out.push_back(mutate(data[g() % data.size()], edits, g));
    } else {
        std::set<std::string> seen;
        while (out.size() < count) {
            std::string s;
            if (kind == "words") {
                s = make_word(g, 1, 5);
                if (s.size() < 3 || s.size() > 20) continue;
            } else if (kind == "names") {
                s = cap(make_word(g, 1, 3)) + " " + cap(make_word(g, 1, 3));
            } else if (kind == "titles") {
                int nw = 3 + (int)(g() % 6);
                for (int i = 0; i < nw; ++i) { if (i) s += " "; s += make_word(g, 1, 3); }
            } else {
                std::cerr << "unknown kind\n";
                return 2;
            }
            if (seen.insert(s).second) out.push_back(s);
        }
    }
    std::ofstream f(a.get("out").c_str());
    if (!f) { std::cerr << "cannot write " << a.get("out") << "\n"; return 2; }
    double total = 0;
    for (size_t i = 0; i < out.size(); ++i) { f << out[i] << "\n"; total += (double)out[i].size(); }
    std::cout << "wrote " << out.size() << " strings to " << a.get("out") << " (avg length "
              << (out.empty() ? 0.0 : total / (double)out.size()) << ")\n";
    return 0;
}
