// ed_tool：單獨測試編輯距離
//   ed_tool "kitten" "sitting"               精確編輯距離
//   ed_tool "kitten" "sitting" --tau 2       帶門檻驗證 (verify_ed)
//   ed_tool "kitten" "sitting" --table       印出 DP 矩陣
//   ed_tool --file data/sample_pairs.txt     批次檢查 (每行: s1<TAB>s2<TAB>預期距離)
#include <iostream>
#include <sstream>

#include "Cli.h"
#include "EditDistance.h"

using namespace bed;

int main(int argc, char** argv) {
    Args a(argc, argv);
    if (a.has("file")) {
        std::vector<std::string> lines = load_lines(a.get("file"));
        int bad = 0, total = 0;
        for (size_t i = 0; i < lines.size(); ++i) {
            if (lines[i][0] == '#') continue;
            std::vector<std::string> f;
            std::stringstream ss(lines[i]);
            std::string t;
            while (std::getline(ss, t, '\t')) f.push_back(t);
            if (f.size() < 3) continue;
            int expect = std::atoi(f[2].c_str());
            int d = edit_distance(f[0], f[1]);
            // 同時驗證帶門檻版本在各種 tau 下與精確版一致
            bool ok = (d == expect);
            for (int tau = 0; tau <= 6; ++tau) {
                int v = verify_ed(f[0], f[1], tau);
                if (v != (d <= tau ? d : tau + 1)) ok = false;
            }
            ++total;
            if (!ok) ++bad;
            std::cout << (ok ? "PASS " : "FAIL ") << "\"" << f[0] << "\" vs \"" << f[1]
                      << "\"  ed=" << d << " expected=" << expect << "\n";
        }
        std::cout << (total - bad) << "/" << total << " passed\n";
        return bad ? 1 : 0;
    }
    if (a.pos.size() < 2) {
        std::cerr << "usage: ed_tool <s1> <s2> [--tau T] [--table] | ed_tool --file pairs.txt\n";
        return 2;
    }
    const std::string& s1 = a.pos[0];
    const std::string& s2 = a.pos[1];
    std::cout << "s1 = \"" << s1 << "\"\ns2 = \"" << s2 << "\"\n";
    std::cout << "edit_distance = " << edit_distance(s1, s2) << "\n";
    if (a.has("tau")) {
        int tau = a.geti("tau", 1);
        int v = verify_ed(s1, s2, tau);
        if (v <= tau) std::cout << "verify_ed(tau=" << tau << ") = " << v << "  (within threshold)\n";
        else std::cout << "verify_ed(tau=" << tau << ") = " << v << "  (exceeds threshold, i.e. ed > " << tau << ")\n";
    }
    if (a.has("table")) {
        std::vector<std::vector<int> > d = edit_distance_matrix(s1, s2);
        std::cout << "\nDP matrix (rows: s1 prefix, cols: s2 prefix)\n      ";
        std::cout << "  -";
        for (size_t j = 0; j < s2.size(); ++j) std::cout << "  " << s2[j];
        std::cout << "\n";
        for (size_t i = 0; i < d.size(); ++i) {
            std::cout << "  " << (i == 0 ? std::string("-") : std::string(1, s1[i - 1])) << " ";
            for (size_t j = 0; j < d[i].size(); ++j) {
                std::string x = itos(d[i][j]);
                std::cout << (x.size() < 2 ? " " : "") << " " << x;
            }
            std::cout << "\n";
        }
    }
    return 0;
}
