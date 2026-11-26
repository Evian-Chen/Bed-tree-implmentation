// 命令列小工具共用：參數解析、讀資料檔、計時
#ifndef BED_CLI_H
#define BED_CLI_H

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace bed {

// 用法：--key value 或 --flag；其餘為位置參數 (pos)
struct Args {
    std::map<std::string, std::string> kv;
    std::vector<std::string> pos;
    Args(int argc, char** argv) {
        for (int i = 1; i < argc; ++i) {
            std::string a = argv[i];
            if (a.size() > 2 && a[0] == '-' && a[1] == '-') {
                std::string key = a.substr(2);
                if (i + 1 < argc && std::string(argv[i + 1]).compare(0, 2, "--") != 0) kv[key] = argv[++i];
                else kv[key] = "1";
            } else {
                pos.push_back(a);
            }
        }
    }
    bool has(const std::string& k) const { return kv.count(k) > 0; }
    std::string get(const std::string& k, const std::string& def = "") const {
        std::map<std::string, std::string>::const_iterator it = kv.find(k);
        return it == kv.end() ? def : it->second;
    }
    int geti(const std::string& k, int def) const { return has(k) ? std::atoi(get(k).c_str()) : def; }
    double getd(const std::string& k, double def) const { return has(k) ? std::atof(get(k).c_str()) : def; }
};

// 每行一個字串；去掉 \r；略過空行
inline std::vector<std::string> load_lines(const std::string& path) {
    std::vector<std::string> out;
    std::ifstream f(path.c_str());
    if (!f) { std::cerr << "cannot open file: " << path << "\n"; std::exit(2); }
    std::string line;
    while (std::getline(f, line)) {
        while (!line.empty() && (line.back() == '\r' || line.back() == '\n')) line.pop_back();
        if (!line.empty()) out.push_back(line);
    }
    return out;
}

struct Timer {
    std::chrono::steady_clock::time_point t0;
    Timer() : t0(std::chrono::steady_clock::now()) {}
    double seconds() const {
        return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    }
};

inline std::string itos(long long v) {  // 舊版 MinGW 的 printf 不支援 %lld，也沒有 std::to_string
    if (v == 0) return "0";
    bool neg = v < 0;
    unsigned long long u = neg ? (unsigned long long)(-(v + 1)) + 1ULL : (unsigned long long)v;
    std::string r;
    while (u) { r.insert(r.begin(), (char)('0' + (int)(u % 10))); u /= 10; }
    return neg ? "-" + r : r;
}

inline std::vector<int> parse_int_list(const std::string& s) {
    std::vector<int> v;
    size_t i = 0;
    while (i < s.size()) {
        size_t j = s.find(',', i);
        if (j == std::string::npos) j = s.size();
        v.push_back(std::atoi(s.substr(i, j - i).c_str()));
        i = j + 1;
    }
    return v;
}

}  // namespace bed

#endif
