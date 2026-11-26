// 編輯距離 (Levenshtein: insert / delete / substitute 各 cost 1)
#ifndef BED_EDIT_DISTANCE_H
#define BED_EDIT_DISTANCE_H

#include <algorithm>
#include <cstdlib>
#include <string>
#include <vector>

namespace bed {

// 標準 O(|a|*|b|) DP，回傳精確編輯距離。用於測試與小資料。
inline int edit_distance(const std::string& a, const std::string& b) {
    const size_t n = a.size(), m = b.size();
    std::vector<int> prev(m + 1), cur(m + 1);
    for (size_t j = 0; j <= m; ++j) prev[j] = (int)j;
    for (size_t i = 1; i <= n; ++i) {
        cur[0] = (int)i;
        for (size_t j = 1; j <= m; ++j) {
            int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            cur[j] = std::min(std::min(prev[j] + 1, cur[j - 1] + 1), prev[j - 1] + cost);
        }
        prev.swap(cur);
    }
    return prev[m];
}

// 回傳完整 DP 矩陣 (教學用，ed_tool --table 使用)
inline std::vector<std::vector<int> > edit_distance_matrix(const std::string& a,
                                                           const std::string& b) {
    std::vector<std::vector<int> > d(a.size() + 1, std::vector<int>(b.size() + 1, 0));
    for (size_t i = 0; i <= a.size(); ++i) d[i][0] = (int)i;
    for (size_t j = 0; j <= b.size(); ++j) d[0][j] = (int)j;
    for (size_t i = 1; i <= a.size(); ++i)
        for (size_t j = 1; j <= b.size(); ++j) {
            int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            d[i][j] = std::min(std::min(d[i - 1][j] + 1, d[i][j - 1] + 1), d[i - 1][j - 1] + cost);
        }
    return d;
}

// 帶門檻的驗證 (論文 Algorithm 1 的精神)：只算寬度 2*tau+1 的對角帶，
// 整列最小值超過 tau 就提早結束。
// 回傳：distance <= tau 時為精確距離，否則回傳 tau+1。
inline int verify_ed(const std::string& a, const std::string& b, int tau) {
    const int n = (int)a.size(), m = (int)b.size();
    if (std::abs(n - m) > tau) return tau + 1;
    if (n == 0) return m;  // 此時 m <= tau
    if (m == 0) return n;
    const int INF = tau + 1;
    static std::vector<int> prev, cur;  // 單執行緒使用，重複利用緩衝區
    prev.assign(m + 2, INF);
    cur.assign(m + 2, INF);
    for (int j = 0; j <= std::min(m, tau); ++j) prev[j] = j;
    for (int i = 1; i <= n; ++i) {
        const int lo = std::max(1, i - tau);
        const int hi = std::min(m, i + tau);
        int rowmin = INF;
        if (lo == 1) {
            cur[0] = (i <= tau) ? i : INF;
            rowmin = cur[0];
        } else {
            cur[lo - 1] = INF;
        }
        for (int j = lo; j <= hi; ++j) {
            int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            int v = std::min(std::min(prev[j] + 1, cur[j - 1] + 1), prev[j - 1] + cost);
            if (v > INF) v = INF;
            cur[j] = v;
            if (v < rowmin) rowmin = v;
        }
        if (hi + 1 <= m) cur[hi + 1] = INF;
        if (rowmin > tau) return tau + 1;
        prev.swap(cur);
    }
    return prev[m] <= tau ? prev[m] : tau + 1;
}

}  // namespace bed

#endif
