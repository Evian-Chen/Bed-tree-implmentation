# Bed-tree 與 Bed⁺-tree（two-resolution gram signatures）實作

- **Bed-tree（baseline）**：B+-tree 的 key 是「L 維 n-gram 計數向量的 z-order」，內部節點用 key 區間推出的上下界算編輯距離下界來剪枝。
- **Bed⁺-tree（proposal）**：排序仍用低維 key（L=4），但每個 leaf entry 額外存一個高維 signature（S=64 個 bucket 的計數），進入昂貴的 DP 之前先用它算下界過濾。
- 兩者共用同一份程式（`BedTree` 類別），差別只有 `sig` 參數：`sig=0` 是 Bed-tree，`sig>0` 是 Bed⁺-tree。

實驗設計、預期與實際結果見 [docs/EXPERIMENTS.md](docs/EXPERIMENTS.md)。

## 使用的資料庫系統與資料集

- **資料庫系統：沒有使用任何資料庫系統（無 SQLite / PostgreSQL / MySQL 等）。** B+-tree 是自己寫的、純記憶體內的結構（[include/BPlusTree.h](include/BPlusTree.h)）。沒有磁碟檔案、沒有 buffer pool；「I/O」只是用「讀到的節點數（假設一個節點 = 一個 4096-byte page）」來估計。
- **資料集：目前只有合成資料與手寫小樣本，沒有用論文的真實資料集。**

| 檔案 | 內容 | 用途 |
|---|---|---|
| `data/sample_words.txt` | 手寫 61 個英文單字 / 人名 / 常見拼錯字 | 教學、手動測試 |
| `data/sample_queries.txt` | 7 個手寫 query | 教學 |
| `data/sample_pairs.txt` | 13 組字串對與正確編輯距離 | 測 `ed_tool` |
| `data/words_syn.txt` | **合成**單字 370,105 筆（`gen_data words --count 370105 --seed 1`，平均長度 13.1，最長 20） | 全部正式實驗 |

`words_syn.txt` 是用「音節」隨機組成的假單字，筆數刻意對齊 proposal 提到的 370,105，但**不是真實詞庫**，分佈與真實資料不同（見實驗報告的限制）。要換成真實資料：準備一個「每行一個字串」的文字檔，所有工具的 `--data` 都可以直接用。

## 環境與建置

需要 g++（C++14）。我在 Windows 11 + MinGW g++ 6.3.0 上開發與驗證；該版本不支援完整 C++17，所以程式用 C++14。**不需要 make 或 cmake。**

Git Bash / Linux：

```bash
bash scripts/build.sh
```

PowerShell：

```powershell
.\scripts\build.ps1
```

執行檔會放在 `build/`（Windows 上是 `.exe`，以下指令省略副檔名即可在 Git Bash 執行；PowerShell 用 `.\build\xxx.exe`）。

`CMakeLists.txt` 也有提供，但**這台機器沒有 make / MSVC，所以 CMake 沒有實際驗證過**，請以建置腳本為準。

## 每個小功能的單獨執行方式

以下範例輸出都是實際執行過的結果（雜湊值依本專案固定的 FNV-1a 演算法，跨平台一致）。

### 1. 編輯距離 `ed_tool`

```bash
./build/ed_tool kitten sitting                 # 精確編輯距離 = 3
./build/ed_tool kitten sitting --tau 2         # 帶門檻驗證：超過門檻回傳 tau+1，並提示 ed > 2
./build/ed_tool kitten sitting --table         # 印出完整 DP 矩陣
./build/ed_tool --file data/sample_pairs.txt   # 批次對照預期值，最後印 13/13 passed
```

`--file` 每行格式：`字串1<TAB>字串2<TAB>預期距離`，`#` 開頭為註解。這個模式同時檢查「帶門檻的 `verify_ed`」在 tau=0..6 下與精確版一致。

### 2. n-gram 與 bucket 向量 `gram_tool`

```bash
./build/gram_tool "Jim Gray"                          # n=2, L=4
./build/gram_tool "Jim Gray" --n 3 --L 8 --sig 64     # 同時顯示 64 維 signature
```

會列出每個 n-gram（頭尾補 `#`，共 |s|+n-1 個）、它落在哪個 bucket、L 維計數向量、z-order key，以及 signature 的非零 bucket。signature 的總和必等於 n-gram 數。

### 3. z-order 與下界 `zorder_tool`

```bash
./build/zorder_tool encode --B 2 3 2 1 3              # -> 11011011（論文 Figure 2）
./build/zorder_tool decode --L 4 --B 2 11011011       # -> 3 2 1 3
./build/zorder_tool bounds --L 4 --B 2 11011011 11011110   # 每個 bucket 的 [lb,ub]（論文 Figure 3）
./build/zorder_tool lb --q "Jam Gray" --s1 "Jim Gray" --s2 "Jim Grey" --L 4
```

`bounds` 的輸出中，第 1 個 bucket 為 `[3, 3]`，與論文一致。若上界等於該維能表示的最大值，會標示 saturated：代表「≥ 該值」，下界計算會把它當作無上界，以保證不會漏答案。
`lb` 會印出 query 到「以 s1、s2 的 key 為端點的區間」的下界，並列出兩端字串的真實距離作對照（下界必須 ≤ 區間內任何字串的真實距離）。

### 4. B+-tree `bptree_tool`（key 用整數，方便觀察）

```bash
./build/bptree_tool insert 5 3 9 1 7 8 2 6 4 --leaf-cap 3 --internal-cap 3       # 依序插入並印出樹與不變量檢查
./build/bptree_tool erase 5 3 9 1 7 --leaf-cap 3 --internal-cap 3 --del 3,9       # 刪除
./build/bptree_tool bulk 20 --leaf-cap 4 --internal-cap 3                         # bulk load 1..20
./build/bptree_tool random --count 20000 --leaf-cap 8 --internal-cap 6            # 隨機插入/刪除，對照 std::set
```

`random` 會每 1000 次操作檢查一次結構不變量（排序、容量、各子樹 lo/hi 正確、所有 leaf 同深度），最後比對內容與 `std::set` 是否一致。

### 5. 資料產生 `gen_data`

```bash
./build/gen_data words  --count 370105 --seed 1 --out data/words_syn.txt
./build/gen_data names  --count 100000 --seed 1 --out data/names_syn.txt
./build/gen_data titles --count 50000  --seed 1 --out data/titles_syn.txt
./build/gen_data queries --data data/sample_words.txt --count 20 --edits 2 --seed 7 --out data/q_demo.txt
```

`queries` 會從資料集中抽字串並隨機編輯 `--edits` 次。固定 seed 可重現。

### 6. 建索引並查詢 `bedtree_cli`

```bash
# 三種方法同時跑並比對結果是否一致
./build/bedtree_cli --data data/sample_words.txt --q recieve --tau 2 --compare

# 只跑一種：scan / bed / bedplus
./build/bedtree_cli --data data/sample_words.txt --q "Jon Smyth" --tau 2 --mode bedplus --L 4 --sig 64

# 一次跑一個 query 檔
./build/bedtree_cli --data data/words_syn.txt --queries data/sample_queries.txt --tau 1 --compare
```

選項：`--n` gram 長度、`--L` 排序維度、`--sig` signature 維度、`--page` 模擬 page bytes、`--show` 最多列出幾筆結果。輸出的統計欄位：

| 欄位 | 意義 |
|---|---|
| `nodes` (internal / leaf) | 讀到的節點數，視為 page 讀取數 |
| `entries_scanned` | 在讀到的 leaf 內檢查的 entry 數 |
| `pruned_key` | 被 L 維 key 向量下界濾掉的數量 |
| `pruned_sig` | 被 signature 下界濾掉的數量（只有 Bed⁺） |
| `len_skip` | 長度差 > tau 直接略過 |
| `DP` | 真正跑 DP 驗證的次數 |

### 7. 測試 `test_all`

```bash
./build/test_all            # 全部 (約 6 秒，目前 220,000+ 個檢查)
./build/test_all ed         # 只測編輯距離
./build/test_all gram       # n-gram / bucket 向量
./build/test_all zorder     # z-order 編碼、解碼、區間界
./build/test_all bptree     # B+-tree 對照 std::set 與不變量
./build/test_all bounds     # 下界永遠 <= 真實編輯距離
./build/test_all bedtree    # Bed-tree / Bed⁺-tree 與暴力搜尋結果完全相同 (含插入刪除後)
```

最重要的是 `bounds` 與 `bedtree`：前者驗證所有下界都合法，後者驗證端到端不漏、不多答案，涵蓋多種 (n, L, signature 維度, page 大小, 字串長度) 組合，包括字串很長造成計數飽和的情況。
我也做過變異測試：故意把下界公式改壞，這兩組測試會大量失敗，所以它們確實有抓錯能力。

### 8. 實驗 `bench`

```bash
./build/bench e2 --data data/words_syn.txt --qcount 100 --reps 2      # 主要比較
./build/bench all --data data/words_syn.txt --qcount 100 --reps 2     # 全部 (約 8 分鐘)
./build/bench e1 --data data/words_syn.txt --limit 50000 --qcount 50  # 小資料快速試跑
```

子命令 `e1`…`e6` 說明見 [docs/EXPERIMENTS.md](docs/EXPERIMENTS.md)。CSV 會寫到 `--out` 指定的資料夾（預設 `results/`，需已存在）。**每個實驗都會對照暴力搜尋檢查結果**，結尾會印 `TOTAL queries whose result differed from brute force`，不是 0 時程式回傳非 0。

## 程式結構

```
include/EditDistance.h   編輯距離：完整 DP 與帶門檻的 band DP (verify_ed)
include/StringMapper.h   n-gram、64-bit hash、L 維 / S 維 bucket 向量 (src/StringMapper.cpp)
include/ZOrder.h         z-order 編碼/解碼、key 區間的每維上下界、Eq.(4) 下界
include/BPlusTree.h      自製 B+-tree：insert / erase / bulk load / 不變量檢查
include/BedTree.h        BedTree 類別：sig=0 → Bed-tree，sig>0 → Bed+-tree；range_query / scan_query
include/Cli.h            參數解析、讀檔、計時 (工具共用)
src/tools/*.cpp          各個命令列工具與 bench
tests/test_all.cpp       單元測試與暴力搜尋對照
scripts/                 建置腳本
```

## 重要實作細節

- **下界**：一次編輯最多改變 n 個 n-gram，所以 `ED ≥ ceil( max(Σ_{q>ub}(q−ub), Σ_{lb>q}(lb−q)) / n )`。
- **z-order key**：最多 256 bit；每維 bit 數預設 `min(8, 256/L)`（L=64 時為 4 bit）。計數超過上限時「飽和」，飽和維度的上界視為無限大，才不會漏答案。
- **兩組向量用同一個 hash**：L 整除 S 時，L 維向量就是 S 維 signature 的摺疊。
- **signature 存成 uint8**：字串的 n-gram 數 > 255 時該 entry 標記為不可用（跳過 signature 過濾），避免飽和造成誤判。
- **節點容量**依 page 大小與 entry 大小計算（entry = key bytes + id + 估計字串長度 + signature bytes），所以 Bed⁺ 的 leaf 容量較小、頁數較多，這就是 proposal 提到的額外儲存成本。
- **刪除是 lazy 的**：只移除空節點，不做 merge / redistribute。
- 驗證步驟 `verify_ed` 與暴力搜尋都先做長度差檢查（長度差 > tau 直接排除），三種方法一致，所以比較是公平的，但這也讓 scan 比「純 DP」快很多。
