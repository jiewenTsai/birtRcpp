# birtRcpp 待辦（2026-10-07）

姊妹套件 birt 的待辦在 `../birt/TODO.md`。

## 定位

- **主菜：ECM。** 用 `block_model()` 宣告模型，用 `ecm()` 估計（Meng & Rubin, 1993）：
  - 每一塊參數各自做條件極大化；
  - logistic 部分用 Pólya–Gamma 下界；
  - 自適應或 Bock–Aitkin 積分，常態線性的潛在變項精確積分；
  - 標準誤用 Louis 公式。
- **同一組步驟：** `gibbs()`（Gibbs；共軛更新、MALA）、`vi()`（CAVI）。`mml()` 用 Newton 直接極大化同一個目標。
- **檢查：** `check_ecm()`、`check_sampler()`（含 `step =` 手寫抽樣器）、`check_vi()`；`algorithm()` 寫出每一步。
- **內建模型：** `mlirt()`、`rtirt_*()`（含 `quantile =`），以 C++ Gibbs 移植 ExtendedRtIrtModeling.jl；`ecm()` 也可以估計內建的常態模型。
- **不做：** lavaan 語法、JAGS、RTMB（屬於 birt）；HMC（MALA 可以）。
- **保留：** Julia 的 camelCase 別名（設計原則）；和 birt 共用泛型的轉接層。

2026-10-06 的變動：
- `em()` 改名 `ecm()`，`check_em()` 改名 `check_ecm()`，`init = "ecm"`，配適結果裡的 `$em` 改成 `$ecm`。
- Newton／quasi-Newton 從 ECM 分出來成為 `mml()`。
- 積木層改用和內建模型同一組動詞（`ecm()`、`gibbs()`、`vi()`）。
- 其他清理：見 NEWS。

## 1. ECM 主菜要補的

1. [ ] 內建的 `ecm()`（內建模型）和 `ecm.block_model()` 只剩一個引擎：讓內建模型改為宣告式模型的捷徑。
   - 目前內建版快一倍（1000 × 20：0.9 秒對 2.1 秒；ba81 1.7 秒）。
   - 要先補上速度，見 4.1 的 C++ 逐格累加。
2. [x] `quadrature = "ba81"` 配 `analytic = FALSE` 時，在後驗很窄的變數上不準（cross 模型的 log-likelihood 差 9）：自動警告，或改用自適應積分。
   - 2026-10-06：既有的 `adaptive_gap` 檢查會觸發（1000 × 20 cross：ba81 + analytic = FALSE 的 log-likelihood −32543.4，自適應 −32364.4；gap −25.9 → 警告）。警告改為寫出差距並建議 `analytic = TRUE` 或 `quadrature = "adaptive"`；測試已補。
3. [x] Louis 資訊矩陣：有不積掉的 `[person]` 參數時，目前退回數值微分。
   - 2026-10-07：改成解析的 Louis 公式（分數的共變異依人分塊，逐人計算）。和分數的數值微分相差 ≤ 9e-7（相對；標準誤 3e-8）；200 × 6 的配適 1.5 秒 → 0.3 秒，其他配適的估計值與標準誤完全不變。
4. [x] 宣告式分數檢定（Glas & van der Linden, 2010 的 LM 檢定）：在虛無模型的 ECM 估計值上，用 Fisher identity 算分數、用 Louis 資訊算 W；需要「固定為 0 的參數」的語法。
   - 2026-10-07：`score_test(fit, add = Y ~ delta[item] * z)`。500 × 8、100 次：LM 與 LR 相關 0.9995（無 DIF）、0.997–0.999（有 DIF）；虛無下拒絕率 5.0%（整體）、4.7%（逐題），KS p 0.62、0.38；DIF 0.5／0.8 的檢定力 57%／92%；重現 `ci_test()`（差 ≤ 4e-4）與內建 `score_test()` 對 b 的 LM（相關 1.000）。被虛無參數吸收的方向（題目整體平移 + 潛在平均）由逐人分數找出並排除。
5. [ ] 混合 MCEM／SAEM 引擎（規則 R1–R6，`../results/62_fuzz/hybrid_report.md`），尚未獲同意。
6. [x] 教學與 vignette 改成以 `ecm()` 為主線（2026-10-06；入門教學、積木教學、vignette、ex01–ex13）。
7. [x] **ECM 的速度（主菜）：** 2026-10-06 查清楚並加速（i、ii、iii 已做：無起點 44 秒 → 5–7 秒；預設 tol 改為 1e-10）。
   - 同一份資料（1000 × 20 cross）、同一個模型、同樣的起點、`prior = "ml"`：內建 `ecm()` 1.3 秒；宣告式 `ecm()` 自適應 9.7 秒（25 次迭代）；`quadrature = "ba81"` 3.0 秒；`mml()` Newton 20.1 秒，quasi-Newton 19.0 秒。log-likelihood 都是 −32364.08。所以 ECM 比 Newton 快 2–7 倍。
   - 「ECM 53 秒、mml 21 秒」是在 CPU 滿載（load 約 10）、用預設設定時量的。預設 `prior = "map"`、沒給起點時，ECM 需要 56 次迭代（16.3 秒），`mml()` 17.5 秒。
   - 和舊紀錄 2.1 秒的差距不是今天清理造成的：今天早上的備份和現在的程式並排量，時間相同（自適應 8.8–9.2 秒，ba81 2.8–3.0 秒）。舊紀錄的確切條件無法重建。
   - **時間花在哪：** 每一步 ECM 0.144 秒。logistic 部分在「人 × 21 個節點 × 題目」（42 萬格）上評估，a 和 b 各算一次（0.053、0.049 秒），佔約 70%。時間部分在每人 4 個點上，只要 0.007 秒。內建版快，是因為 logistic 部分彙整成共同節點 × 題目（41 × 20）。
   - **可做的加速：**
     - (i) a、b 共用同一次的 eta 與 Pólya–Gamma 權重，或合成一個 (a, b) 區塊：logistic 部分約減半。
     - (ii) 刪掉權重可忽略（< 1e-12）的節點。
     - (iii) 預設改用依資料的起點：迭代數約減半。
     - (iv) 預設積分改成 ba81，並保留和自適應積分的差距檢查。2026-10-06 使用者決定預設不動。

8. [x] ECM 的預設起始值（2026-10-06 已做 `data_start()`；Bolsinova 模型待重測）：宣告式 Bolsinova 條件相依模型從預設起點跑了 8.5 分鐘，停在不好的點（cc = −2.1，log posterior −17739.6）；用依資料的起點 4.5 分鐘到合理解（−16396.2）。考慮預設改用依資料的起點，或多起點。
9. [x] `check_ecm()` 接受內建模型（2026-10-06）：用 em.R 的 em_setup／em_estep／em_mstep／em_score 做三項檢查；mlirt、null、latreg、latent、cross（含 fixed、1PL；ml 與 map）全部通過；故意改錯的 M 步（sigma2t × 1.1、rho + 0.05）與改錯的 score 都被抓到（測試）。
10. [x] LLTM 的 ECM 估計 vb = 0.054（真值 0.09）：教學以「右偏後驗的眾數 vs 平均」解釋，沒有引用；查證或改用模擬說明。
   - 2026-10-06：40 次模擬（教學的設計）：ECM 平均 0.065，Gibbs 後驗眾數 0.086、平均 0.110。原解釋不對：ECM 是 (b, vb) 的聯合眾數（b 沒積掉；O'Hagan, 1976），加上 3 個係數用掉自由度（Harville, 1977）。教學已改寫並重新產生 html。

## 2. 積木層的擴充

1. [ ] **具名維度**（類似 PyMC 的 `dims=`）：`tau[person, testlet]`、`d[item, category]`。需要改形狀檢查、加總規則、C++ kernel、ECM／CAVI 的 E 步；是大重構，要先規劃。
2. [ ] **缺的分配與積木：**
   - probit（Albert & Chib, 1993）、多元常態區塊。
   - 混合模型：首要目標是 Wang & Xu (2015, BJMSP 68, 456–477)，scripts/66 有原型。
3. [ ] **Pólya–Gamma 擴展到其他資料**（需要 2.1 的具名維度）：
   - 二項、負二項（Polson et al., 2013；Zhou et al., 2012）；
   - 多類別（Linderman et al., 2015）；
   - GPCM 的相鄰類別 logit（要查證能否保持條件共軛）。
4. [ ] `vi()`：沒有 MALA 的對應（候選：Laplace／NCVMP）。
   - [x] 2026-10-07：`quadrature = "ba81"` 與 2L 點（含 ECM 的 affine 精確步驟）。同一個不動點（1e-13），ELBO 不下降；1000 × 20：203 秒 → 56 秒（ba81 21 秒；load 1.7–2.9）。
5. [x] `algorithm()` 的 Gibbs 與 ECM 改成和 CAVI 一樣的 Unicode／LaTeX 格式。
   - 2026-10-07：模型 p，再列每一步的完整條件分配（Gibbs，寫出 p 與 m）或條件眾數（ECM：目標值、E 步、M 步的期望值）；`format = "latex"`。順便修正運算式少掉的括號（CAVI 的殘差）。

## 3. 內部整理（2026-10-06 稽核，未處理）

0. [x] ex03：q = 0.9 的 `rtirt_cross` 在偏態資料上，各鏈停在不同的眾數（R-hat 22–29；`collapse = TRUE`、換 seed 都一樣）。
   - 2026-10-06 查明：ALD 工作後驗真的多峰，不是抽樣器或起點的問題。θ 變成第二個速度因子（a ≈ 0.05–0.37，|rho| ≈ 0.55，各題正負號不同）；各鏈停在不同的正負號組合，邊際 ALD log-likelihood 相近，比 rho ≈ 0 的解高約 190。從 ECM 起點開始也一樣；20000 次迭代也不換峰。只有 q = 0.9 如此（q = 0.1–0.75 的 R-hat ≤ 1.04）。ex03 的註解已改。

1. [x] **重複的程式**（2026-10-06）：
   - 速度常態積分：只剩 `normal_given_sums()`（em_estep、marg_ll_draw、calc_normal_given_nodes 共用）。
   - 2PL log-likelihood：只剩 `calc_loglik_nodes()`（softplus）；`marginal_loglik()`、`dic()` 在 eta ≈ 2400 時也有限（測試）。
   - 節點只用 `calc_nodes()`；`lse_cols` 刪掉（用 `calc_node_weights()`）；R 版 `log_ald_rt` 移到測試當參考。
   - `pg_mean()` 改成 C++（向量化，保留維度），C++ 只剩一個 `pg_mean1()`。
   - Fisher-scoring MH 只剩 `item_mh_step()`（form 2 = (a, d)），`draw_ab_joint` 與 `.blk_items_mh` 都呼叫它；新舊抽樣值相同到 2e-15，同 seed 的 Geweke p 值新舊一致；cross 模型 16 seeds × 3 設定的 Geweke 每種設定至多 1 次未過（Holm 0.01，與預設抽樣器相同），SBC 通過。
   - 格式化：`fmt_num()`／`fmt_trunc()` 取代 al_*／bt_*；`var_text` 改用 `bt_var`；所有 algorithm()／equations() 文字逐字相同。
   - `runSimulation()` 呼叫 `run_simulation()`（新增 `speed_var`），結果與舊版相同（identical）。
   - 數值：ECM 估計值改變 ≤ 6e-8（z ≤ 4e-7，log-likelihood ≤ 2e-12），在收斂容許度（參數 1e-5）之內；原因是速度積分的捨入順序不同。ALD 模型 `collapse = TRUE` 的鏈因最後一位的差異在數十次迭代後分開（同一個抽樣器，不同的實現）。
2. [x] 文件：WORDLIST 補新詞（101 → 244 詞，`spelling::spell_check_package()` 無錯；behaviour → behavior）；ex07 已無 `:::`；積木教學 html 已重新產生。
3. [x] 已決定保留（2026-10-06）：`dic`、`calc_artificial_data`、`calc_normal_given_nodes` 繼續匯出；`intercept`／`itemtype` 可以在建構子和 `gibbs()` 兩處設定。

## 4. 速度（1000 × 20，cross 模型）

| | 時間 |
|---|---|
| `ecm()`（宣告式） | 2.1 秒（ba81：1.7 秒） |
| `ecm()`（內建） | 0.9 秒 |
| `gibbs()`（宣告式） | 23 秒（內建 11 秒；Pólya–Gamma 抽樣約佔 44%） |
| `vi()` | 56 秒（ba81：21 秒；2026-10-07，之前 203 秒） |

1. [ ] C++ 逐格累加（方案 A）：ECM 主菜的速度。
2. [ ] `vi()` 的候選做法：
   - 對 CAVI 套 SQUAREM（要查證 ELBO 的單調性）；
   - ~~改用 ECM 的 2L 點與 ba81 彙整~~（2026-10-07 已做：203 → 56 秒，ba81 21 秒）；
   - ~~把 `pg_mean()` 改成 C++~~（2026-10-06 已做）。

## 5. 研究線（以 birtRcpp 的工具為主）

- ECM 的理論出處：`../results/62_fuzz/ecm_theory_sources.md`。兩點目前查不到直接出處，是我們自己推論的：
  - 速度精確積分是我們的化簡；
  - 固定中心使 EM 單調，是依 Dempster et al. (1977) Theorem 1 推得。
- 分位數的 Gibbs（ALD）：`rtirt_latent()`／`rtirt_cross()` 的 `quantile =`，作為 birt 研究線 5a 的第 1、2、8 項的對照。
- Zotero 建議補上：Wang & Xu (2015)、Ulitzsch, von Davier & Pohl (2020)、Molenaar & De Boeck (2018)，以及 `ecm_theory_sources.md` 末尾列的 15 篇。
