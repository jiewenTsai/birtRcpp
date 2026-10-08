> English version: [README.md](README.md)

# birtRcpp

**BIvariate Response Time Modeling with 'Rcpp'**

從完整條件分配寫出估計法的套件：

- **積木層（主軸）：** 用 `block_model()` 宣告模型。主菜是 `ecm()`（ECM：每一塊參數各自做條件極大化；Pólya–Gamma 下界、自適應或 Bock–Aitkin 積分、常態線性潛在變項精確積分、Louis 標準誤）。同一組步驟也給 `gibbs()`（Gibbs；共軛更新、MALA）和 `vi()`（CAVI）；`mml()` 用 Newton 直接極大化同一個目標。每個引擎都有檢查（`check_ecm()`、`check_sampler()`、`check_vi()`），`algorithm()` 寫出每一步。
- **內建模型：** ExtendedRtIrtModeling.jl（Tsai, 2025 博士論文第 3–5 章）的 R 版本，C++（RcppArmadillo）寫的 Gibbs 抽樣器，每個模型一個迴圈；Julia 的函數名稱都保留為別名。
- **分工：** lavaan 語法、JAGS、RTMB 的自動化流程在姊妹套件 birt。

教學文件有兩份（原始檔在 `inst/tutorial/`）：

- `birtRcpp_blocks.qmd`：積木，自己宣告模型、組抽樣器與 EM；
- `birtRcpp_intro.qmd`：用內建模型分析資料。

## 用法

函數、物件、欄位一律用 snake_case；Julia 的原名都保留為別名（沿用 Julia 的參數順序與預設值）。

```r
library(birtRcpp)
raw <- read.csv("timss2019math.csv")
items <- names(raw)[6:19]; times <- paste0(items, "_S"); covs <- names(raw)[44:53]

dat <- input_data(resp = items, time = times, cov = covs, id = "IDSTUD", data = raw)
fit <- gibbs(rtirt_cross(dat, quantile = 0.25), n_iter = 10000, n_chain = 3, seed = 1)

fit                    # 摘要：信度、DIC、R-hat 警告
summary(fit)           # 摘要表物件：$fit、$items、$se、$parameters（含 prior）、$covariance、$reliability
coef(fit)              # 後驗平均數（具名向量）
estimates(fit)         # 後驗摘要（資料框：est, sd, q025, q975, rhat, ess；欄名與 birt 相同）
scores(fit)            # 每人的 ability / speed 後驗平均數與 SD（帶 id）
reliability(fit); convergence(fit); fit_indices(fit)
plot(fit)              # 題目參數與 95% 區間；plot(fit, "trace") 軌跡圖

qs <- gibbs(rtirt_cross(dat, quantile = c(.1, .25, .5, .75, .9)), n_iter = 10000, n_chain = 3)
coef(qs); plot(qs)     # rho[item] 對分位數的曲線
```

- **範例：** `system.file("examples", package = "birtRcpp")` 有 13 個完整腳本：ex01–ex06 是內建模型，以 `ecm()` 為主（快速開始、五個模型、分位數（Gibbs）、`ex04_ecm.R`：ECM、`prior = "map"` 與 `init = "ecm"`、模型比較、摘要表與 `equations()`），ex07–ex13 是積木與延伸模型（testlet、mirt 加 RT、條件相依、LNIRT、變動速度、預知題目）。
- **摘要表：** `summary()` 一次算好所有表格，回傳 S3 物件（`rtirt_summary`），版面參考 tppcm 的 `irt_pars()` 與 blavaan。
  - `s$items`、`s$se`：每題一列的寬表（a、b、λ、σ²、ρ 與其 SD 或 SE）。
  - `s$parameters`：每個參數一列，含 Prior 欄；`as.data.frame(s)` 也會得到它。
  - `print(s, tables = "parameters")`：只印部分表格。
- **模型式子：** `equations(model)` 輸出測量模型、人參數的分配與先驗。抽樣前也能用；預設輸出 Unicode，`format = "latex"` 輸出 LaTeX。
- **資料：** `resp`、`time`、`cov` 可以是矩陣／資料框，或配合 `data =` 給欄名。`time` 是秒數（套件取 log）；已取 log 的時間用 `log_time =`。
- **防呆：**
  - 時間看起來已經取過 log（全部 < 15），或 `log_time` 看起來是秒數時會警告；
  - 0/1 以外的編碼（如 1/2）、遺漏值會報錯並說明怎麼處理；
  - 共變數沒有置中時會警告（潛在回歸沒有截距）；
  - 模型用不到 `cov` 時（`rtirt_cross`、`rtirt_null`）會警告。
- **維度：** 由資料推得，`set_cond()` 只在模擬研究需要。
- **MCMC 設定：** 是 `gibbs()` 的參數（`n_iter`、`n_chain`、`n_burnin`、`n_thin`、`seed`）。
- **分位數：** 是模型的參數：`rtirt_cross(dat, quantile = .25)`、`rtirt_latent(dat, quantile = .25)`；給向量就得到一組 `rtirt_qset`。
- **`speed_var`：** `"free"`（預設）估計速度變異數；`"fixed"` 固定為 1（Julia 的 `cov2one = true`）。log 秒數的速度變異數遠小於 1（TIMSS 約 0.09–0.12），第 5 章固定 Var(speed) = 1 會過度限制：以 RTMB 驗證，−2logL 差約 870–930。
- **參數名稱：**

  | 模型 | 名稱與意義 |
  |---|---|
  | `rtirt_latreg` | `cor_ability_speed`：ability 與 speed 的相關 |
  | `rtirt_latent` | `b_ability`：speed 對 ability 的迴歸係數；`beta[x]`：共變數效果 |
  | `rtirt_cross` | `rho[題名]`：θ 與各題 log RT 的交叉關係 |
  | 所有 RT 模型 | `var_speed`：速度變異數 |

| Julia | birtRcpp | Julia 名別名 |
|---|---|---|
| `setCond(nSubj=…, qRt=…)` | `set_cond(...)`（只在模擬需要） | `setCond()` |
| `InputData(Y=, T=, X=)` | `input_data(resp = , time = , cov = , id = , data = )` | `InputData()`、`InputData4R()` |
| `GibbsMlIrt`、`GibbsRtIrtNull`、`GibbsRtIrt` | `mlirt(dat)`、`rtirt_null(dat)`、`rtirt_latreg(dat)` | 同 Julia 名 |
| `GibbsRtIrtLatent`、`…Qr` | `rtirt_latent(dat, quantile = NULL / q)` | 同 Julia 名（`Qr` 用 `Cond$q_rt`） |
| `GibbsRtIrtCross`、`…Qr` | `rtirt_cross(dat, quantile = NULL / q)` | 同 Julia 名 |
| `sample!(MCMC)` | `fit <- gibbs(model, n_iter = , n_chain = , seed = )` | `sampleGibbs()` |
| `coef(MCMC)`、`precis(MCMC)` | `summary(fit)`、`precis(fit)` / `estimates(fit)` | — |
| `getDic`、`checkConvergence` | `dic()` / `fit_indices()`、`convergence()` | 同 Julia 名 |
| `comparePara`、`runSimulation` | `compare_para`、`run_simulation` | 同 Julia 名 |
| `getBias`、`getRmse`、`getCorr` | `get_bias`、`get_rmse`、`get_corr` | 同 Julia 名 |
| `setTruePara<模型>`、`setData<模型>` | `sim_para(cond, model)`、`sim_data(cond, para, model)` | 同 Julia 名 |
| `MCMC.Post.mean.β` | `fit$post$mean$beta` | — |

**S3 類別：**
- 模型物件的類別是 `c(<模型>, "rtirt")`，例如 `c("rtirt_cross", "rtirt")`；均值與分位數版本共用同一個類別，以 `fit$qr` 和 `fit$q` 區分。
- 多個分位數的類別是 `rtirt_qset`。
- 其他物件的類別：`sim_cond`、`input_data`、`input_para`。

**與 birt 共用的存取函數**（`estimates`、`scores`、`reliability`、`fit_indices`、`convergence`）：兩個套件同時載入時，各自的方法會互相登錄到對方的泛型函數，所以不論誰後載入都能正確派發。

## 部分積分的抽樣器：`gibbs(collapse = TRUE)`

`collapse = TRUE` 把潛在變數積掉再抽（van Dyk & Park, 2008 的 partially collapsed Gibbs）：

- **speed：** 抽 ability、RT 題目參數 (λ, ρ) 與第 4 章的結構迴歸時積掉 speed，再從完整條件分配重抽 speed；`mlirt()` 的 β 則積掉 ability。
- **Pólya–Gamma 變數：** 2PL 的 (a, b) 積掉 PG 變數後，在截距參數化 (a, d = −ab) 下一起抽（Fisher scoring 提議的 Metropolis–Hastings），之後照常做 PG 步驟。

```r
fit <- gibbs(rtirt_cross(dat, quantile = 0.1), n_iter = 12000, collapse = TRUE)
```

- 後驗與原本的抽樣器相同（模擬與 TIMSS 上，所有參數的後驗平均數差異都在蒙地卡羅誤差內）。
- TIMSS 上各模型的最小 ESS 增為 1.9–7 倍，原本最慢的是困難題的 a、b。換算成每秒是 1.3–5 倍，因為每次迭代多花約 35–50% 時間。
- 沒有變快的情況：1PL、`rtirt_latreg(speed_var = "fixed")`。預設仍是原本的抽樣器。

## 計算步驟：`algorithm()`

```r
algorithm(fit)                                  # Gibbs：每一步的完整條件分配（先驗的數值已代入）
algorithm(model, collapse = TRUE)               # 部分積分的版本
algorithm(em_fit)                               # EM：E-step、M-step、SQUAREM、標準誤
algorithm(fit, format = "latex")                # LaTeX aligned 區塊
```

每一步右邊註明執行它的積木函數（`draw_*()`），可以照著在 R 裡組裝新的抽樣器（積木教學有完整範例）。

## 檢查抽樣器：`check_sampler()`

```r
check_sampler(model)                       # Geweke (2004) 的聯合分配檢定
check_sampler(model, method = "sbc")       # 模擬式校準（Talts et al., 2018）
check_sampler(model, collapse = TRUE)      # 部分積分的抽樣器
```

兩種方法都只用先驗與概似，檢查的是真正的 C++ 抽樣器：
- Geweke：「抽樣一輪 → 依目前參數重新模擬資料」交替進行，參數鏈的分配必須等於先驗；
- SBC：從先驗抽真值、模擬資料、跑抽樣器，真值在後驗抽樣中的排名必須均勻。

先驗必須是固定、不依賴資料的（`rtirt_priors_julia()` 不行）。鏈要夠長，ESS 太小時會提醒。

## 反應時間的檢定

```r
pf <- person_fit(gibbs_fit)          # 作答時間的 person fit l^t（Marianti et al., 2014）
ci_test(em_fit)                      # 作答與時間的條件獨立 LM 檢定（van der Linden & Glas, 2010）
preknowledge_test(em_fit, compromised = items)   # 洩題偵測（Sinharay, 2020）
convergence(fit, diagnostics = "rank")   # rank-normalized R̂、bulk / tail ESS（Vehtari et al., 2021）
```

`inst/examples/ex11_lnirt.R` 用積木組出 LNIRT 的模式（probit、時間鑑別度、題目層多元常態），在公開的認證考試資料上重現 Fox & Marianti (2017) 的結果：person fit 不適配人數 124 人與原文相同，每秒有效樣本數約為 LNIRT 的 17 倍。`ex12_variable_speed.R` 重現 Fox & Marianti (2016) 在西洋棋資料上的可變速度模式（ρ(θ, ζ₀) = .72），`ex13_preknowledge.R` 重現 Sinharay (2020) 的 Table 2。

## 先驗：`rtirt_priors()`

```r
gibbs(model, priors = rtirt_priors(a = c(mean = 1, sd = 0.5), rho = c(sd = 0.3)))
ecm(model, prior = "map", priors = rtirt_priors(a = c(1, 0.5)))      # 題目參數的先驗
```

預設是不依賴資料的一般常用先驗，所有完整條件分配仍是共軛的（純 Gibbs）：

| 參數 | 預設 | Gibbs 步驟 |
|---|---|---|
| a | N⁺(1, 1) | 截斷常態 |
| b | N(0, 3²)，不截斷 | 常態 |
| λ | N(0, 10²)，不截斷 | 常態 |
| RT 殘差 SD、speed SD | half-t(3, 1) | 兩步 IG（Huang & Wand, 2013 的輔助變數） |
| 迴歸係數、ρ、c | N(0, 1) | 常態 |

half-t 先驗寫成 σ² | u ~ IG(ν/2, ν/u)、u ~ IG(1/2, 1/A²)；每輪先抽 u | σ² ~ IG((ν+1)/2, ν/σ² + 1/A²)，再抽 σ² | u, 資料 ~ IG(ν/2 + n/2, ν/u + SS/2)。

`rtirt_priors_julia()` 是 Julia 版的先驗（b 截在 [−4, 4]、λ 以 log 時間的平均數與標準差為中心、變異數 IG(0.001, 0.001)），和舊版的抽樣逐位元相同，用來重現博論結果。變異數也可以直接指定 IG：`sigma2t = c(shape = 1, scale = 1)`。摘要表的 Prior 欄與 `equations()` 顯示實際用到的先驗。

## Score-based 的參數不變性檢定：`score_test()`

```r
fit <- ecm(rtirt_cross(dat), se = FALSE)
score_test(fit, raw$gender, by_item = TRUE)       # 性別 DIF，逐題並做 Holm 校正
score_test(fit, rowSums(dat$log_t))               # 沿總作答時間（聯合模型中對 a、b 有效）
```

- 只需要一次 `ecm()` 配適。每人的 score 由 `estfun()` 給出，也可以直接交給 strucchange 或 partykit（`mob()`）。
- 調節變項和能力有關時（impact），請把它放進潛在迴歸（`rtirt_latreg(input_data(..., cov = z))`）再檢定；函數會提醒。
- DIF 樹在 birt：`birt::dif_tree()`（partykit 的 `mob()`，節點用 RTMB 最大概似）。

## 積木：宣告模型，用 Gibbs、ECM、CAVI 估計（`?block_model`）

分工：birt 是整套自動估計（JAGS、RTMB、Laplace、ELGM）；birtRcpp 的積木是「逐塊用完整條件分配」的估計法，模型宣告一次，三種估計法共用同一份步驟清單：

```r
m <- block_model(Y = Y,
  theta[person] ~ normal(0, 1),
  a[item] ~ normal(1, 1, lower = 0),
  b[item] ~ normal(0, 3),
  Y ~ bernoulli_logit(a * (theta - b)))
steps <- list(step_mala(a, b), step_gibbs(theta, a, b), step_shift(theta = 1, b = 1))
fit_g <- gibbs(m, steps)          # Gibbs：後驗樣本
fit_e <- ecm(m, steps)              # ECM：後驗眾數（prior = "ml"：最大概似）與標準誤
fit_v <- vi(m, list(step_gibbs(theta, a, b)))   # CAVI：每塊一個近似分配 q
algorithm(fit_e)                          # 每一步的式子
check_sampler(m, steps = steps); check_ecm(m, steps); check_vi(m)   # 三種檢查（Geweke 用小資料較快）
```

- **同一個條件分配，三種用法：** Gibbs 抽一個值；ECM 移到眾數（θ 用自適應 Gauss–Hermite 積掉，ω 換成 E[ω]）；CAVI 設 q 為同型分配、係數換成期望值（Durante & Rigon, 2019）。
- **系統處理的事：** Pólya–Gamma 變數（需要時自動重抽，「ω 過期」不會發生）、log 尺度的 Jacobian、截斷、維度（`[person]`、`[item]`）、先驗只寫一次。
- **cut：** `step_gibbs(theta, cut = "logT")`，θ 只由作答決定（Plummer, 2015）。
- **LM（分數）檢定：** `score_test(fit_e, add = Y ~ delta[item] * z)` 只用虛無模型的 ECM 估計，檢定加在 logit（或常態的平均數）上的一項，例如沿組別 z 的 DIF（Glas, 1998；Glas & van der Linden, 2010）；分數用 Fisher identity，資訊用 Louis 公式；整體與逐題（Holm）。
- **只出現在常態部分、而且線性的人層變項**（例如速度 ζ），ECM 與 CAVI 在給定 θ 時用公式精確積分。
- **計算：** 式子編成多線性形式，由 C++（RcppArmadillo）計算。
- **驗證：** 宣告的交叉關係模型和 `gibbs(rtirt_cross())` 的後驗相同（41 個參數，最大差 2.0 個 Monte Carlo 標準誤），和 `ecm(rtirt_cross())` 的估計相同（差 1e-4 以內，標準誤相同）；58 組隨機堆疊的步驟全部正確；CAVI 的平均數接近後驗，但 SD 偏小（平均場的已知性質）。
- **手寫每一步**（`run_sampler()`、`draw_*()`、`run_em()`、`estep_irt()`/`mstep_items()`）仍保留，用於宣告寫法還沒有的結構；`inst/examples/ex07`–`ex13` 是手寫的文獻範例。
- 積木教學逐行解釋：第 3 節三種估計法、第 4 節三種檢查、第 6 節作答時間與 cut、第 9 節手寫。

## 最大概似：`ecm()`

常態模型（`mlirt`、`rtirt_null`、`rtirt_latreg`、`rtirt_latent`、`rtirt_cross`，不含 `quantile`）可以用 Bock–Aitkin EM 估計邊際最大概似：

```r
fit <- ecm(rtirt_cross(dat))      # 估計值、標準誤、logLik、AIC、BIC
fit; summary(fit); estimates(fit); scores(fit); fit_indices(fit); logLik(fit); vcov(fit)
g <- gibbs(fit, init = "ecm")     # 從 EM 解出發抽樣，縮短 burn-in
estimates(g, method = "em")      # 抽樣後仍可取出 EM 的表
```

**何時用 EM、何時用 Gibbs：**

| 需求 | `ecm()` | `gibbs()` |
|---|---|---|
| 估計值與 SE | 1–2 秒 | 數十秒 |
| AIC、BIC、概似比檢定 | 有 | 沒有，改用邊際 DIC、WAIC、LOO |
| 小樣本、參數接近邊界 | 可能發散，可用 `prior = "map"` | 先驗正則化，較穩定 |
| 分位數（ALD）模型 | 不支援 | 支援 |
| 衍生量的完整不確定性 | Wald 近似 | 後驗抽樣 |

**診斷：**
- `fit$ecm$trace`：每次迭代的 log-likelihood，應該單調上升到平台。
- `fit$ecm$converged`：是否收斂。
- 用 `ecm(model, nodes = 61)` 重估一次，log-likelihood 應該幾乎不變，以此檢查積分精度。
- 有巢套關係的模型（例如 `rtirt_null` 與 `rtirt_latreg`）可以用 `2 * (logLik(m1) - logLik(m0))` 做概似比檢定。

例子與小樣本設定見入門教學第 8 節；同樣的 EM 也可以用宣告的模型寫（積木教學第 3 節）。

**作法：**
- **E 步：**
  - θ 放在 `marginal_loglik()` 用的 Gauss–Hermite 節點上（N(μθ, 1)，預設 41 點）；
  - 給定 θ 時，ζ 的後驗是常態，以封閉式積掉，並得到 E[ζ]、E[ζ²]。
- **M 步：**
  - 題目的 logistic 部分用 Pólya–Gamma 下界：E[ω] = tanh(η/2)/(2η)（Polson, Scott & Windle, 2013），每題 (a, d) 是一個 2×2 加權最小平方；這是 EM 裡的 MM 步，log-likelihood 不會下降；
  - RT 參數（λ、σ²、ρ）與結構參數（β、相關、變異數、b_ability）都有封閉解；`speed_var = "fixed"` 時相關用一維搜尋。
- **加速：** SQUAREM（Varadhan & Roland, 2008）；外插會降低概似時改回一般 EM 步。
- **標準誤：** 觀察訊息矩陣；以 Fisher identity 算出解析 score，再對它做數值微分，最後用 delta method 換回原尺度。Wald 區間在 a、σ²、var_speed 取 log 尺度，相關取 Fisher-z 尺度。
- **識別：** 與 `gibbs()` 相同（θ 變異數 1、迴歸沒有截距）；`intercept = TRUE` 不可用，因為沒有先驗時截距不可識別。
- **`prior = "map"`：** 題目參數加上 `gibbs()` 的先驗（取後驗眾數），結構參數仍是最大概似；適合小樣本。
- 分位數（ALD）模型請用 `gibbs()`。
- **共同節點的摺疊：** 沒有潛在回歸時每人共用節點，E 步與 M 步化為 Bock–Aitkin 人工資料上的矩陣乘法（迭代完全相同；5000 人 × 60 題：配適 1 秒、含標準誤 20 秒）。
- **組件（`?em_components`）：** `ecm()` 由 `calc_loglik_nodes()`、`calc_node_weights()`、`calc_item_sums()`、`update_items_logistic()`、`em_squarem()`、`calc_hessian()` 等組成，可以用來組自己的 EM。

**驗證**（專案資料夾的 `scripts/40_em_check.R` → `results/em_check.txt`，不隨套件發布；N = 1000、10 題）：

| 模型 | logLik 與 birt RTMB（AGHQ）差距 | 參數最大差距 | SE 最大相對差距 | 與 Gibbs 後驗平均差距（後驗 SD 倍數） | EM 秒數（含 SE） | RTMB 秒數 |
|---|---|---|---|---|---|---|
| mlirt | 4e-7 | 2e-6 | 1e-6 | 0.09 | 0.8 | 2.9 |
| rtirt_null | 7e-6 | 2e-5 | 2e-5 | 0.19 | 1.3 | 18.9 |
| rtirt_latreg | 2e-5 | 1e-5 | 2e-5 | 0.13 | 1.4 | 22.2 |
| rtirt_latent | 3e-5 | 6e-6 | 9e-6 | 0.14 | 1.3 | 24.7 |
| rtirt_cross | 2e-6 | 3e-5 | 2e-5 | 0.14 | 2.0 | 25.2 |

- 一般 EM（不加速）每一步的 log-likelihood 都上升。
- EM 的 logLik 與 `marginal_loglik()` 在同一組估計值下的計算相差 < 2e-6。

## 與 Julia 原版的差異（修正）

| 編號 | 問題（Julia） | birtRcpp |
|---|---|---|
| B1 | 第 4、5 章 θ 的條件分配只用作答正確性 | 加入 RT 端（第 5 章）／結構回歸端（第 4 章） |
| B2 | 第 3 章 θ、ζ 各自用邊際事前；IW 後強制改成相關矩陣 | θ \| ζ、ζ \| θ 用同一個共變異；Var(θ)=1 精確成立。`speed_var = "free"` 時 (c, v) 共軛；`"fixed"` 時 ρ 用 Metropolis |
| B3 | 第 4 章 QR 的 β 用未加權 OLS 點值 | 以 ALD 權重 1/(k₂sν) 聯合抽樣 |
| B4 | 多條「鏈」共用同一個狀態 | 真正獨立的鏈，各自隨機起始值 |
| B5 | DIC 含 burn-in | 只用 burn-in 之後（仍是條件化於人參數與 ν 的完整資料 DIC） |
| B6 | β 的事前精確度 `1/σ² .+ 矩陣` 加到每一格 | I/σ² |
| B7 | 第 4 章 QR 的尺度更新用矩陣除法 `/` | 逐元素 |
| — | b 被 clamp 到 [−4, 4] | 截尾常態 |
| — | 無截距時先聯合抽樣再把截距設 0 | 設計矩陣本身不含截距 |
| — | 模擬資料：logT 截在 0 以上；交叉關係資料沒用到真值 σ²t | 與配適的模型完全一致 |
| 新增 | λ 與 ζ 的共同位移（似然只看 λ − ζ）混合很慢 | 每次迭代多一步精確的位置移動（平移群上的 Gibbs，Liu & Sabatti, 2000），不改變目標分配 |
| 新增 | 第 5 章：(ρ_i − d, ζ_j + dθ_j) 不改變似然，d 只由先驗識別，混合很慢 | 精確的剪切移動（Jacobian 1）。極端分位數的 R̂ 從 1.42 降到 1.01 |

## 驗證

1. **抽樣器單元測試：**
   - Pólya-Gamma：核對解析平均數與變異數；
   - 截尾常態：含遠尾；
   - 逆高斯。
2. **模擬回收：** 每個模型都有，見 `tests/testthat`。R CMD check：Status OK。
3. **博論資料（N = 631，14 題）三方比對：**
   - 比對對象：修正後的 Julia（專案資料夾的 `julia_check/RtIrtFixed.jl`，不隨套件發布）、birtRcpp、似然與先驗完全相同的 JAGS 模型。
   - 範圍：11 個模型（第 3 章 ×2、第 4 章 MR 與 p = .1/.5/.9、第 5 章 MR 與 p = .1/.5/.9、第 5 章 speed 變異數自由），每個自由參數都比較。
   - z 的定義：兩方法後驗平均數之差 ÷ 合併蒙地卡羅誤差。
   - **birtRcpp vs Julia：** 所有參數 |z| ≤ 3.4，差距 ≤ 0.37 個後驗 SD，SD 比 0.87–1.23。
   - **birtRcpp vs JAGS、Julia vs JAGS：** 除了第 5 章 p = .1，|z| ≤ 3.9。
   - **第 5 章 p = .1 的例外：** JAGS 的所有 ρ_i 一起偏 −0.02。原因是 JAGS 在只由先驗識別的「共同 ρ 位移」方向混合不良（JAGS 自己的 R̂ 1.06–1.18）。birtRcpp 有精確的群移動步驟，R̂ 1.00。
   - 細節：專案資料夾的 `julia_check/verify_samplers.csv`、`julia_check/verify_samplers.txt`（不隨套件發布）。
4. **速度：** 12000 次迭代 × 3 條鏈，約 2–3 分鐘；JAGS 需 8–40 分鐘。

## 安裝

```r
install.packages("birtRcpp")                       # CRAN 上架之後
# install.packages("remotes")
remotes::install_github("jiewenTsai/birtRcpp")     # GitHub 上的版本（需要 C++ 編譯器）
```

教學與範例：[birtExamples](https://github.com/jiewenTsai/birtExamples)。
