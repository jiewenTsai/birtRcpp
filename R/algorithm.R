# algorithm(): the computation steps of a fit (Gibbs full conditionals or EM steps) ----------------
#
# The steps are written once in LaTeX and converted to Unicode by tex_to_unicode(); the prior
# constants of the fit are substituted into the formulas. Persons j = 1..N, items i = 1..K, as in
# equations(). Each step names the C++ block that performs it and its R wrapper (?blocks).

#' Computation steps of a fit
#'
#' Lists the steps of the algorithm behind a fit: for [gibbs()] every full conditional
#' distribution in the order of one iteration (with `collapse = TRUE` the partially collapsed
#' version), for [ecm()] the E- and M-steps. The prior constants of the fit are substituted into
#' the formulas, and each step names the building block that draws it ([blocks]), so a step can
#' be run, replaced or checked from R. Works for unfitted models too.
#' @param object A model object (see [models]), sampled or not, or an `rtirt_qset`.
#' @param method `"gibbs"` or `"ecm"`; by default the method of the fit (`"gibbs"` for an unfitted
#'   model).
#' @param collapse For `"gibbs"`: the partially collapsed sampler; by default the setting of the
#'   fit (`FALSE` for an unfitted model).
#' @param format `"unicode"` or `"latex"` (an `aligned` block).
#' @param ... Unused.
#' @return An `rtirt_algorithm` object (printed); `as.character()` gives the lines.
#' @seealso [equations()] for the model, [blocks] for the steps as R functions.
#' @examples
#' cond <- set_cond(n_subj = 100, n_item = 5)
#' m <- rtirt_cross(sim_data(cond, sim_para(cond, "cross"), "cross"))
#' algorithm(m)                       # Gibbs full conditionals
#' algorithm(m, collapse = TRUE)      # partially collapsed Gibbs
#' algorithm(m, method = "ecm")       # ECM steps
#' algorithm(m, format = "latex")
#' @export
algorithm <- function(object, ...) UseMethod("algorithm")

#' @rdname algorithm
#' @export
algorithm.rtirt_qset <- function(object, method = NULL, collapse = NULL, format = c("unicode", "latex"), ...) {
  out <- algorithm(object$models[[1]], method = method, collapse = collapse, format = format, ...)
  out$sections <- c(out$sections, list(list(title = "Quantiles", lines = list(
    al_line("q", "\\in", sprintf("\\{%s\\}", paste(format(object$quantile), collapse = ", ")), "one model and sampler per quantile", out$format)))))
  out
}

#' @rdname algorithm
#' @export
algorithm.rtirt <- function(object, method = NULL, collapse = NULL, format = c("unicode", "latex"), ...) {
  format <- match.arg(format)
  fitted_em <- is.null(object$post) && !is.null(object$ecm)
  method <- match.arg(method %||% if (fitted_em) "ecm" else "gibbs", c("gibbs", "ecm"))
  if (method == "ecm" && isTRUE(object$qr)) stop("ecm() fits the normal models; quantile models use gibbs()", call. = FALSE)
  collapse <- collapse %||% isTRUE(object$post$settings$collapse)
  ctx <- al_context(object)
  secs <- if (method == "gibbs") al_gibbs(ctx, collapse) else al_em(ctx, object)
  secs <- lapply(Filter(function(s) length(s$lines) > 0, secs), function(s) {
    s$lines <- lapply(s$lines, function(l) al_render(l, format))
    if (format == "latex") s$title <- al_tex_text(s$title)
    s })
  structure(list(sections = secs, format = format), class = c("rtirt_algorithm", "rtirt_equations"))
}

#' @rdname algorithm
#' @param x An `rtirt_algorithm` object.
#' @export
print.rtirt_algorithm <- function(x, ...) { cat(as.character(x), sep = "\n"); invisible(x) }

# ---- context: model type, options and prior constants ---------------------------------------------
al_context <- function(object) {
  info <- .model_info[[class(object)[1]]]; S <- object$settings; PS <- object$post$settings
  itemtype <- PS$itemtype %||% S$itemtype; speed_var <- PS$speed_var %||% S$speed_var
  X <- PS$X %||% (if (info$X) cbind(if (isTRUE(S$intercept)) 1, object$data$X) else NULL)
  pr <- fit_priors(object)
  q <- if (isTRUE(object$qr)) object$q else NA_real_
  list(eng = info$engine, rt = info$rt, qr = isTRUE(object$qr), q = q, one_pl = identical(itemtype, "1pl"),
       fixed = identical(speed_var, "fixed") && info$engine %in% c("rtirt", "cross"), P = NCOL(X) * !is.null(X),
       pr = pr, lp = if (info$rt) lambda_prior(pr, object$data$log_t) else NULL, cond = object$cond,
       fit = object$post$settings, N = nrow(object$data$Y), K = ncol(object$data$Y))
}

# ---- formatting helpers ------------------------------------------------------------------------
# precision 1 / sd^2 and the term mean / sd^2 of a normal prior (NULL when the mean is 0)
al_prec <- function(sd, n = NULL) { v <- fmt_num(sd^2); if (is.null(n)) { if (sd == 1) "1" else sprintf("1/%s", v) } else if (sd == 1) n else sprintf("%s/%s", n, v) }
al_mterm <- function(mean, sd, n = NULL) {
  if (mean == 0) return(NULL)
  m <- if (is.null(n)) fmt_num(mean) else sprintf("%s%s", n, fmt_num(mean))
  if (sd == 1) m else sprintf("%s/%s", m, fmt_num(sd^2))
}
al_sum <- function(...) { t <- unlist(list(...)); t <- t[nzchar(t)]; if (!length(t)) return("0"); gsub("\\+ -", "- ", paste(t, collapse = " + ")) }

# a line: LaTeX source, rendered at the end
al_line <- function(lhs, rel = "", rhs = "", note = "", format = NULL) {
  l <- list(lhs = lhs, rel = rel, rhs = rhs, note = note)
  if (is.null(format)) l else al_render(l, format)
}
L <- al_line
al_render <- function(l, format) {
  if (format == "latex") return(list(lhs = l$lhs, rel = l$rel, rhs = l$rhs, note = if (nzchar(l$note)) sprintf("\\text{%s}", al_tex_text(l$note)) else ""))
  list(lhs = tex_to_unicode(l$lhs), rel = tex_to_unicode(l$rel), rhs = tex_to_unicode(l$rhs), note = l$note)
}
# plain text inside \\text{}: escape the special characters
al_tex_text <- function(x) gsub("~", "\\\\textasciitilde{}", gsub("\\^", "\\\\textasciicircum{}", tex_esc(x)))
al_sec <- function(title, ...) list(title = title, lines = Filter(Negate(is.null), list(...)))

# the LaTeX subset used here, converted to Unicode
tex_to_unicode <- function(x) {
  if (!nzchar(x)) return(x)
  gr <- c(alpha = "\u03b1", beta = "\u03b2", gamma = "\u03b3", delta = "\u03b4", epsilon = "\u03b5", zeta = "\u03b6",
          eta = "\u03b7", theta = "\u03b8", kappa = "\u03ba", lambda = "\u03bb", mu = "\u03bc", nu = "\u03bd", xi = "\u03be",
          pi = "\u03c0", rho = "\u03c1", sigma = "\u03c3", tau = "\u03c4", phi = "\u03c6", chi = "\u03c7", psi = "\u03c8",
          omega = "\u03c9", Sigma = "\u03a3", Delta = "\u0394")
  sym <- c(infty = "\u221e", sum = "\u03a3", prod = "\u220f", mid = "|", sim = "\u223c", cdot = "\u00b7", le = "\u2264",
           ge = "\u2265", leftarrow = "\u2190", times = "\u00d7", ldots = "\u2026", nabla = "\u2207", propto = "\u221d",
           approx = "\u2248", top = "\u1d40", `in` = "\u2208", arg = "arg", ell = "\u2113", partial = "\u2202", quad = "   ", min = "min", max = "max", log = "log", exp = "exp", tanh = "tanh",
           atanh = "atanh", left = "", right = "")
  sub <- c(`0` = "\u2080", `1` = "\u2081", `2` = "\u2082", `3` = "\u2083", `4` = "\u2084", `5` = "\u2085", `6` = "\u2086",
           `7` = "\u2087", `8` = "\u2088", `9` = "\u2089", `+` = "\u208a", `-` = "\u208b", `(` = "\u208d", `)` = "\u208e",
           a = "\u2090", e = "\u2091", o = "\u2092", x = "\u2093", h = "\u2095", k = "\u2096", l = "\u2097", m = "\u2098",
           n = "\u2099", p = "\u209a", s = "\u209b", t = "\u209c", i = "\u1d62", j = "\u2c7c", r = "\u1d63", u = "\u1d64",
           v = "\u1d65", `,` = ",", "\u03b2" = "\u1d66", "\u03b3" = "\u1d67", "\u03c1" = "\u1d68", "\u03c6" = "\u1d69", "\u03c7" = "\u1d6a")
  sup <- c(`0` = "\u2070", `1` = "\u00b9", `2` = "\u00b2", `3` = "\u00b3", `4` = "\u2074", `5` = "\u2075", `6` = "\u2076",
           `7` = "\u2077", `8` = "\u2078", `9` = "\u2079", `+` = "\u207a", `-` = "\u207b", `(` = "\u207d", `)` = "\u207e",
           n = "\u207f", i = "\u2071", t = "\u1d57", r = "\u02b3", "\u1d40" = "\u1d40", `*` = "*")
  rep_fun <- function(x, pattern, f) {
    repeat {
      m <- regexpr(pattern, x, perl = TRUE)
      if (m < 0) return(x)
      st <- attr(m, "capture.start"); ln <- attr(m, "capture.length")
      caps <- substring(x, st, st + ln - 1)
      x <- paste0(substr(x, 1, m - 1), f(caps), substr(x, m + attr(m, "match.length"), nchar(x)))
    }
  }
  x <- gsub("\\\\\\{", "\002", gsub("\\\\\\}", "\003", x))
  x <- rep_fun(x, "\\\\(?:text|mathrm|operatorname)\\{([^{}]*)\\}", function(c) gsub("-", "\001", c))
  x <- gsub("\\\\mathcal\\{N\\}", "N", x); x <- gsub("\\\\mathbf\\{1\\}", "\U0001d7cf", x)
  for (n in names(gr)) x <- gsub(sprintf("\\\\%s(?![A-Za-z])", n), gr[[n]], x, perl = TRUE)
  for (n in names(sym)) x <- gsub(sprintf("\\\\%s(?![A-Za-z])", n), sym[[n]], x, perl = TRUE)
  x <- rep_fun(x, "\\\\tilde\\s*\\{?([^{}\\s])\\}?", function(c) paste0(c, "\u0303"))
  x <- rep_fun(x, "\\\\hat\\s*\\{?([^{}\\s])\\}?", function(c) paste0(c, "\u0302"))
  x <- rep_fun(x, "\\\\bar\\s*\\{?([^{}\\s])\\}?", function(c) paste0(c, "\u0304"))
  conv <- function(c, map, mark) {
    ch <- strsplit(c, "")[[1]]
    if (length(ch) && all(ch %in% names(map))) paste(map[ch], collapse = "")
    else if (!grepl("[ +,/-]", c)) paste0(mark, c) else sprintf("%s(%s)", mark, c)    # mark: placeholder, restored below
  }
  x <- rep_fun(x, "_\\{([^{}]*)\\}", function(c) conv(c, sub, "\004"))
  x <- rep_fun(x, "\\^\\{([^{}]*)\\}", function(c) conv(c, sup, "\005"))
  x <- rep_fun(x, "_([^{}\\s\\\\(),|/])", function(c) conv(c, sub, "\004"))
  x <- rep_fun(x, "\\^([^{}\\s\\\\(),|/])", function(c) conv(c, sup, "\005"))
  paren <- function(s) if (grepl("[-+ /]", s) && !grepl("^\\(.*\\)$", s)) sprintf("(%s)", s) else s
  x <- gsub("\\\\\\|", "\u2016", x)
  x <- rep_fun(x, "\\\\t?frac\\{([^{}]*)\\}\\{([^{}]*)\\}", function(c) paste0(paren(c[1]), "/", paren(c[2])))
  x <- rep_fun(x, "\\\\sqrt\\{([^{}]*)\\}", function(c) paste0("\u221a", if (nchar(c) > 1) sprintf("(%s)", c) else c))
  x <- gsub("\\\\[,;!]", " ", x); x <- gsub("\\\\ ", " ", x)
  x <- gsub("[{}]", "", x)
  x <- gsub("(?<![A-Za-z])-|-(?![A-Za-z])", "\u2212", x, perl = TRUE)
  x <- gsub("\001", "-", x); x <- gsub("\002", "{", x); x <- gsub("\003", "}", x)
  x <- gsub("\004", "_", x); x <- gsub("\005", "^", x)
  gsub(" {2,}", " ", gsub("^ +| +$", "", x))
}

# a variance (or ALD scale) step: likelihood var^-alpha exp(-beta / var) and the prior of the fit
al_var <- function(v, aux, alpha, beta, val, fam, note) {
  if (fam == "inverse-gamma")
    return(list(L(sprintf("%s \\mid \\cdot", v), "\\sim", sprintf("\\mathrm{IG}(%s + %s,\\; %s + %s)", fmt_num(val[["shape"]]), alpha, fmt_num(val[["scale"]]), beta),
                  sprintf("%s (conjugate inverse gamma)", note))))
  df <- val[["df"]]; A <- val[["scale"]]
  list(L(sprintf("%s \\mid %s", aux, v), "\\sim", sprintf("\\mathrm{IG}(%s,\\; %s/%s + %s)", fmt_num((df + 1) / 2), fmt_num(df), v, al_prec(A)),
         sprintf("%s: half-t(%s, %s) prior through the auxiliary %s (Huang & Wand, 2013)", note, fmt_num(df), fmt_num(A), "variable")),
       L(sprintf("%s \\mid %s, \\cdot", v, aux), "\\sim", sprintf("\\mathrm{IG}(%s + %s,\\; %s/%s + %s)", fmt_num(df / 2), alpha, fmt_num(df), aux, beta), ""))
}

# ---- Gibbs ---------------------------------------------------------------------------------------
al_gibbs <- function(C, collapse) {
  pr <- C$pr; eng <- C$eng; qr <- C$qr; P <- C$P; rt <- C$rt
  ab_joint <- collapse && !C$one_pl
  # RT weights and the ALD shift of the measurement model (cross only) or of the structural model (latent)
  w <- if (eng == "cross") "w_{ij}" else "\\sigma_i^{-2}"
  k1nu <- if (qr && eng == "cross") " + k_1\\nu_{ij}" else ""
  rho_t <- function(s = "+") if (eng == "cross") sprintf(" %s \\rho_i\\theta_j", s) else ""
  xb <- function(b) if (P) sprintf("x_j^\\top\\%s", b) else NULL
  wz <- if (qr) "\\tilde w_j" else "s^{-1}"                     # prior precision of speed (latent)
  shift <- if (qr && eng == "latent") " + k_1\\tilde\\nu_j" else ""
  sb <- pr$beta[["sd"]]

  notation <- al_sec("Notation",
    L("j = 1, \\ldots, N;\\; i = 1, \\ldots, K", "", "", sprintf("persons j and items i (N = %d, K = %d); y_ij: item i, person j", C$N, C$K)),
    L(if (rt) "\\theta_j,\\; \\zeta_j" else "\\theta_j", "", if (rt) "\\text{ability, speed of person } j" else "\\text{ability of person } j", ""),
    L("\\kappa_{ij}", "=", "y_{ij} - 1/2", "Polya-Gamma: logistic likelihood exp(kappa eta - omega eta^2 / 2) given omega"),
    if (eng == "cross") L("w_{ij}", "=", if (qr) "1/(k_2\\sigma_i\\nu_{ij})" else "1/\\sigma_i^2", "precision of log t_ij"),
    if (qr) L("k_1,\\; k_2", "=", sprintf("(1 - 2q)/(q(1 - q)) = %s,\\; 2/(q(1 - q)) = %s", fmt_num((1 - 2 * C$q) / (C$q * (1 - C$q))), fmt_num(2 / (C$q * (1 - C$q)))),
              sprintf("ALD(q = %s) as a normal mixture (Kozumi & Kobayashi, 2011)", fmt_num(C$q))),
    if (qr && eng == "latent") L("\\tilde w_j", "=", "1/(k_2 s\\tilde\\nu_j)", "precision of speed in the quantile regression"),
    if (eng == "rtirt") L("\\mu_{1j},\\; \\mu_{2j}", "=", if (P) "x_j^\\top\\beta_1,\\; x_j^\\top\\beta_2" else "0,\\; 0",
                          "speed_j = mu_2j + c (theta_j - mu_1j) + N(0, v)"),
    L("\\cdot", "", "\\text{all other parameters and the data}", ""))

  # accuracy: (joint MH), Polya-Gamma, b, a
  pa <- pr$a; pb <- pr$b
  acc <- c(
    if (ab_joint) list(
      L("(a_i, d_i),\\; d_i = -a_i b_i", "", "\\text{Metropolis-Hastings, } \\omega \\text{ integrated out}", "draw_items_joint(): first, along the a-b ridge"),
      L("f(a_i, d_i)", "=", sprintf("\\sum_j [y_{ij}\\eta_{ij} - \\log(1 + \\exp \\eta_{ij})] + \\log p(a_i) + \\log p(b_i) - \\log a_i,\\; \\eta_{ij} = a_i\\theta_j + d_i"),
        "log posterior in (a, d); -log a_i is the Jacobian of b = -d/a"),
      L("G_i", "=", sprintf("\\sum_j \\pi_{ij}(1 - \\pi_{ij})(\\theta_j, 1)^\\top(\\theta_j, 1) + \\mathrm{diag}(%s,\\; %s)", al_prec(pa[2]), if (pb[["sd"]] == 1) "1/a_i^2" else sprintf("1/(%s a_i^2)", fmt_num(pb[["sd"]]^2))),
        "Fisher information + prior precision"),
      L("x^*", "\\sim", "\\mathcal{N}(x + G_i^{-1}\\nabla f(x),\\; G_i^{-1}),\\; x = (a_i, d_i)", "Fisher-scoring proposal"),
      L("\\text{accept}", "", "\\min\\{1,\\; \\exp[f(x^*) - f(x)]\\, q(x \\mid x^*)/q(x^* \\mid x)\\}", "then the Polya-Gamma steps below")),
    list(L("\\omega_{ij} \\mid \\cdot", "\\sim", "\\mathrm{PG}(1,\\; a_i(\\theta_j - b_i))", "draw_omega()"),
         L("b_i \\mid \\omega, a, \\theta", "\\sim", paste0("\\mathcal{N}(m_i/p_i,\\; 1/p_i)", fmt_trunc(pb[["lower"]], pb[["upper"]])), "draw_items_pg(): b first"),
         L("p_i", "=", al_sum(al_prec(pb[["sd"]]), "a_i^2 \\sum_j \\omega_{ij}")),
         L("m_i", "=", al_sum(al_mterm(pb[["mean"]], pb[["sd"]]), "a_i \\sum_j (a_i\\omega_{ij}\\theta_j - \\kappa_{ij})"))),
    if (C$one_pl) list(L("a_i", "=", "1", "1PL")) else list(
      L("a_i \\mid \\omega, b, \\theta", "\\sim", "\\mathcal{N}^{+}(m_i/p_i,\\; 1/p_i)", "then a (truncated normal)"),
      L("p_i", "=", al_sum(al_prec(pa[2]), "\\sum_j \\omega_{ij}(\\theta_j - b_i)^2")),
      L("m_i", "=", al_sum(al_mterm(pa[1], pa[2]), "\\sum_j \\kappa_{ij}(\\theta_j - b_i)"))))
  irt <- L("P_j,\\; h_j", "=", "\\sum_i \\omega_{ij}a_i^2,\\; \\sum_i a_i(\\kappa_{ij} + \\omega_{ij}a_ib_i)", "irt_part(): accuracy part of theta_j")

  # theta (not collapsed)
  theta <- switch(eng,
    mlirt = list(L("p_j", "=", "1 + P_j"), L("m_j", "=", al_sum("h_j", xb("beta")))),
    rtirt = list(L("p_j", "=", "1 + P_j + c^2/v"),
                 L("m_j", "=", if (P) "h_j + \\mu_{1j} + (c/v)(\\zeta_j - \\mu_{2j} + c\\mu_{1j})" else "h_j + (c/v)\\zeta_j")),
    latent = list(L("p_j", "=", sprintf("1 + P_j + \\gamma^2 %s", wz)),
                  L("m_j", "=", sprintf("h_j + \\gamma %s(\\zeta_j%s%s)", wz, if (P) " - x_j^\\top\\beta" else "", sub("\\+", "-", shift)))),
    cross = list(L("p_j", "=", "1 + P_j + \\sum_i w_{ij}\\rho_i^2"),
                 L("m_j", "=", sprintf("h_j + \\sum_i w_{ij}\\rho_i(\\lambda_i%s - \\log t_{ij} - \\zeta_j)", k1nu))))
  theta_sec <- c(list(irt, L("\\theta_j \\mid \\cdot", "\\sim", "\\mathcal{N}(m_j/p_j,\\; 1/p_j)", "draw_gauss()")), theta)

  if (!rt) {
    beta_full <- list(L("\\beta \\mid \\theta", "\\sim", "\\mathcal{N}(Q^{-1}h,\\; Q^{-1})", "multivariate normal"),
                      L("Q,\\; h", "=", sprintf("%s + X^\\top X,\\; X^\\top\\theta", al_prec(sb, "I"))))
    beta_coll <- list(L("\\beta \\mid \\omega, a, b", "\\sim", "\\mathcal{N}(Q^{-1}h,\\; Q^{-1})", "theta integrated out (after the accuracy step)"),
                      L("\\hat\\theta_j", "=", "h_j/P_j \\sim \\mathcal{N}(x_j^\\top\\beta,\\; 1 + 1/P_j)", "accuracy part as a pseudo-observation"),
                      L("Q,\\; h", "=", sprintf("%s + \\sum_j \\bar w_j x_j x_j^\\top,\\; \\sum_j \\bar w_j x_j\\hat\\theta_j,\\; \\bar w_j = P_j/(1 + P_j)", al_prec(sb, "I"))))
    steps <- list(if (P && !collapse) c("Latent regression", beta_full), c("Accuracy items", acc),
                  if (P && collapse) c("Latent regression (theta integrated out)", c(list(irt), beta_coll[c(1, 2, 3)])),
                  c("Ability", theta_sec))
    return(c(list(al_settings(C, collapse), notation), al_number(steps)))
  }

  # RT items: lambda, sigma2 (not collapsed)
  lam <- list(L("\\lambda_i \\mid \\cdot", "\\sim", paste0("\\mathcal{N}(m_i/p_i,\\; 1/p_i)", fmt_trunc(C$lp[3], Inf)), "draw_lambda()"),
              L("p_i", "=", al_sum(al_prec(C$lp[2]), if (eng == "cross") "\\sum_j w_{ij}" else "N/\\sigma_i^2")),
              L("m_i", "=", al_sum(al_mterm(C$lp[1], C$lp[2]),
                                   if (eng == "cross") sprintf("\\sum_j w_{ij}(\\log t_{ij} + \\zeta_j%s%s)", rho_t(), if (qr) " - k_1\\nu_{ij}" else "")
                                   else "\\sum_j (\\log t_{ij} + \\zeta_j)/\\sigma_i^2")))
  r_ij <- sprintf("\\log t_{ij} - \\lambda_i + \\zeta_j%s", rho_t())
  sig <- if (qr && eng == "cross")
    c(list(L("r_{ij}", "=", r_ij)), al_var("\\sigma_i", "u_i", "3N/2", "\\sum_j (r_{ij} - k_1\\nu_{ij})^2/(2k_2\\nu_{ij}) + \\sum_j \\nu_{ij}", pr$sigma2t, pr$sigma2t_family, "draw_sigma2(): ALD scale"))
    else c(list(L("S_i", "=", sprintf("\\sum_j (%s)^2", r_ij))), al_var("\\sigma_i^2", "u_i", "N/2", "S_i/2", pr$sigma2t, pr$sigma2t_family, "draw_sigma2()"))
  # speed (not collapsed): prior mean mz and precision of the structural model
  zeta <- switch(eng,
    rtirt = list(L("p_j", "=", "1/v + \\sum_i \\sigma_i^{-2}"),
                 L("m_j", "=", sprintf("%s/v + \\sum_i (\\lambda_i - \\log t_{ij})/\\sigma_i^2", if (P) "(\\mu_{2j} + c(\\theta_j - \\mu_{1j}))" else "c\\theta_j"))),
    latent = list(L("p_j", "=", sprintf("%s + \\sum_i \\sigma_i^{-2}", wz)),
                  L("m_j", "=", sprintf("%s(%s) + \\sum_i (\\lambda_i - \\log t_{ij})/\\sigma_i^2", wz, al_sum(xb("beta"), "\\gamma\\theta_j", if (qr) "k_1\\tilde\\nu_j")))),
    cross = list(L("p_j", "=", sprintf("%s + \\sum_i w_{ij}", if (C$fixed) "1" else "1/s")),
                 L("m_j", "=", sprintf("\\sum_i w_{ij}(\\lambda_i - \\rho_i\\theta_j%s - \\log t_{ij})", k1nu))))
  zeta_sec <- c(list(L("\\zeta_j \\mid \\cdot", "\\sim", "\\mathcal{N}(m_j/p_j,\\; 1/p_j)", "draw_gauss()")), zeta)
  # location move
  # prior of speed in the move: mean mz and precision (a constant variance "s", "v", "1", or the ALD weights)
  zv <- switch(eng, rtirt = "v", latent = if (qr) NA else "s", cross = if (C$fixed) "1" else "s")
  zm <- switch(eng, rtirt = if (P) "\\mu_{2j} + c(\\theta_j - \\mu_{1j})" else "c\\theta_j",
               latent = al_sum(xb("beta"), "\\gamma\\theta_j", if (qr) "k_1\\tilde\\nu_j"), cross = "0")
  sl <- C$lp[2]; dl <- if (sl == 1) "" else sprintf("/%s", fmt_num(sl^2))
  lam_m <- if (C$lp[1] == 0) sprintf("-\\sum_i \\lambda_i%s", dl) else sprintf("\\sum_i (%s - \\lambda_i)%s", fmt_num(C$lp[1]), dl)
  zeta_p <- if (is.na(zv)) "\\sum_j \\tilde w_j" else if (zv == "1") "N" else sprintf("N/%s", zv)
  zeta_m <- if (is.na(zv)) sprintf("\\sum_j \\tilde w_j(%s - \\zeta_j)", zm)
            else if (zm == "0") sprintf("-\\sum_j \\zeta_j%s", if (zv == "1") "" else paste0("/", zv))
            else sprintf("\\sum_j (%s - \\zeta_j)%s", zm, if (zv == "1") "" else paste0("/", zv))
  loc <- list(L("(\\lambda_i, \\zeta_j)", "\\leftarrow", "(\\lambda_i + \\delta,\\; \\zeta_j + \\delta)", "exact shift (as draw_shift()): leaves the likelihood unchanged (Liu & Sabatti, 2000); a shift that would cross the lower bound of lambda keeps the state"),
              L("\\delta \\mid \\cdot", "\\sim", "\\mathcal{N}(m/p,\\; 1/p)"),
              L("p", "=", al_sum(al_prec(sl, "K"), zeta_p)),
              L("m", "=", gsub("\\+ -", "- ", paste(lam_m, "+", zeta_m))),
              if (is.finite(C$lp[3])) L("", "", sprintf("\\text{no move if } \\min_i \\lambda_i + \\delta \\le %s", fmt_num(C$lp[3]))))
  sr <- pr$rho[["sd"]]
  ss <- if (C$fixed) "1" else "s"
  cross_move <- list(L("(\\rho_i, \\zeta_j)", "\\leftarrow", "(\\rho_i - \\delta,\\; \\zeta_j + \\delta\\theta_j)", "leaves the likelihood unchanged"),
                     L("\\delta \\mid \\cdot", "\\sim", "\\mathcal{N}(m/p,\\; 1/p)"),
                     L("p", "=", sprintf("%s + \\sum_j \\theta_j^2/%s", al_prec(sr, "K"), ss)),
                     L("m", "=", sprintf("\\sum_i \\rho_i%s - \\sum_j \\zeta_j\\theta_j/%s", if (sr == 1) "" else sprintf("/%s", fmt_num(sr^2)), ss)))
  nu <- list(L("r_{ij}", "=", r_ij), L("\\nu_{ij}^{-1} \\mid \\cdot", "\\sim", "\\mathrm{InvGauss}(\\sqrt{\\psi_i/\\chi_{ij}},\\; \\psi_i)", "draw_nu()"),
             L("\\chi_{ij},\\; \\psi_i", "=", "r_{ij}^2/(k_2\\sigma_i),\\; (k_1^2 + 2k_2)/(k_2\\sigma_i)"))
  rho <- list(L("\\rho_i \\mid \\cdot", "\\sim", "\\mathcal{N}(m_i/p_i,\\; 1/p_i)", "normal regression of the log times on theta"),
              L("p_i", "=", sprintf("%s + \\sum_j w_{ij}\\theta_j^2", al_prec(sr))),
              L("m_i", "=", sprintf("\\sum_j w_{ij}\\theta_j(\\lambda_i - \\zeta_j%s - \\log t_{ij})", k1nu)))
  svar_cross <- al_var("s", "u_s", "N/2", "\\sum_j \\zeta_j^2/2", pr$var_speed, pr$var_speed_family, "draw_var(): speed variance")

  # collapsed steps (zeta integrated out)
  coll_theta <- function() {
    q <- switch(eng,
      cross = list("1 + P_j + \\sum_i w_{ij}\\rho_i^2", "\\sum_i w_{ij}\\rho_i", sprintf("%s + \\sum_i w_{ij}", if (C$fixed) "1" else "1/s"),
                   "h_j + \\sum_i w_{ij}\\rho_i y_{ij}", "\\sum_i w_{ij}y_{ij}", sprintf("y_{ij} = \\lambda_i%s - \\log t_{ij}", k1nu)),
      rtirt = list("1 + P_j + c^2/v", "-c/v", "1/v + \\sum_i \\sigma_i^{-2}", if (P) "h_j + \\mu_{1j} - c\\, m_{0j}/v" else "h_j",
                   sprintf("\\sum_i (\\lambda_i - \\log t_{ij})/\\sigma_i^2%s", if (P) " + m_{0j}/v" else ""), if (P) "m_{0j} = \\mu_{2j} - c\\mu_{1j}"),
      latent = {
        m0 <- P || qr
        list(sprintf("1 + P_j + \\gamma^2 %s", wz), sprintf("-\\gamma %s", wz), sprintf("%s + \\sum_i \\sigma_i^{-2}", wz),
             if (m0) sprintf("h_j - \\gamma %s m_{0j}", wz) else "h_j",
             sprintf("\\sum_i (\\lambda_i - \\log t_{ij})/\\sigma_i^2%s", if (m0) sprintf(" + %s m_{0j}", wz) else ""),
             if (m0) sprintf("m_{0j} = %s", al_sum(xb("beta"), if (qr) "k_1\\tilde\\nu_j")))
      })
    list(irt,
         L("(\\theta_j, \\zeta_j) \\mid \\cdot", "\\sim", "\\mathcal{N}_2 \\text{ with precision } Q \\text{ and linear term } (h_\\theta, h_\\zeta)", "draw_theta_collapsed()"),
         L("Q_{\\theta\\theta},\\; Q_{\\theta\\zeta},\\; Q_{\\zeta\\zeta}", "=", paste(q[[1]], q[[2]], q[[3]], sep = ",\\; ")),
         L("h_\\theta,\\; h_\\zeta", "=", paste(q[[4]], q[[5]], sep = ",\\; "), ""),
         if (length(q) == 6 && !is.null(q[[6]])) L("", "", q[[6]]),
         L("\\theta_j \\mid \\cdot \\text{ (speed integrated out)}", "\\sim", "\\mathcal{N}(m_j/p_j,\\; 1/p_j)"),
         L("p_j,\\; m_j", "=", "Q_{\\theta\\theta} - Q_{\\theta\\zeta}^2/Q_{\\zeta\\zeta},\\; h_\\theta - Q_{\\theta\\zeta}h_\\zeta/Q_{\\zeta\\zeta}", "Schur complement"))
  }
  coll_items <- function() {
    tau <- switch(eng, cross = if (C$fixed) "1" else "s", rtirt = "v", latent = "1/\\tilde w_j")
    if (eng == "latent" && !qr) tau <- "s"
    u <- switch(eng, cross = sprintf("\\log t_j%s", if (qr) " - k_1\\nu_j" else ""),
                rtirt = sprintf("\\log t_j + (%s)\\mathbf{1}", if (P) "\\mu_{2j} + c(\\theta_j - \\mu_{1j})" else "c\\theta_j"),
                latent = sprintf("\\log t_j + (%s)\\mathbf{1}", al_sum(xb("beta"), "\\gamma\\theta_j", if (qr) "k_1\\tilde\\nu_j")))
    d <- if (eng == "cross") "w_{ij}" else "\\sigma_i^{-2}"
    list(L(if (eng == "cross") "(\\lambda, \\rho) \\mid \\theta, \\cdot" else "\\lambda \\mid \\theta, \\cdot", "\\sim",
           "\\mathcal{N}(Q^{-1}h,\\; Q^{-1})", "draw_rt_items_collapsed(): speed integrated out"),
         L("u_j", "=", sprintf("%s \\sim \\mathcal{N}_K(%s,\\; \\Sigma_j),\\; \\Sigma_j = D_j^{-1} + %s\\,\\mathbf{1}\\mathbf{1}^\\top", u,
                               if (eng == "cross") "\\lambda - \\theta_j\\rho" else "\\lambda", tau), sprintf("D_j = diag of the item precisions")),
         L("\\Sigma_j^{-1}", "=", sprintf("D_j - g_j d_j d_j^\\top,\\; d_{ij} = %s,\\; g_j = 1/(1/%s + \\sum_i d_{ij})", d, tau), "Sherman-Morrison"),
         L("Q,\\; h", "=", sprintf("\\sum_j A_j^\\top\\Sigma_j^{-1}A_j + %s,\\; \\sum_j A_j^\\top\\Sigma_j^{-1}u_j%s",
                                   if (eng == "cross") sprintf("\\mathrm{diag}(%s, %s)", al_prec(C$lp[2], "I"), al_prec(sr, "I")) else al_prec(C$lp[2], "I"),
                                   if (C$lp[1] == 0) "" else if (eng == "cross") sprintf(" + ((%s)\\mathbf{1}, 0)", al_mterm(C$lp[1], C$lp[2]))
                                   else sprintf(" + (%s)\\mathbf{1}", al_mterm(C$lp[1], C$lp[2]))),
           if (eng == "cross") "A_j = [I, -theta_j I]" else "A_j = I"),
         if (is.finite(C$lp[3])) L("", "", sprintf("\\lambda_i > %s \\text{ by rejection}", fmt_num(C$lp[3])), "after 50 failures one Gibbs sweep through the block"))
  }
  zeta_redraw <- c(list(L("\\zeta_j \\mid \\cdot", "\\sim", "\\mathcal{N}(m_j/p_j,\\; 1/p_j)", "as the full conditional of speed below; drawn right after the steps that integrate it out")), zeta)

  steps <- switch(eng,
    cross = list(
      if (qr) c("ALD mixing variables", nu),
      if (!collapse) c("Cross-relations", rho),
      if (!C$fixed) c("Speed variance", svar_cross),
      c("Accuracy items", acc),
      if (collapse) c("Ability (speed integrated out)", coll_theta()) else c("Ability", theta_sec),
      if (collapse) c("RT items (lambda, rho) jointly, speed integrated out", coll_items()) else c("Time intensities", lam),
      if (collapse) c("Speed", zeta_redraw) else c(if (qr) "ALD scales" else "Residual variances", sig),
      if (collapse) c(if (qr) "ALD scales" else "Residual variances", sig) else c("Speed", zeta_sec),
      if (!collapse) c("Location move", loc),
      if (!collapse) c("Cross-relation move", cross_move)),
    rtirt = {
      sc <- pr$cov[["sd"]]
      cov <- if (C$fixed) list(
        L("c", "=", "\\tanh(z)", "c = Corr(ability, speed), prior U(-1, 1); Var(speed | x) = 1"),
        L("z^*", "=", "z + 0.05\\,\\epsilon,\\; \\epsilon \\sim \\mathcal{N}(0, 1)", "5 random-walk Metropolis steps per iteration"),
        L("\\text{target}", "", "(1 - c^2)\\prod_j \\mathcal{N}(r_{2j};\\; c\\,e_{1j},\\; 1 - c^2)", "(1 - c^2): Jacobian of tanh"),
        L("v", "=", "1 - c^2"))
        else c(list(L("c \\mid \\cdot", "\\sim", "\\mathcal{N}(m/p,\\; 1/p)", "regression of the speed residual on the ability residual"),
                    L("p,\\; m", "=", sprintf("%s + \\sum_j e_{1j}^2/v,\\; \\sum_j e_{1j}r_{2j}/v", al_prec(sc)))),
               al_var("v", "u_v", "N/2", "\\sum_j (r_{2j} - c\\,e_{1j})^2/2", pr$var_speed, pr$var_speed_family, "draw_var(): residual speed variance"))
      list(
        if (P) c("Regression of ability", list(L("\\beta_1 \\mid \\cdot", "\\sim", "\\mathcal{N}(Q^{-1}h,\\; Q^{-1})", "multivariate normal"),
                 L("Q,\\; h", "=", sprintf("%s + (1 + c^2/v)X^\\top X,\\; X^\\top\\theta - (c/v)X^\\top(\\zeta - X\\beta_2 - c\\theta)", al_prec(sb, "I"))))),
        if (P) c("Regression of speed", list(L("\\beta_2 \\mid \\cdot", "\\sim", "\\mathcal{N}(Q^{-1}h,\\; Q^{-1})"),
                 L("Q,\\; h", "=", sprintf("%s + X^\\top X/v,\\; X^\\top(\\zeta - c\\,e_1)/v", al_prec(sb, "I"))))),
        c("Ability-speed covariance", c(list(L("e_{1j},\\; r_{2j}", "=", sprintf("\\theta_j%s,\\; \\zeta_j%s", if (P) " - \\mu_{1j}" else "", if (P) " - \\mu_{2j}" else ""))), cov)),
        c("Accuracy items", acc),
        if (collapse) c("Ability (speed integrated out)", coll_theta()) else c("Ability", theta_sec),
        if (collapse) c("Time intensities (speed integrated out)", coll_items()) else c("Time intensities", lam),
        if (collapse) c("Speed", zeta_redraw) else c("Residual variances", sig),
        if (collapse) c("Residual variances", sig) else c("Speed", zeta_sec),
        if (!collapse) c("Location move", loc))
    },
    latent = {
      nut <- list(L("r_j", "=", "\\zeta_j - z_j^\\top(\\beta, \\gamma),\\; z_j = (x_j, \\theta_j)"),
                  L("\\tilde\\nu_j^{-1} \\mid \\cdot", "\\sim", "\\mathrm{InvGauss}(\\sqrt{\\psi/\\chi_j},\\; \\psi)", "draw_nu()"),
                  L("\\chi_j,\\; \\psi", "=", "r_j^2/(k_2 s),\\; (k_1^2 + 2k_2)/(k_2 s)"))
      reg <- list(L("(\\beta, \\gamma) \\mid \\cdot", "\\sim", "\\mathcal{N}(Q^{-1}h,\\; Q^{-1})", "weighted regression of speed on z_j = (x_j, theta_j)"),
                  L("Q,\\; h", "=", sprintf("%s + \\sum_j %s z_j z_j^\\top,\\; \\sum_j %s z_j(\\zeta_j%s)", al_prec(sb, "I"), wz, wz, if (qr) " - k_1\\tilde\\nu_j" else "")))
      vz <- if (qr) "1/\\tilde w_j" else "s"                     # variance of speed given z_i
      reg_c <- list(L("\\hat\\zeta_j", "=", "\\sum_i (\\lambda_i - \\log t_{ij})\\sigma_i^{-2}/P_r,\\; P_r = \\sum_i \\sigma_i^{-2}", "RT part of speed as a pseudo-observation"),
                    L("\\hat\\zeta_j", "\\sim", sprintf("\\mathcal{N}(z_j^\\top(\\beta, \\gamma)%s,\\; %s + 1/P_r)", if (qr) " + k_1\\tilde\\nu_j" else "", vz), "speed integrated out"),
                    L("(\\beta, \\gamma) \\mid \\cdot", "\\sim", "\\mathcal{N}(Q^{-1}h,\\; Q^{-1})"),
                    L("Q,\\; h", "=", sprintf("%s + \\sum_j \\hat w_j z_j z_j^\\top,\\; \\sum_j \\hat w_j z_j%s,\\; \\hat w_j = 1/(%s + 1/P_r)", al_prec(sb, "I"), if (qr) "(\\hat\\zeta_j - k_1\\tilde\\nu_j)" else "\\hat\\zeta_j", vz)),
                    L("\\zeta_j \\mid \\cdot", "\\sim", sprintf("\\mathcal{N}((%s m_j + P_r\\hat\\zeta_j)/(%s + P_r),\\; 1/(%s + P_r)),\\; m_j = %s", wz, wz, wz,
                                                              al_sum("z_j^\\top(\\beta, \\gamma)", if (qr) "k_1\\tilde\\nu_j")), "then speed from its full conditional"))
      svar <- if (qr) al_var("s", "u_s", "3N/2", "\\sum_j (r_j - k_1\\tilde\\nu_j)^2/(2k_2\\tilde\\nu_j) + \\sum_j \\tilde\\nu_j", pr$var_speed, pr$var_speed_family, "draw_var(): ALD scale of speed")
              else al_var("s", "u_s", "N/2", "\\sum_j r_j^2/2", pr$var_speed, pr$var_speed_family, "draw_var(): residual speed variance")
      list(
        if (qr) c("ALD mixing variables of speed", nut),
        if (collapse) c("Structural regression (speed integrated out)", reg_c) else c("Structural regression", reg),
        c(if (qr) "ALD scale of speed" else "Residual speed variance", c(if (!qr) list(L("r_j", "=", "\\zeta_j - z_j^\\top(\\beta, \\gamma),\\; z_j = (x_j, \\theta_j)")), svar)),
        c("Accuracy items", acc),
        if (collapse) c("Ability (speed integrated out)", coll_theta()) else c("Ability", theta_sec),
        if (collapse) c("Time intensities (speed integrated out)", coll_items()) else c("Time intensities", lam),
        if (collapse) c("Speed", zeta_redraw) else c("Residual variances", sig),
        if (collapse) c("Residual variances", sig) else c("Speed", zeta_sec),
        if (!collapse) c("Location move", loc))
    })
  c(list(al_settings(C, collapse), notation), al_number(steps))
}

# number the steps: each element is c(title, line, line, ...)
al_number <- function(steps) {
  steps <- Filter(Negate(is.null), steps)
  lapply(seq_along(steps), function(k) list(title = sprintf("Step %d. %s", k, steps[[k]][[1]]), lines = Filter(Negate(is.null), steps[[k]][-1])))
}

al_settings <- function(C, collapse) {
  f <- C$fit; cond <- C$cond
  run <- if (!is.null(f)) sprintf("%s chains x %s iterations, burn-in %s, thinning %s", cond$n_chain, cond$n_iter, cond$n_burnin, cond$n_thin) else "not sampled yet"
  al_sec("Gibbs sampler",
    L("\\text{one iteration}", "", "\\text{the steps below, in this order}", run),
    L("\\text{sampler}", "", if (collapse) "\\text{partially collapsed Gibbs (collapse = TRUE)}" else "\\text{Gibbs (collapse = FALSE)}",
      if (collapse) "speed integrated out where marked and redrawn before it is used (van Dyk & Park, 2008)" else ""))
}

# ---- EM ------------------------------------------------------------------------------------------
al_em <- function(C, object) {
  eng <- C$eng; P <- C$P; rt <- C$rt; f <- object$ecm
  map <- isTRUE(f$map); Q <- f$Q %||% 41
  pr <- C$pr
  head <- al_sec("EM algorithm (Bock & Aitkin, 1981)",
    L("\\text{fit}", "", if (is.null(f)) "\\text{not fitted yet}" else sprintf("\\text{%d iterations, %s, log-likelihood %s}", f$iterations,
      if (isTRUE(f$converged)) "converged" else "not converged", fmt_num(f$logLik)), if (map) "prior = \"map\": priors on the item parameters" else "maximum marginal likelihood"),
    L("\\theta_{jq}", "=", sprintf("%sz_q,\\; q = 1, \\ldots, %d", if (P && eng %in% c("mlirt", "rtirt")) sprintf("x_j^\\top\\beta%s + ", if (eng == "rtirt") "_1" else "") else "", Q),
      "Gauss-Hermite nodes z_q and weights g_q of N(0, 1)"))
  ll_rt <- if (rt) L("\\ell_{jq}", "\\leftarrow", "\\ell_{jq} + \\log \\mathcal{N}_K(\\log t_j \\mid \\theta_{jq})", "speed integrated in closed form (normal given ability)")
  speed <- if (rt) list(
    L("\\zeta_j \\mid \\theta_{jq}, t_j", "\\sim", "\\mathcal{N}(M_{jq},\\; 1/P_\\zeta),\\; P_\\zeta = 1/v_\\zeta + \\sum_i \\sigma_i^{-2}", "posterior of speed at each node"),
    L("M_{jq}", "=", {
      m <- switch(eng, rtirt = if (P) "(\\mu_{2j} + c(\\theta_{jq} - \\mu_{1j}))" else "c\\theta_{jq}",
                  latent = if (P) "(x_j^\\top\\beta + \\gamma\\theta_{jq})" else "\\gamma\\theta_{jq}", cross = NULL)
      sprintf("[%s\\sum_i (\\log t_{ij} - \\lambda_i%s)/\\sigma_i^2]/P_\\zeta", if (is.null(m)) "-" else sprintf("%s/v_\\zeta - ", m),
              if (eng == "cross") " + \\rho_i\\theta_{jq}" else "")
    },
      sprintf("v_zeta = %s", switch(eng, rtirt = "v", latent = "s", cross = if (C$fixed) "1" else "s"))))
  estep <- al_sec("E-step",
    L("\\ell_{jq}", "=", "\\log g_q + \\sum_i [y_{ij}\\eta_{ijq} - \\log(1 + \\exp \\eta_{ijq})],\\; \\eta_{ijq} = a_i(\\theta_{jq} - b_i)", if (P && eng %in% c("mlirt", "rtirt")) "calc_loglik_nodes(): nodes differ between persons (C++ kernel)" else "calc_loglik_nodes(): common nodes, two matrix products (Bock-Aitkin collapsing)"),
    ll_rt,
    L("\\pi_{jq}", "=", "\\exp \\ell_{jq}/\\sum_{q'} \\exp \\ell_{jq'}", "calc_node_weights()"),
    if (rt) speed[[1]], if (rt) speed[[2]],
    L("\\text{moments}", "", if (rt) "E[\\theta_j],\\; E[\\theta_j^2],\\; E[\\zeta_j],\\; E[\\zeta_j^2],\\; E[\\theta_j\\zeta_j]" else "E[\\theta_j],\\; E[\\theta_j^2]", "sums over the nodes with weights pi_jq"))
  pa <- pr$a; pb <- pr$b
  items <- al_sec("M-step: accuracy items (Polya-Gamma minorizer)",
    L("\\bar\\omega_{ijq}", "=", "\\tanh(\\eta_{ijq}/2)/(2\\eta_{ijq})", "E[omega] of PG(1, eta): quadratic minorizer of the logistic log-likelihood (Jaakkola & Jordan, 2000)"),
    L("S_{ri}", "=", "\\sum_j\\sum_q \\pi_{jq}\\bar\\omega_{ijq}\\theta_{jq}^r,\\; r = 0, 1, 2", "calc_item_sums(); with common nodes sum_q n_iq omega_iq z_q^r, n_iq = sum_j pi_jq"),
    L("T_{0i},\\; T_{1i}", "=", "\\sum_j \\kappa_{ij},\\; \\sum_j \\kappa_{ij}E[\\theta_j]"),
    if (!map && !C$one_pl) L("(a_i, a_i b_i)", "=", "(S_{0i}T_{1i} - S_{1i}T_{0i},\\; S_{1i}T_{1i} - S_{2i}T_{0i})/(S_{2i}S_{0i} - S_{1i}^2)",
                             "update_items_logistic(): 2 x 2 weighted least squares; step halving keeps a_i > 0"),
    if (!map && C$one_pl) L("b_i", "=", "(S_{1i} - T_{0i})/S_{0i}", "1PL"),
    if (map) L("b_i", "=", sprintf("\\frac{a_i(a_i S_{1i} - T_{0i})%s}{%s + a_i^2 S_{0i}}%s", if (pb[["mean"]] != 0) sprintf(" + %s", al_mterm(pb[["mean"]], pb[["sd"]])) else "",
                                   al_prec(pb[["sd"]]), if (is.finite(pb[["lower"]]) || is.finite(pb[["upper"]])) sprintf(" \\text{ clamped to } [%s, %s]", fmt_num(pb[["lower"]]), fmt_num(pb[["upper"]])) else ""),
              "conditional mode with the prior of b"),
    if (map && !C$one_pl) L("a_i", "=", sprintf("\\max\\{(%s)/(%s),\\; 0.001\\}", al_sum(al_mterm(pa[1], pa[2]), "T_{1i} - b_iT_{0i}"),
                                               al_sum(al_prec(pa[2]), "S_{2i} - 2b_iS_{1i} + b_i^2S_{0i}")), "conditional mode with the N+ prior of a"),
    L("", "", "\\text{the log-likelihood never decreases (MM step inside EM)}"))
  rtm <- if (rt) {
    ss <- L("S_i", "=", sprintf("\\sum_j E[(\\log t_{ij} + \\zeta_j - \\lambda_i%s)^2]", if (eng == "cross") " + \\rho_i\\theta_j" else ""))
    lam <- if (eng == "cross") L("(\\lambda_i, \\rho_i)", "", "\\text{2 x 2 least squares of } E[\\log t_{ij} + \\zeta_j] \\text{ on } (1, -E[\\theta_j])",
                                  if (map) "with the priors of lambda and rho, weighted by 1/sigma_i^2" else "expected normal equations")
           else L("\\lambda_i", "=", if (map) sprintf("[%s]/(N/\\sigma_i^2 + %s)", al_sum("\\sum_j(\\log t_{ij} + E[\\zeta_j])/\\sigma_i^2", al_mterm(C$lp[1], C$lp[2])), al_prec(C$lp[2]))
                                     else "\\sum_j(\\log t_{ij} + E[\\zeta_j])/N")
    sig <- if (!map) L("\\sigma_i^2", "=", "S_i/N")
           else if (pr$sigma2t_family == "inverse-gamma") L("\\sigma_i^2", "=", sprintf("(%s + S_i/2)/(%s + N/2 + 1)", fmt_num(pr$sigma2t[["scale"]]), fmt_num(pr$sigma2t[["shape"]])), "mode with the IG prior")
           else L("\\sigma_i^2", "=", sprintf("\\arg\\max_v\\; -\\tfrac{N}{2} \\log v - S_i/(2v) + \\log p(v)"), sprintf("one-dimensional search; half-t(%s, %s) prior on sigma_i", fmt_num(pr$sigma2t[["df"]]), fmt_num(pr$sigma2t[["scale"]])))
    al_sec("M-step: response-time items", lam, ss, sig)
  }
  st <- switch(eng,
    mlirt = if (P) al_sec("M-step: latent regression", L("\\beta", "=", "(X^\\top X)^{-1}X^\\top E[\\theta]")),
    rtirt = al_sec("M-step: person model",
      if (P) L("\\beta_1", "=", "(X^\\top X)^{-1}X^\\top E[\\theta]"),
      L(sprintf("(%s c)", if (P) "\\beta_2 - c\\beta_1,\\; " else ""), "", "\\text{expected normal equations of the regression of } \\zeta \\text{ on } (x, \\theta)"),
      if (C$fixed) L("c", "=", "\\arg\\max_c\\; -\\tfrac{N}{2} \\log(1 - c^2) - \\mathrm{SS}(c)/(2(1 - c^2)),\\; v = 1 - c^2", "one-dimensional search (Var(speed | x) = 1)")
      else L("v", "=", "E[\\text{residual sum of squares}]/N")),
    latent = al_sec("M-step: structural regression",
      L("(\\beta, \\gamma)", "", "\\text{expected normal equations of the regression of } \\zeta \\text{ on } (x, \\theta)"),
      L("s", "=", "E[\\text{residual sum of squares}]/N")),
    cross = if (!C$fixed) al_sec("M-step: speed variance", L("s", "=", "\\sum_j E[\\zeta_j^2]/N")))
  conv <- al_sec("Acceleration, convergence and standard errors",
    L("\\text{SQUAREM}", "", "x \\leftarrow x - 2\\alpha d_1 + \\alpha^2 d_2,\\; \\alpha = \\min(-\\|d_1\\|/\\|d_2\\|,\\; -1)", "em_squarem(): SqS3 (Varadhan & Roland, 2008); a plain EM step when it would lower the log-likelihood"),
    L("\\text{stop}", "", "|\\Delta\\ell| \\le 10^{-10}|\\ell| \\text{ and } \\max|\\Delta x| \\le 10^{-5}", "x: log a, b, lambda, log sigma2, ... (unconstrained)"),
    L("\\mathrm{SE}", "", "\\mathrm{diag}(-H^{-1})^{1/2},\\; H = \\partial s(x)/\\partial x", "calc_hessian(): central differences of the analytic score s (Fisher identity), then the delta method"))
  list(head, estep, items, rtm, st, conv)
}
