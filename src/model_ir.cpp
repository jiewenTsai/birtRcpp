// Multilinear expressions of declared models (block_model()): a sum of monomials
// coef * x1 * x2 * ..., where every x is a number per model, per person (rows), per item
// (columns) or per cell. The value is accumulated cell by cell, without expanding any vector to
// a persons x items matrix and without temporaries.
#include <RcppArmadillo.h>
#include "rng.h"
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]

// dims: 0 = one value, 1 = per person (length nr), 2 = per item (length nc), 3 = per cell (nr x nc)
// [[Rcpp::export(.ml_eval)]]
Rcpp::NumericMatrix ml_eval(const Rcpp::NumericVector& coef, const Rcpp::List& ids, const Rcpp::List& vals,
                            const Rcpp::IntegerVector& dims, int nr, int nc) {
  for (int k = 0; k < vals.size(); ++k) {                      // every value must have the length its dimension says
    const R_xlen_t need = dims[k] == 0 ? 1 : dims[k] == 1 ? nr : dims[k] == 2 ? nc : (R_xlen_t) nr * nc;
    if (Rf_xlength(vals[k]) != need) Rcpp::stop("internal: a value of a compiled expression has length %d, %d needed", (int) Rf_xlength(vals[k]), (int) need);
  }
  Rcpp::NumericMatrix out(nr, nc);
  double* o = out.begin();
  const int nm = coef.size();
  std::vector<double> pv(nr), iv(nc);
  for (int m = 0; m < nm; ++m) {
    Rcpp::IntegerVector id = ids[m];
    double sc = coef[m];
    bool has_p = false, has_i = false;
    std::vector<const double*> cells;
    for (int k = 0; k < id.size(); ++k) {
      Rcpp::NumericVector v = vals[id[k]];
      const double* x = v.begin();
      switch (dims[id[k]]) {
      case 0: sc *= x[0]; break;
      case 1: if (has_p) { for (int j = 0; j < nr; ++j) pv[j] *= x[j]; } else { std::copy(x, x + nr, pv.begin()); has_p = true; } break;
      case 2: if (has_i) { for (int i = 0; i < nc; ++i) iv[i] *= x[i]; } else { std::copy(x, x + nc, iv.begin()); has_i = true; } break;
      default: cells.push_back(x);
      }
    }
    for (int i = 0; i < nc; ++i) {
      const double ci = sc * (has_i ? iv[i] : 1.0);
      if (ci == 0.0) continue;
      double* oc = o + (size_t) i * nr;
      const size_t off = (size_t) i * nr;
      if (cells.empty()) {
        if (has_p) for (int j = 0; j < nr; ++j) oc[j] += ci * pv[j];
        else for (int j = 0; j < nr; ++j) oc[j] += ci;
      } else {
        for (int j = 0; j < nr; ++j) {
          double t = ci * (has_p ? pv[j] : 1.0);
          for (const double* c : cells) t *= c[off + j];
          oc[j] += t;
        }
      }
    }
  }
  return out;
}

// A compiled expression bound to values, evaluated one cell at a time.
struct MLExpr {
  int nr = 0;
  std::vector<double> coef;
  std::vector<std::vector<int>> dim;
  std::vector<std::vector<const double*>> ptr;
  bool empty = true;
  MLExpr() {}
  MLExpr(const Rcpp::List& pk, int nr_, int nc_) : nr(nr_) {
    if (pk.size() == 0) return;
    Rcpp::NumericVector cf = pk[0]; Rcpp::List ids = pk[1], vals = pk[2]; Rcpp::IntegerVector dims = pk[3];
    for (int k = 0; k < vals.size(); ++k) {
      const R_xlen_t need = dims[k] == 0 ? 1 : dims[k] == 1 ? nr : dims[k] == 2 ? (R_xlen_t) nc_ : (R_xlen_t) nr * nc_;
      if (Rf_xlength(vals[k]) != need) Rcpp::stop("internal: a value of a compiled expression has length %d, %d needed", (int) Rf_xlength(vals[k]), (int) need);
    }
    for (int m = 0; m < cf.size(); ++m) {
      coef.push_back(cf[m]);
      Rcpp::IntegerVector id = ids[m];
      std::vector<int> d; std::vector<const double*> p;
      for (int k = 0; k < id.size(); ++k) { Rcpp::NumericVector v = vals[id[k]]; d.push_back(dims[id[k]]); p.push_back(v.begin()); }
      dim.push_back(d); ptr.push_back(p);
    }
    empty = false;
  }
  inline double at(int j, int i) const {
    double s = 0.0;
    for (size_t m = 0; m < coef.size(); ++m) {
      double t = coef[m];
      const std::vector<int>& d = dim[m]; const std::vector<const double*>& p = ptr[m];
      for (size_t k = 0; k < d.size(); ++k)
        t *= (d[k] == 0) ? p[k][0] : (d[k] == 1) ? p[k][j] : (d[k] == 2) ? p[k][i] : p[k][(size_t) i * nr + j];
      s += t;
    }
    return s;
  }
};

// The normal part (precision, linear term) of one parameter from one part of the model, summed to
// the dimension of the parameter (red: 0 = one value, 1 = per row, 2 = per column), without
// temporaries. c1 = d(expression)/d(parameter), r0 = the expression at parameter = 0.
// mode 0: normal part,   weight 1/sd2 (third expression), P += w c1^2,  H -= w c1 r0
// mode 1: logistic part, omega given (Gibbs),             P += o c1^2,  H += c1 (kappa - o r0)
// mode 2: logistic part, E[omega] at eta (third expr.),   as mode 1 with o = pg_mean(eta)
// mode 3: sums for a variance: P += weight, H += weight r0^2 (r0 = the residual)
// mask and row weights multiply every cell (NULL: 1); cells with weight below 1e-13 (far quadrature
// nodes) are skipped. keep_om (mode 2): also return omega = mask * E[omega] per cell, so that the next
// parameter of the same part can use it in mode 1 (one minorizer per M-step, Polson et al., 2013).
// [[Rcpp::export(.ml_part)]]
Rcpp::List ml_part(const Rcpp::List& c1, const Rcpp::List& r0, const Rcpp::List& third, SEXP om_, SEXP mask_, SEXP kappa_, SEXP roww_,
                   int mode, int red, int nr, int nc, bool keep_om = false) {
  MLExpr C(c1, nr, nc), R(r0, nr, nc), T(third, nr, nc);
  const R_xlen_t ncell = (R_xlen_t) nr * nc;                  // cell-wise inputs have one value per cell, row weights one per row
  if (!Rf_isNull(om_) && Rf_xlength(om_) != ncell) Rcpp::stop("internal: omega has length %d, %d needed", (int) Rf_xlength(om_), (int) ncell);
  if (!Rf_isNull(mask_) && Rf_xlength(mask_) != ncell) Rcpp::stop("internal: the mask has length %d, %d needed", (int) Rf_xlength(mask_), (int) ncell);
  if (!Rf_isNull(kappa_) && Rf_xlength(kappa_) != ncell) Rcpp::stop("internal: kappa has length %d, %d needed", (int) Rf_xlength(kappa_), (int) ncell);
  if (!Rf_isNull(roww_) && Rf_xlength(roww_) != nr) Rcpp::stop("internal: the row weights have length %d, %d needed", (int) Rf_xlength(roww_), nr);
  const double* om = Rf_isNull(om_) ? nullptr : REAL(om_);
  const double* mk = Rf_isNull(mask_) ? nullptr : REAL(mask_);
  const double* ka = Rf_isNull(kappa_) ? nullptr : REAL(kappa_);
  const double* rw = Rf_isNull(roww_) ? nullptr : REAL(roww_);
  const int n = red == 0 ? 1 : (red == 1 ? nr : nc);
  Rcpp::NumericVector P(n), H(n);
  const bool keep = keep_om && mode == 2;
  Rcpp::NumericVector omv(keep ? ncell : 0);
  for (int i = 0; i < nc; ++i) for (int j = 0; j < nr; ++j) {
    const size_t cell = (size_t) i * nr + j;
    double w = (mk ? mk[cell] : 1.0) * (rw ? rw[j] : 1.0);
    if (w < 1e-13) continue;
    const int k = red == 0 ? 0 : (red == 1 ? j : i);
    if (mode == 3) { const double r = R.at(j, i); P[k] += w; H[k] += w * r * r; continue; }
    const double c = C.at(j, i), r = R.at(j, i);
    if (mode == 0) { const double ww = w / T.at(j, i); P[k] += ww * c * c; H[k] -= ww * c * r; }
    else {
      double o;
      if (mode == 1) o = (om ? om[cell] : 0.0) * (rw ? rw[j] : 1.0);
      else { const double pg = pg_mean1(T.at(j, i)); o = w * pg; if (keep) omv[cell] = (mk ? mk[cell] : 1.0) * pg; }
      const double kap = (ka ? ka[cell] : 0.0) * (rw ? rw[j] : 1.0);
      P[k] += o * c * c; H[k] += c * kap - o * c * r;
    }
  }
  if (keep) return Rcpp::List::create(P, H, omv);
  return Rcpp::List::create(P, H);
}
