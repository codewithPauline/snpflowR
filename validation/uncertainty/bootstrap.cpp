// Research-only acceleration for protocol v1; not part of the R package API.
// Independent implementation of the documented biallelic statistic.
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

struct Component {
  double ua = 0, va = 0, ub = 0, vb = 0;
  int used = 0, low = 0, no_shared = 0, shared_fixed = 0;
};

struct Genotypes {
  int n, k;
  std::vector<int> x;
  Genotypes(NumericMatrix g) : n(g.nrow()), k(g.ncol()), x(g.size()) {
    if (n < 2 || k < 1) stop("Need at least two individuals and one block.");
    for (int i = 0; i < g.size(); ++i) {
      if (NumericVector::is_na(g[i])) x[i] = -1;
      else if (g[i] == 0 || g[i] == 1 || g[i] == 2) x[i] = int(g[i]);
      else stop("Only diploid genotypes 0, 1, 2, or NA are supported.");
    }
  }
};

std::vector<Component> components(const Genotypes& a, const Genotypes& b,
                                const std::vector<int>& wa, const std::vector<int>& wb) {
  std::vector<Component> out(a.k);
  for (int k = 0; k < a.k; ++k) {
    double aa = 0, ca = 0, ab = 0, cb = 0;
    for (int i = 0; i < a.n; ++i) {
      int v = a.x[k * a.n + i];
      if (v >= 0 && wa[i]) { aa += wa[i] * v; ca += 2 * wa[i]; }
    }
    for (int i = 0; i < b.n; ++i) {
      int v = b.x[k * b.n + i];
      if (v >= 0 && wb[i]) { ab += wb[i] * v; cb += 2 * wb[i]; }
    }
    Component& z = out[k];
    if (ca < 2 || cb < 2) { z.low = 1; continue; }
    double p = aa / ca, q = ab / cb;
    if (p == q && (p == 0 || p == 1)) { z.shared_fixed = 1; continue; }
    double sa = std::sqrt(p * q), sr = std::sqrt((1-p) * (1-q));
    if (sa + sr == 0) { z.no_shared = 1; continue; }
    double r = sa / (sa + sr);
    z.ua = (p-r)*(p-r)/2; z.va = (p+r)*(1-(p+r)/2);
    z.ub = (q-r)*(q-r)/2; z.vb = (q+r)*(1-(q+r)/2);
    z.used = 1;
  }
  return out;
}

NumericVector aggregate_components(const std::vector<Component>& z, const std::vector<int>& w) {
  double ua=0, va=0, ub=0, vb=0, used=0, low=0, ns=0, sf=0;
  for (unsigned k = 0; k < z.size(); ++k) {
    ua += w[k]*z[k].ua; va += w[k]*z[k].va;
    ub += w[k]*z[k].ub; vb += w[k]*z[k].vb;
    used += w[k]*z[k].used; low += w[k]*z[k].low;
    ns += w[k]*z[k].no_shared; sf += w[k]*z[k].shared_fixed;
  }
  double ga = va > 0 ? ua/va : NA_REAL, gb = vb > 0 ? ub/vb : NA_REAL;
  return NumericVector::create(_["gst_a"] = ga, _["gst_b"] = gb,
    _["delta"] = gb-ga, _["used"] = used, _["low"] = low,
    _["no_shared"] = ns, _["shared_fixed"] = sf);
}

void draw_weights(std::vector<int>& w, int size) {
  std::fill(w.begin(), w.end(), 0);
  for (int i = 0; i < size; ++i) ++w[int(unif_rand() * w.size())];
}

double contrast(const std::vector<Component>& z, const std::vector<int>& w) {
  double ua=0, va=0, ub=0, vb=0;
  for (unsigned k = 0; k < z.size(); ++k) {
    ua += w[k]*z[k].ua; va += w[k]*z[k].va;
    ub += w[k]*z[k].ub; vb += w[k]*z[k].vb;
  }
  return va > 0 && vb > 0 ? ub/vb - ua/va : NA_REAL;
}

void check_weights(IntegerVector w, int n) {
  if (w.size() != n) stop("Weight length mismatch.");
  bool positive = false;
  for (int i = 0; i < n; ++i) {
    if (w[i] == NA_INTEGER || w[i] < 0) stop("Weights must be nonnegative integers.");
    positive = positive || w[i] > 0;
  }
  if (!positive) stop("At least one weight must be positive.");
}

// [[Rcpp::export]]
NumericVector uncertainty_weighted_cpp(NumericMatrix ga, NumericMatrix gb,
                                      IntegerVector wa, IntegerVector wb, IntegerVector wk) {
  Genotypes a(ga), b(gb);
  if (a.k != b.k) stop("Block counts must match.");
  check_weights(wa, a.n); check_weights(wb, b.n); check_weights(wk, a.k);
  return aggregate_components(components(a, b, as<std::vector<int>>(wa), as<std::vector<int>>(wb)),
                              as<std::vector<int>>(wk));
}

// [[Rcpp::export]]
List uncertainty_bootstrap_cpp(NumericMatrix ga, NumericMatrix gb, int copies, int n_boot) {
  Genotypes a(ga), b(gb);
  if (a.k != b.k) stop("Block counts must match.");
  if ((copies != 1 && copies != 4) || n_boot < 1) stop("Unsupported protocol dimensions.");
  int k = a.k;
  std::vector<int> wa(a.n, 1), wb(b.n, 1), ones(k, 1), blocks(k), loci(k);
  std::vector<Component> original = components(a, b, wa, wb);
  NumericVector point = aggregate_components(original, ones);
  NumericMatrix draws(n_boot, 4);
  colnames(draws) = CharacterVector::create("locus", "block", "individual", "crossed");
  for (int r = 0; r < n_boot; ++r) {
    if (r % 64 == 0) checkUserInterrupt();
    draw_weights(wa, a.n); draw_weights(wb, b.n);
    draw_weights(blocks, k);
    if (copies == 1) loci = blocks;
    else draw_weights(loci, k * copies);
    std::vector<Component> sampled = components(a, b, wa, wb);
    draws(r, 0) = contrast(original, loci);
    draws(r, 1) = contrast(original, blocks);
    draws(r, 2) = contrast(sampled, ones);
    draws(r, 3) = contrast(sampled, blocks);
  }
  return List::create(_["point"] = point, _["draws"] = draws);
}
