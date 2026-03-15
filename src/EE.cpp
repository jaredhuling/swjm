// EE.cpp
// Expose am1(), am2(), temLog(), temLog2(), reRate() to R with Rcpp/RcppArmadillo.
// Also exposes fast JFM EE functions: jfm_s0s1_cpp(), jfm_s0s1_cf_cpp(),
// jfm_score_cpp().

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <algorithm>   // std::lower_bound
#include <iterator>    // std::distance
#include <set>
#include <utility>

// type alias used by am1(), am2(), and reRate()
using cmp_par = std::pair<double, arma::uword>;

// -------------------- am1 --------------------
// [[Rcpp::export(rng = false)]]
arma::rowvec am1(const arma::vec& a,
                 const arma::vec& T,
                 const arma::vec& Y,
                 const arma::vec& W,
                 const arma::mat& X,
                 const arma::vec& m) {

  arma::uword const nm = accu(m);
  arma::uword const n  = X.n_rows;
  arma::uword const p  = X.n_cols;
  arma::vec m2 = cumsum(m);

  arma::mat Xi(nm, p, arma::fill::zeros);
  arma::vec Yi(nm, arma::fill::zeros);
  arma::vec T0 = log(Y) + X * a;
  arma::uword const mn = m.n_elem;
  for (arma::uword i = 0; i < mn; i++) {
    if (i == 0 && m(i) > 0) {
      Yi.subvec(0, m2(i) - 1).fill(Y(i));
      Xi.submat(0, 0, m2(i) - 1, p - 1) = repmat(X.row(i), m(i), 1);
    }
    if (i > 0 && m(i) > 0) {
      Yi.subvec(m2(i - 1), m2(i) - 1).fill(Y(i));
      Xi.submat(m2(i - 1), 0, m2(i) - 1, p - 1) = repmat(X.row(i), m(i), 1);
    }
  }
  arma::vec texa = log(T) + Xi * a;
  arma::vec yexa = log(Yi) + Xi * a;
  arma::vec Lam(n, arma::fill::zeros);
  arma::vec de(nm, arma::fill::zeros);
  arma::uvec const idx = arma::sort_index(texa);
  auto cmp = [](cmp_par const &x, cmp_par const &y){
    return x.first <= y.first;
  };
  std::set<cmp_par, decltype(cmp)> indices(cmp);
  double w_sum{};

  {
    auto const idx_i = idx[0];
    indices.emplace(yexa[idx_i], idx_i);
    w_sum += 1;
    de(idx_i) = w_sum;
  }
  auto indices_head = indices.begin();
  for(arma::uword i = 1; i < nm; ++i) {
    auto const idx_i = idx[i];
    indices.emplace(yexa[idx_i], idx_i);
    if(yexa[idx_i] > indices_head->first)
      while(indices_head->first < texa[idx_i]){
        w_sum -= 1;
        ++indices_head;
      }
    else
      --indices_head;
    w_sum += 1;
    de(idx_i) = w_sum;
    for(arma::uword j = 1; j <= i; ++j) {
      if (T[idx_i] == T(idx[i - j]))
        de(idx[i - j]) = w_sum;
      else break;
    }
  }

  // --- Optimised suffix-cumsum approach: O((n+nm) log nm) ---
  arma::vec inv_de(nm);
  for (arma::uword i = 0; i < nm; i++)
    inv_de(i) = (de(i) > 0.0) ? 1.0 / de(i) : 0.0;

  // Sort texa ascending (reuse idx which is already argsort(texa) ascending)
  arma::vec texa_sorted = texa(idx);
  arma::vec inv_de_sorted = inv_de(idx);
  // Suffix cumulative sum: suf_cum[j] = sum_{l=j}^{nm-1} inv_de_sorted[l]
  arma::vec suf_cum = arma::flipud(arma::cumsum(arma::flipud(inv_de_sorted)));

  for (arma::uword k = 0; k < n; k++) {
    auto it = std::lower_bound(texa_sorted.begin(), texa_sorted.end(), T0[k]);
    arma::uword pos = (arma::uword)std::distance(texa_sorted.begin(), it);
    Lam[k] = (pos < nm) ? suf_cum(pos) : 0.0;
  }

  Lam = exp(-Lam);

  arma::vec R = m / Lam;

  return ((W % R - mean(W % R))).t() * X / n;
}

// -------------------- temLog --------------------
// [[Rcpp::export(rng = false)]]
arma::rowvec temLog(const arma::vec& a,
                    const arma::vec& b,
                    const arma::mat& X,
                    const arma::vec& Y,
                    const arma::vec& Z,
                    const arma::vec& D,
                    const arma::vec& W) {
  int n = Y.n_elem;
  int p = X.n_cols;
  arma::vec yexa = Y % exp(X * a);
  arma::vec ebaxZ = Z % exp(X * (b - a));
  arma::uvec ind = stable_sort_index(yexa, "descend");
  arma::vec ordD = D(ind);
  arma::vec ordW = W(ind);
  arma::mat xz = X % repmat(ebaxZ, 1, p);
  xz = cumsum(xz.rows(ind), 0);
  arma::mat c1 = X.rows(ind);
  arma::vec tmp = cumsum(ebaxZ(ind));
  arma::mat r = c1 - xz / repmat(tmp, 1, p);
  r.replace(arma::datum::nan, 0);
  return sum(repmat(ordW % ordD, 1, p) % r, 0) / n;
}

// -------------------- temLog2 --------------------
// [[Rcpp::export(rng = false)]]
arma::rowvec temLog2(const arma::vec& yexa,
                    const arma::mat& X,
                    const arma::vec& Z,
                    const arma::vec& D,
                    const arma::vec& W) {
  int n = yexa.n_elem;
  int p = X.n_cols;

  arma::uvec ind = stable_sort_index(yexa, "descend");
  arma::vec ordD = D(ind);
  arma::vec ordW = W(ind);
  arma::mat xz = X % repmat(Z, 1, p);
  xz = cumsum(xz.rows(ind), 0);
  arma::mat c1 = X.rows(ind);
  arma::vec tmp = cumsum(Z(ind));
  arma::mat r = c1 - xz / repmat(tmp, 1, p);
  r.replace(arma::datum::nan, 0);
  return sum(repmat(ordW % ordD, 1, p) % r, 0) / n;
}

// -------------------- reRate --------------------
// [[Rcpp::export(rng = false)]]
arma::vec reRate(const arma::vec& T,
                 const arma::vec& Y,
                 const arma::vec& W,
                 const arma::vec& T0) {
  arma::uword const n = Y.n_elem;
  arma::uword const m = T0.n_elem;
  arma::vec out(m, arma::fill::zeros);
  arma::vec de(n, arma::fill::zeros);

  arma::uvec const idx = arma::sort_index(T);
  auto cmp = [](cmp_par const &x, cmp_par const &y){
    return x.first <= y.first;
  };
  std::set<cmp_par, decltype(cmp)> indices(cmp);
  double w_sum{};

  {
    auto const idx_i = idx[0];
    indices.emplace(Y[idx_i], idx_i);
    w_sum += W(idx_i);
    de(idx_i) = w_sum;
  }
  auto indices_head = indices.begin();
  for(arma::uword i = 1; i < n; ++i) {
    auto const idx_i = idx[i];
    indices.emplace(Y[idx_i], idx_i);
    if(Y[idx_i] > indices_head->first)
      while(indices_head->first < T[idx_i]){
        w_sum -= W(indices_head->second);
        ++indices_head;
      }
    else
      --indices_head;
    w_sum += W(idx_i);
    de(idx_i) = w_sum;
    for(arma::uword j = 1; j <= i; ++j) {
      if (T[idx_i] == T(idx[i - j]))
        de(idx[i - j]) = w_sum;
      else break;
    }
  }

  // --- Optimised suffix-cumsum approach: O((n+m) log n) ---
  arma::vec g(n);
  for (arma::uword i = 0; i < n; i++)
    g(i) = (de(i) > 0.0) ? W(i) / de(i) : 0.0;

  // idx is already sort_index(T) ascending
  arma::vec T_sorted = T(idx);
  arma::vec g_sorted = g(idx);
  arma::vec suf_cum = arma::flipud(arma::cumsum(arma::flipud(g_sorted)));

  for (arma::uword k = 0; k < m; k++) {
    auto it = std::lower_bound(T_sorted.begin(), T_sorted.end(), T0[k]);
    arma::uword pos = (arma::uword)std::distance(T_sorted.begin(), it);
    out[k] = (pos < n) ? suf_cum(pos) : 0.0;
  }
  return out;
}

// -------------------- am2 --------------------
// [[Rcpp::export(rng = false)]]
arma::rowvec am2(
    const arma::vec& texa,
    const arma::vec& yexa,
    const arma::vec& T0,
    const arma::vec& T,
    const arma::vec& W,
    const arma::mat& X,
    const arma::vec& m
) {
  arma::uword const n  = X.n_rows;
  arma::uword const nm = texa.n_elem;

  arma::vec Lam(n, arma::fill::zeros);
  arma::vec de(nm, arma::fill::zeros);

  arma::uvec const idx = arma::sort_index(texa);

  auto cmp = [](cmp_par const& x, cmp_par const& y){
    return x.first <= y.first;
  };
  std::set<cmp_par, decltype(cmp)> indices(cmp);

  double w_sum{};

  {
    auto const idx_i = idx[0];
    indices.emplace(yexa[idx_i], idx_i);
    w_sum += 1;
    de(idx_i) = w_sum;
  }

  auto indices_head = indices.begin();

  for (arma::uword i = 1; i < nm; ++i) {
    auto const idx_i = idx[i];
    indices.emplace(yexa[idx_i], idx_i);

    if (yexa[idx_i] > indices_head->first)
      while (indices_head->first < texa[idx_i]) {
        w_sum -= 1;
        ++indices_head;
      }
    else
      --indices_head;

    w_sum += 1;
    de[idx_i] = w_sum;

    for (arma::uword j = 1; j <= i; ++j) {
      if (T[idx_i] == T(idx[i - j]))
        de(idx[i - j]) = w_sum;
      else break;
    }
  }

  // --- Optimised suffix-cumsum approach: O((n+nm) log nm) ---
  arma::vec inv_de(nm);
  for (arma::uword i = 0; i < nm; i++)
    inv_de(i) = (de(i) > 0.0) ? 1.0 / de(i) : 0.0;

  // idx is already sort_index(texa) ascending
  arma::vec texa_sorted = texa(idx);
  arma::vec inv_de_sorted = inv_de(idx);
  arma::vec suf_cum = arma::flipud(arma::cumsum(arma::flipud(inv_de_sorted)));

  for (arma::uword k = 0; k < n; ++k) {
    auto it = std::lower_bound(texa_sorted.begin(), texa_sorted.end(), T0[k]);
    arma::uword pos = (arma::uword)std::distance(texa_sorted.begin(), it);
    Lam[k] = (pos < nm) ? suf_cum(pos) : 0.0;
  }

  Lam = exp(-Lam);

  arma::vec R = m / Lam;

  return ((W % R - mean(W % R))).t() * X / n;
}

// ==================== JFM fast EE functions ====================
//
// jfm_s0s1_cpp  — combined S0(t) and S1(t) for one sub-model in one C++ pass.
//   Replaces four separate R calls per iteration (S0t_death, S1t_death,
//   S0t_recurrent, S1t_recurrent).  The same function handles both death
//   and recurrent sub-models; the caller passes the appropriate arguments.
//
// jfm_s0s1_cf_cpp — cross-fitted version: each subject uses the coefficient
//   vector for its held-out fold (cv_fold[j]).
//
// jfm_score_cpp — vectorised score (gradient) given already-computed S0t/S1t.

// Combined S0(t) + S1(t) for one JFM sub-model (death or recurrent).
// Y: n observed times; wt: n×ne weights; t_events: ne sorted event times;
// idx_mat: n×ne 1-based row indices into Z_pseudo (double matrix — R default).
// coef: p-vec.  Returns list(S0t = ne-vec, S1t = p×ne), divided by n.
// [[Rcpp::export(rng = false)]]
Rcpp::List jfm_s0s1_cpp(
    const arma::vec& Y,
    const arma::mat& wt,
    const arma::vec& t_events,
    const arma::mat& idx_mat,
    const arma::mat& Z_pseudo,
    const arma::vec& coef
) {
  arma::uword n  = Y.n_elem;
  arma::uword ne = t_events.n_elem;
  arma::uword p  = coef.n_elem;

  // Precompute one linear predictor per pseudo-entry row: O(npseudo × p) BLAS
  arma::vec lp_all = Z_pseudo * coef;

  arma::vec S0t(ne, arma::fill::zeros);
  arma::mat S1t(p,  ne, arma::fill::zeros);
  arma::vec s1_buf(p);

  arma::uword nrows = Z_pseudo.n_rows;
  for (arma::uword i = 0; i < ne; i++) {
    double s0 = 0.0;
    s1_buf.zeros();
    double t_i = t_events(i);

    for (arma::uword j = 0; j < n; j++) {
      if (Y(j) >= t_i) {
        arma::uword k = (arma::uword)(idx_mat(j, i) - 1.0);  // 0-based
        double w = wt(j, i) * std::exp(lp_all(k));
        s0 += w;
        // Accumulate w * Z_pseudo[k, :] into s1_buf (column-major access)
        const double* zk = Z_pseudo.memptr() + k;
        for (arma::uword l = 0; l < p; l++)
          s1_buf(l) += w * zk[l * nrows];
      }
    }
    S0t(i)    = s0;
    S1t.col(i) = s1_buf;
  }

  double inv_n = 1.0 / (double)n;
  arma::vec s0t_out = S0t * inv_n;
  return Rcpp::List::create(
    Rcpp::Named("S0t") = Rcpp::NumericVector(s0t_out.begin(), s0t_out.end()),
    Rcpp::Named("S1t") = S1t * inv_n
  );
}

// Cross-fitted S0(t)+S1(t): subject j uses coef_mat.row(cv_fold[j]).
// coef_mat: K×p; cv_fold: n-vec of 0-based fold indices (double for R compat).
// [[Rcpp::export(rng = false)]]
Rcpp::List jfm_s0s1_cf_cpp(
    const arma::vec& Y,
    const arma::mat& wt,
    const arma::vec& t_events,
    const arma::mat& idx_mat,
    const arma::mat& Z_pseudo,
    const arma::mat& coef_mat,
    const arma::vec& cv_fold
) {
  arma::uword n       = Y.n_elem;
  arma::uword ne      = t_events.n_elem;
  arma::uword p       = coef_mat.n_cols;
  arma::uword K       = coef_mat.n_rows;
  arma::uword npseudo = Z_pseudo.n_rows;

  // Precompute lp for each fold: lp_by_fold[k, fold] = Z_pseudo[k,:] . coef_mat[fold,:]
  // Shape: npseudo × K
  arma::mat lp_by_fold(npseudo, K);
  for (arma::uword f = 0; f < K; f++)
    lp_by_fold.col(f) = Z_pseudo * coef_mat.row(f).t();

  arma::vec S0t(ne, arma::fill::zeros);
  arma::mat S1t(p,  ne, arma::fill::zeros);
  arma::vec s1_buf(p);

  arma::uword nrows = Z_pseudo.n_rows;
  for (arma::uword i = 0; i < ne; i++) {
    double s0 = 0.0;
    s1_buf.zeros();
    double t_i = t_events(i);

    for (arma::uword j = 0; j < n; j++) {
      if (Y(j) >= t_i) {
        arma::uword k    = (arma::uword)(idx_mat(j, i) - 1.0);
        arma::uword fold = (arma::uword)cv_fold(j);
        double w = wt(j, i) * std::exp(lp_by_fold(k, fold));
        s0 += w;
        const double* zk = Z_pseudo.memptr() + k;
        for (arma::uword l = 0; l < p; l++)
          s1_buf(l) += w * zk[l * nrows];
      }
    }
    S0t(i)    = s0;
    S1t.col(i) = s1_buf;
  }

  double inv_n = 1.0 / (double)n;
  arma::vec s0t_out = S0t * inv_n;
  return Rcpp::List::create(
    Rcpp::Named("S0t") = Rcpp::NumericVector(s0t_out.begin(), s0t_out.end()),
    Rcpp::Named("S1t") = S1t * inv_n
  );
}

// Score U = sum_i [z_{event_id[i],i} - S1t[:,i]/S0t[i]].
// idx_mat: n×ne 1-based (double); event_id: ne 0-based subject indices (double).
// [[Rcpp::export(rng = false)]]
arma::vec jfm_score_cpp(
    const arma::mat& idx_mat,
    const arma::vec& event_id,
    const arma::mat& Z_pseudo,
    const arma::mat& S1t,
    const arma::vec& S0t
) {
  arma::uword ne = S0t.n_elem;
  arma::uword p  = S1t.n_rows;
  arma::vec U(p, arma::fill::zeros);

  arma::uword nrows = Z_pseudo.n_rows;
  for (arma::uword i = 0; i < ne; i++) {
    arma::uword subj = (arma::uword)event_id(i);
    arma::uword k    = (arma::uword)(idx_mat(subj, i) - 1.0);
    double inv_s0    = 1.0 / S0t(i);
    const double* zk = Z_pseudo.memptr() + k;
    const double* s1 = S1t.colptr(i);
    double* u_ptr    = U.memptr();
    for (arma::uword l = 0; l < p; l++)
      u_ptr[l] += zk[l * nrows] - s1[l] * inv_s0;
  }
  return U;
}
