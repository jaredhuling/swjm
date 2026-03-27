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

// ==================== Fast O(n_pseudo) JFM EE functions ====================
//
// These replace the O(n × ne) jfm_s0s1_cpp with a timeline-based incremental
// scan.  The timeline is precomputed once (data-dependent, not
// coefficient-dependent) and reused at every stagewise iteration.

// Precompute sorted timeline of enter/exit/query events.
//
// entry_times: n_pseudo entry (left-truncation) times
// exit_times:  n_pseudo exit times (next entry start or Y_i)
// is_last:     n_pseudo int (1 if last entry for its subject, 0 otherwise)
// query_times: ne sorted event times where S0t/S1t are needed
// query_before_enter: int, if 1 then QUERY has priority over ENTER at same
//     time (recurrent sub-model: uses L < t, strict).  If 0, ENTER has
//     priority over QUERY (death sub-model: uses L <= t, non-strict).
//
// Returns list(type = int vec, idx = int vec, size = int).
// Output type values:  1 = ENTER, 2 = QUERY, 0 or 3 = EXIT.
//
// [[Rcpp::export(rng = false)]]
Rcpp::List jfm_precompute_timeline(
    const arma::vec& entry_times,
    const arma::vec& exit_times,
    const arma::ivec& is_last,
    const arma::vec& query_times,
    int query_before_enter = 0
) {
  arma::uword n_pseudo = entry_times.n_elem;
  arma::uword ne = query_times.n_elem;
  arma::uword total = 2 * n_pseudo + ne;

  // Build unsorted timeline: (time, sort_priority, operation_type, index)
  // sort_priority determines order at ties; operation_type is used in the scan.
  std::vector<double> tl_time(total);
  std::vector<int>    tl_prio(total);   // for sorting only
  std::vector<int>    tl_type(total);   // operation type for scan
  std::vector<int>    tl_idx(total);

  // Priority assignments depend on query_before_enter:
  //   death (qbe=0): EXIT_nonlast(0) < ENTER(1) < QUERY(2) < EXIT_last(3)
  //   recur (qbe=1): QUERY(0)  < EXIT_nonlast(1) < ENTER(2) < EXIT_last(3)
  int prio_enter, prio_query, prio_exit_nonlast, prio_exit_last;
  if (query_before_enter) {
    prio_query        = 0;
    prio_exit_nonlast = 1;
    prio_enter        = 2;
    prio_exit_last    = 3;
  } else {
    prio_exit_nonlast = 0;
    prio_enter        = 1;
    prio_query        = 2;
    prio_exit_last    = 3;
  }

  arma::uword pos = 0;

  // ENTER events
  for (arma::uword j = 0; j < n_pseudo; j++) {
    tl_time[pos] = entry_times(j);
    tl_prio[pos] = prio_enter;
    tl_type[pos] = 1;   // ENTER
    tl_idx[pos]  = (int)j;
    pos++;
  }

  // EXIT events
  for (arma::uword j = 0; j < n_pseudo; j++) {
    tl_time[pos] = exit_times(j);
    tl_prio[pos] = is_last(j) ? prio_exit_last : prio_exit_nonlast;
    tl_type[pos] = is_last(j) ? 3 : 0;
    tl_idx[pos]  = (int)j;
    pos++;
  }

  // QUERY events
  for (arma::uword k = 0; k < ne; k++) {
    tl_time[pos] = query_times(k);
    tl_prio[pos] = prio_query;
    tl_type[pos] = 2;   // QUERY
    tl_idx[pos]  = (int)k;
    pos++;
  }

  // Stable sort by (time, priority)
  std::vector<arma::uword> order(total);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(),
    [&](arma::uword a, arma::uword b) {
      if (tl_time[a] != tl_time[b]) return tl_time[a] < tl_time[b];
      return tl_prio[a] < tl_prio[b];
    });

  // Write sorted output (operation type, not priority)
  Rcpp::IntegerVector out_type(total);
  Rcpp::IntegerVector out_idx(total);
  for (arma::uword i = 0; i < total; i++) {
    out_type[i] = tl_type[order[i]];
    out_idx[i]  = tl_idx[order[i]];
  }

  return Rcpp::List::create(
    Rcpp::Named("type") = out_type,
    Rcpp::Named("idx")  = out_idx,
    Rcpp::Named("size") = (int)total
  );
}


// Build timeline with weight-update events at specified times.
// Same as jfm_precompute_timeline but adds type-4 WEIGHT_UPDATE events.
//
// wt_update_times: n_wt sorted times where weights change (typically death times)
// query_before_enter: same as jfm_precompute_timeline
//
// WEIGHT_UPDATE priority: processed just before QUERYs (so queries see updated weights)
//
// [[Rcpp::export(rng = false)]]
Rcpp::List jfm_precompute_timeline_wt(
    const arma::vec& entry_times,
    const arma::vec& exit_times,
    const arma::ivec& is_last,
    const arma::vec& query_times,
    const arma::vec& wt_update_times,
    int query_before_enter = 0
) {
  arma::uword n_pseudo = entry_times.n_elem;
  arma::uword ne = query_times.n_elem;
  arma::uword n_wt = wt_update_times.n_elem;
  arma::uword total = 2 * n_pseudo + ne + n_wt;

  std::vector<double> tl_time(total);
  std::vector<int>    tl_prio(total);
  std::vector<int>    tl_type(total);
  std::vector<int>    tl_idx(total);

  // Priority assignments:
  //   death (qbe=0): EXIT_nonlast(0) < ENTER(1) < WT_UPDATE(2) < QUERY(3) < EXIT_last(4)
  //   recur (qbe=1): WT_UPDATE(0) < QUERY(1) < EXIT_nonlast(2) < ENTER(3) < EXIT_last(4)
  int prio_enter, prio_query, prio_exit_nonlast, prio_exit_last, prio_wt;
  if (query_before_enter) {
    prio_wt           = 0;
    prio_query        = 1;
    prio_exit_nonlast = 2;
    prio_enter        = 3;
    prio_exit_last    = 4;
  } else {
    prio_exit_nonlast = 0;
    prio_enter        = 1;
    prio_wt           = 2;
    prio_query        = 3;
    prio_exit_last    = 4;
  }

  arma::uword pos = 0;

  for (arma::uword j = 0; j < n_pseudo; j++) {
    tl_time[pos] = entry_times(j);
    tl_prio[pos] = prio_enter;
    tl_type[pos] = 1;
    tl_idx[pos]  = (int)j;
    pos++;
  }
  for (arma::uword j = 0; j < n_pseudo; j++) {
    tl_time[pos] = exit_times(j);
    tl_prio[pos] = is_last(j) ? prio_exit_last : prio_exit_nonlast;
    tl_type[pos] = is_last(j) ? 3 : 0;
    tl_idx[pos]  = (int)j;
    pos++;
  }
  for (arma::uword k = 0; k < ne; k++) {
    tl_time[pos] = query_times(k);
    tl_prio[pos] = prio_query;
    tl_type[pos] = 2;
    tl_idx[pos]  = (int)k;
    pos++;
  }
  for (arma::uword k = 0; k < n_wt; k++) {
    tl_time[pos] = wt_update_times(k);
    tl_prio[pos] = prio_wt;
    tl_type[pos] = 4;
    tl_idx[pos]  = (int)k;
    pos++;
  }

  std::vector<arma::uword> order(total);
  std::iota(order.begin(), order.end(), 0);
  std::stable_sort(order.begin(), order.end(),
    [&](arma::uword a, arma::uword b) {
      if (tl_time[a] != tl_time[b]) return tl_time[a] < tl_time[b];
      return tl_prio[a] < tl_prio[b];
    });

  Rcpp::IntegerVector out_type(total);
  Rcpp::IntegerVector out_idx(total);
  for (arma::uword i = 0; i < total; i++) {
    out_type[i] = tl_type[order[i]];
    out_idx[i]  = tl_idx[order[i]];
  }

  return Rcpp::List::create(
    Rcpp::Named("type") = out_type,
    Rcpp::Named("idx")  = out_idx,
    Rcpp::Named("size") = (int)total
  );
}


// Fast O(n_pseudo) S0t/S1t via precomputed timeline.
//
// tl_type, tl_idx: sorted timeline from jfm_precompute_timeline
// tl_size: length of timeline
// Z_pseudo: n_pseudo × p covariate matrix
// coef: p-vector
// ne: number of query events (output size)
// n: number of subjects (for /n normalisation)
//
// Returns list(S0t = ne-vec, S1t = p × ne matrix), divided by n.
//
// [[Rcpp::export(rng = false)]]
Rcpp::List jfm_s0s1_fast_cpp(
    const Rcpp::IntegerVector& tl_type,
    const Rcpp::IntegerVector& tl_idx,
    int tl_size,
    const arma::mat& Z_pseudo,
    const arma::vec& coef,
    int ne,
    int n
) {
  arma::uword p     = coef.n_elem;
  arma::uword nrows = Z_pseudo.n_rows;

  // Precompute exp(z' coef) for every pseudo-entry: O(n_pseudo × p)
  arma::vec exp_lp = arma::exp(Z_pseudo * coef);

  // Running accumulators
  double run_s0 = 0.0;
  arma::vec run_s1(p, arma::fill::zeros);

  // Output
  arma::vec S0t(ne, arma::fill::zeros);
  arma::mat S1t(p, ne, arma::fill::zeros);

  for (int i = 0; i < tl_size; i++) {
    int type = tl_type[i];
    int idx  = tl_idx[i];

    if (type == 1) {
      // ENTER: add pseudo-entry contribution
      double e = exp_lp(idx);
      run_s0 += e;
      const double* zk = Z_pseudo.memptr() + idx;
      for (arma::uword l = 0; l < p; l++)
        run_s1(l) += e * zk[l * nrows];

    } else if (type == 0 || type == 3) {
      // EXIT: remove pseudo-entry contribution
      double e = exp_lp(idx);
      run_s0 -= e;
      const double* zk = Z_pseudo.memptr() + idx;
      for (arma::uword l = 0; l < p; l++)
        run_s1(l) -= e * zk[l * nrows];

    } else {
      // QUERY: record S0t and S1t
      S0t(idx)     = run_s0;
      S1t.col(idx) = run_s1;
    }
  }

  double inv_n = 1.0 / (double)n;
  arma::vec s0t_out = S0t * inv_n;
  return Rcpp::List::create(
    Rcpp::Named("S0t") = Rcpp::NumericVector(s0t_out.begin(), s0t_out.end()),
    Rcpp::Named("S1t") = S1t * inv_n
  );
}


// Cross-fitted fast S0t/S1t: each pseudo-entry uses the coefficient from
// its subject's held-out fold.
//
// subject_of_entry: n_pseudo int vector, 0-based subject index per entry
// cv_fold:          n int vector, 0-based fold index per subject
// coef_mat:         K × p matrix of fold-specific coefficients
//
// [[Rcpp::export(rng = false)]]
Rcpp::List jfm_s0s1_fast_cf_cpp(
    const Rcpp::IntegerVector& tl_type,
    const Rcpp::IntegerVector& tl_idx,
    int tl_size,
    const arma::mat& Z_pseudo,
    const arma::mat& coef_mat,
    const Rcpp::IntegerVector& subject_of_entry,
    const Rcpp::IntegerVector& cv_fold,
    int ne,
    int n
) {
  arma::uword p      = coef_mat.n_cols;
  arma::uword nrows  = Z_pseudo.n_rows;
  arma::uword npseudo = nrows;

  // Precompute exp(z_j' · coef_{fold(subject(j))}) for each pseudo-entry
  arma::vec exp_lp(npseudo);
  for (arma::uword j = 0; j < npseudo; j++) {
    int subj = subject_of_entry[j];
    int fold = cv_fold[subj];
    double lp = arma::dot(Z_pseudo.row(j), coef_mat.row(fold));
    exp_lp(j) = std::exp(lp);
  }

  // Running accumulators and scan (same structure as non-CF version)
  double run_s0 = 0.0;
  arma::vec run_s1(p, arma::fill::zeros);

  arma::vec S0t(ne, arma::fill::zeros);
  arma::mat S1t(p, ne, arma::fill::zeros);

  for (int i = 0; i < tl_size; i++) {
    int type = tl_type[i];
    int idx  = tl_idx[i];

    if (type == 1) {
      double e = exp_lp(idx);
      run_s0 += e;
      const double* zk = Z_pseudo.memptr() + idx;
      for (arma::uword l = 0; l < p; l++)
        run_s1(l) += e * zk[l * nrows];

    } else if (type == 0 || type == 3) {
      double e = exp_lp(idx);
      run_s0 -= e;
      const double* zk = Z_pseudo.memptr() + idx;
      for (arma::uword l = 0; l < p; l++)
        run_s1(l) -= e * zk[l * nrows];

    } else {
      S0t(idx)     = run_s0;
      S1t.col(idx) = run_s1;
    }
  }

  double inv_n = 1.0 / (double)n;
  arma::vec s0t_out = S0t * inv_n;
  return Rcpp::List::create(
    Rcpp::Named("S0t") = Rcpp::NumericVector(s0t_out.begin(), s0t_out.end()),
    Rcpp::Named("S1t") = S1t * inv_n
  );
}


// Weighted O(n_pseudo + n*n_D) S0t/S1t via precomputed timeline.
//
// Same as jfm_s0s1_fast_cpp but with per-subject, time-varying weights.
// The timeline must include WEIGHT_UPDATE events (type 4) at each death
// time where weights change.  At these events the function recomputes all
// active entries' contributions using the new weights from wt_subj.
//
// tl_type, tl_idx: sorted timeline (types: 0/3=EXIT, 1=ENTER, 2=QUERY, 4=WEIGHT_UPDATE)
// tl_size: length of timeline
// Z_pseudo: n_pseudo × p covariate matrix
// entry_subject: n_pseudo int, 0-based subject index for each pseudo-entry
// coef: p-vector
// wt_subj: n × n_wt weight matrix (subjects × weight-update times, column-major)
//          At weight-update event with index k, subject i's new weight is wt_subj(i, k).
// ne: number of QUERY events (output size)
// n: number of subjects
//
// Returns list(S0t = ne-vec, S1t = p × ne matrix), divided by n.
//
// [[Rcpp::export(rng = false)]]
Rcpp::List jfm_s0s1_wt_fast_cpp(
    const Rcpp::IntegerVector& tl_type,
    const Rcpp::IntegerVector& tl_idx,
    int tl_size,
    const arma::mat& Z_pseudo,
    const Rcpp::IntegerVector& entry_subject,
    const arma::vec& coef,
    const arma::mat& wt_subj,
    int ne,
    int n
) {
  arma::uword p     = coef.n_elem;
  arma::uword nrows = Z_pseudo.n_rows;
  arma::uword npseudo = nrows;
  arma::uword nsub  = (arma::uword)n;

  // Precompute exp(z' coef) for every pseudo-entry
  arma::vec exp_lp = arma::exp(Z_pseudo * coef);

  // Per-entry current weight (initialised to 1; updated at WEIGHT_UPDATE events)
  arma::vec entry_wt(npseudo, arma::fill::ones);

  // Track which pseudo-entry is currently active for each subject (-1 = none)
  std::vector<int> active_entry(nsub, -1);

  // Running accumulators
  double run_s0 = 0.0;
  arma::vec run_s1(p, arma::fill::zeros);

  // Output
  arma::vec S0t(ne, arma::fill::zeros);
  arma::mat S1t(p, ne, arma::fill::zeros);

  for (int i = 0; i < tl_size; i++) {
    int type = tl_type[i];
    int idx  = tl_idx[i];

    if (type == 1) {
      // ENTER: add pseudo-entry contribution with current weight
      int subj = entry_subject[idx];
      entry_wt(idx) = (active_entry[subj] >= 0) ?
                       entry_wt(active_entry[subj]) : 1.0;
      // If subject already had an active entry, the EXIT for that entry
      // should have been processed already (priority ordering ensures this)
      active_entry[subj] = idx;
      double e = entry_wt(idx) * exp_lp(idx);
      run_s0 += e;
      const double* zk = Z_pseudo.memptr() + idx;
      for (arma::uword l = 0; l < p; l++)
        run_s1(l) += e * zk[l * nrows];

    } else if (type == 0 || type == 3) {
      // EXIT: remove pseudo-entry contribution
      double e = entry_wt(idx) * exp_lp(idx);
      run_s0 -= e;
      const double* zk = Z_pseudo.memptr() + idx;
      for (arma::uword l = 0; l < p; l++)
        run_s1(l) -= e * zk[l * nrows];
      int subj = entry_subject[idx];
      if (active_entry[subj] == idx) active_entry[subj] = -1;

    } else if (type == 4) {
      // WEIGHT_UPDATE: idx is the column index into wt_subj.
      // Recompute all active entries' contributions with new weights.
      // Reset running sums and recompute from active entries.
      run_s0 = 0.0;
      run_s1.zeros();
      for (arma::uword s = 0; s < nsub; s++) {
        int ae = active_entry[s];
        if (ae < 0) continue;
        double new_w = wt_subj(s, idx);
        entry_wt(ae) = new_w;
        double e = new_w * exp_lp(ae);
        run_s0 += e;
        const double* zk = Z_pseudo.memptr() + ae;
        for (arma::uword l = 0; l < p; l++)
          run_s1(l) += e * zk[l * nrows];
      }

    } else {
      // QUERY: record S0t and S1t
      S0t(idx)     = run_s0;
      S1t.col(idx) = run_s1;
    }
  }

  double inv_n = 1.0 / (double)n;
  arma::vec s0t_out = S0t * inv_n;
  return Rcpp::List::create(
    Rcpp::Named("S0t") = Rcpp::NumericVector(s0t_out.begin(), s0t_out.end()),
    Rcpp::Named("S1t") = S1t * inv_n
  );
}


// Simplified score using precomputed event pseudo-entry indices.
//
// event_pseudo_idx: ne-vector of 0-based pseudo-entry row indices
// Z_pseudo: n_pseudo × p covariate matrix
// S1t: p × ne matrix
// S0t: ne vector
//
// Returns p-vector: U = sum_k [z_{event_pseudo_idx[k]} - S1t[:,k]/S0t[k]]
//
// Batched cross-fitted CV score norms: evaluates score norms at ALL lambda
// points in one C++ call.  Avoids the R loop overhead that dominates CV time.
//
// tl_type, tl_idx, tl_size: precomputed timeline
// Z_pseudo: n_pseudo × p
// coef_path: p × n_lambda coefficient matrix (columns = lambda points)
// subject_of_entry: n_pseudo int, 0-based subject per entry
// cv_fold: n int, 0-based fold per subject
// event_pseudo_idx: ne int, 0-based entry per event
// ne, n: event count and subject count
//
// Returns: n_lambda-vector of score L2 norms
//
// [[Rcpp::export(rng = false)]]
arma::vec jfm_cv_scores_batch_cpp(
    const Rcpp::IntegerVector& tl_type,
    const Rcpp::IntegerVector& tl_idx,
    int tl_size,
    const arma::mat& Z_pseudo,
    const arma::mat& coef_path,
    const Rcpp::IntegerVector& subject_of_entry,
    const Rcpp::IntegerVector& cv_fold,
    const Rcpp::IntegerVector& event_pseudo_idx,
    int ne,
    int n
) {
  arma::uword p       = Z_pseudo.n_cols;
  arma::uword nrows   = Z_pseudo.n_rows;
  arma::uword n_lam   = coef_path.n_cols;
  arma::uword K       = 0;
  for (int i = 0; i < (int)n; i++)
    if ((arma::uword)cv_fold[i] >= K) K = cv_fold[i] + 1;

  double inv_n = 1.0 / (double)n;
  arma::vec score_norms(n_lam);

  // For each lambda point
  for (arma::uword lam = 0; lam < n_lam; lam++) {
    arma::vec coef = coef_path.col(lam);

    // Build K fold-specific coef vectors (each fold uses a different coef)
    // In cross-fitting, fold k's subjects use the coef trained WITHOUT fold k.
    // coef_path already contains the cross-fitted coef for this lambda.
    // Actually, coef_path.col(lam) is a SINGLE coef vector, but for CF
    // we need K different coefs.  The caller should pass a K×p matrix per lambda.
    // For simplicity, this function handles the common case where all folds
    // use the same coef (non-CF), and we provide a separate batch CF function.

    // Precompute exp(z' coef) for all pseudo-entries
    arma::vec exp_lp = arma::exp(Z_pseudo * coef);

    // Timeline scan
    double run_s0 = 0.0;
    arma::vec run_s1(p, arma::fill::zeros);
    arma::vec S0t(ne); arma::mat S1t(p, ne);

    for (int i = 0; i < tl_size; i++) {
      int type = tl_type[i];
      int idx  = tl_idx[i];
      if (type == 1) {
        double e = exp_lp(idx);
        run_s0 += e;
        const double* zk = Z_pseudo.memptr() + idx;
        for (arma::uword l = 0; l < p; l++) run_s1(l) += e * zk[l * nrows];
      } else if (type == 0 || type == 3) {
        double e = exp_lp(idx);
        run_s0 -= e;
        const double* zk = Z_pseudo.memptr() + idx;
        for (arma::uword l = 0; l < p; l++) run_s1(l) -= e * zk[l * nrows];
      } else {
        S0t(idx)     = run_s0 * inv_n;
        S1t.col(idx) = run_s1 * inv_n;
      }
    }

    // Compute score
    arma::vec U(p, arma::fill::zeros);
    for (arma::uword i = 0; i < (arma::uword)ne; i++) {
      arma::uword k = (arma::uword)event_pseudo_idx[i];
      double inv_s0 = 1.0 / S0t(i);
      const double* zk = Z_pseudo.memptr() + k;
      const double* s1 = S1t.colptr(i);
      for (arma::uword l = 0; l < p; l++)
        U(l) += zk[l * nrows] - s1[l] * inv_s0;
    }
    U *= inv_n;
    score_norms(lam) = std::sqrt(arma::dot(U, U));
  }

  return score_norms;
}


// Batched cross-fitted version: coef_array is K × p × n_lambda (flattened as
// K*p × n_lambda matrix, where rows 0..p-1 = fold 0, rows p..2p-1 = fold 1, etc.)
//
// [[Rcpp::export(rng = false)]]
arma::vec jfm_cv_scores_batch_cf_cpp(
    const Rcpp::IntegerVector& tl_type,
    const Rcpp::IntegerVector& tl_idx,
    int tl_size,
    const arma::mat& Z_pseudo,
    const arma::mat& coef_array,
    int K_folds,
    const Rcpp::IntegerVector& subject_of_entry,
    const Rcpp::IntegerVector& cv_fold,
    const Rcpp::IntegerVector& event_pseudo_idx,
    int ne,
    int n
) {
  arma::uword p       = Z_pseudo.n_cols;
  arma::uword nrows   = Z_pseudo.n_rows;
  arma::uword npseudo = nrows;
  arma::uword n_lam   = coef_array.n_cols;
  arma::uword K       = (arma::uword)K_folds;

  double inv_n = 1.0 / (double)n;
  arma::vec score_norms(n_lam);

  for (arma::uword lam = 0; lam < n_lam; lam++) {
    // Extract K fold-specific coefs from flattened array
    // coef_array rows: [fold0_coef(p), fold1_coef(p), ..., foldK-1_coef(p)]
    arma::mat coef_by_fold(K, p);
    for (arma::uword f = 0; f < K; f++)
      coef_by_fold.row(f) = coef_array.col(lam).subvec(f * p, (f + 1) * p - 1).t();

    // Precompute exp(z_j' coef_{fold(subject(j))})
    arma::vec exp_lp(npseudo);
    for (arma::uword j = 0; j < npseudo; j++) {
      int fold = cv_fold[subject_of_entry[j]];
      exp_lp(j) = std::exp(arma::dot(Z_pseudo.row(j), coef_by_fold.row(fold)));
    }

    // Timeline scan
    double run_s0 = 0.0;
    arma::vec run_s1(p, arma::fill::zeros);
    arma::vec S0t(ne); arma::mat S1t(p, ne);

    for (int i = 0; i < tl_size; i++) {
      int type = tl_type[i];
      int idx  = tl_idx[i];
      if (type == 1) {
        double e = exp_lp(idx);
        run_s0 += e;
        const double* zk = Z_pseudo.memptr() + idx;
        for (arma::uword l = 0; l < p; l++) run_s1(l) += e * zk[l * nrows];
      } else if (type == 0 || type == 3) {
        double e = exp_lp(idx);
        run_s0 -= e;
        const double* zk = Z_pseudo.memptr() + idx;
        for (arma::uword l = 0; l < p; l++) run_s1(l) -= e * zk[l * nrows];
      } else {
        S0t(idx)     = run_s0 * inv_n;
        S1t.col(idx) = run_s1 * inv_n;
      }
    }

    // Score
    arma::vec U(p, arma::fill::zeros);
    for (arma::uword i = 0; i < (arma::uword)ne; i++) {
      arma::uword k = (arma::uword)event_pseudo_idx[i];
      double inv_s0 = 1.0 / S0t(i);
      const double* zk = Z_pseudo.memptr() + k;
      const double* s1 = S1t.colptr(i);
      for (arma::uword l = 0; l < p; l++)
        U(l) += zk[l * nrows] - s1[l] * inv_s0;
    }
    U *= inv_n;
    score_norms(lam) = std::sqrt(arma::dot(U, U));
  }

  return score_norms;
}


// Simplified score using precomputed event pseudo-entry indices.
// [[Rcpp::export(rng = false)]]
arma::vec jfm_score_fast_cpp(
    const Rcpp::IntegerVector& event_pseudo_idx,
    const arma::mat& Z_pseudo,
    const arma::mat& S1t,
    const arma::vec& S0t
) {
  arma::uword ne = S0t.n_elem;
  arma::uword p  = S1t.n_rows;
  arma::vec U(p, arma::fill::zeros);

  arma::uword nrows = Z_pseudo.n_rows;
  for (arma::uword i = 0; i < ne; i++) {
    arma::uword k = (arma::uword)event_pseudo_idx[i];
    double inv_s0 = 1.0 / S0t(i);
    const double* zk = Z_pseudo.memptr() + k;
    const double* s1 = S1t.colptr(i);
    double* u_ptr    = U.memptr();
    for (arma::uword l = 0; l < p; l++)
      u_ptr[l] += zk[l * nrows] - s1[l] * inv_s0;
  }
  return U;
}


// ==================== Fast JFM preprocessing ====================
//
// Builds sorted pseudo-data, exit times, is_last, subject mapping,
// index matrices, weight matrix, and event pseudo indices in one
// efficient O(n_pseudo + n*n_D + n*n_R) pass.
//
// [[Rcpp::export(rng = false)]]
Rcpp::List jfm_build_pseudo_cpp(
    const arma::vec& t_start,
    const arma::vec& subj_id_raw,
    const arma::mat& Z_raw,
    const arma::vec& Y,
    const arma::vec& td,
    const arma::vec& td_id_raw,
    const arma::vec& tr,
    const arma::vec& tr_id_raw,
    const arma::vec& coef_death,
    const arma::vec& lambda0_d,
    double theta
) {
  arma::uword n_pseudo = t_start.n_elem;
  arma::uword n        = Y.n_elem;
  arma::uword n_D      = td.n_elem;
  arma::uword n_R      = tr.n_elem;
  arma::uword p        = Z_raw.n_cols;

  // 1. Sort pseudo-data by t_start
  arma::uvec sort_idx = arma::sort_index(t_start);
  arma::vec  L_sorted = t_start(sort_idx);
  arma::vec  I_sorted(n_pseudo);
  arma::mat  Z_sorted(n_pseudo, p);
  for (arma::uword i = 0; i < n_pseudo; i++) {
    I_sorted(i) = subj_id_raw(sort_idx(i));
    Z_sorted.row(i) = Z_raw.row(sort_idx(i));
  }

  // 2. Compute exit times + is_last via per-subject entry lists
  arma::vec exit_times(n_pseudo);
  arma::ivec is_last_vec(n_pseudo, arma::fill::zeros);
  std::vector<std::vector<arma::uword>> subj_rows(n);
  for (arma::uword i = 0; i < n_pseudo; i++)
    subj_rows[(int)I_sorted(i) - 1].push_back(i);

  for (arma::uword s = 0; s < n; s++) {
    auto& rows = subj_rows[s];
    arma::uword nr = rows.size();
    for (arma::uword j = 0; j + 1 < nr; j++)
      exit_times(rows[j]) = L_sorted(rows[j + 1]);
    exit_times(rows[nr - 1]) = Y(s);
    is_last_vec(rows[nr - 1]) = 1;
  }

  // 3. Sort death times
  arma::uvec td_order = arma::sort_index(td);
  arma::vec td_sorted = td(td_order);
  arma::vec td_id_sorted(n_D);
  for (arma::uword k = 0; k < n_D; k++)
    td_id_sorted(k) = td_id_raw(td_order(k));

  // 4. Sort recurrent times
  arma::uvec tr_order = arma::sort_index(tr);
  arma::vec tr_sorted = tr(tr_order);
  arma::vec tr_id_sorted(n_R);
  for (arma::uword k = 0; k < n_R; k++)
    tr_id_sorted(k) = tr_id_raw(tr_order(k));

  // 5. Build death index matrix + weight matrix: single forward scan
  arma::mat idx_de(n, n_D);
  arma::mat wt_de(n, n_D);
  arma::vec cum_integral(n, arma::fill::zeros);
  std::vector<int> active(n, -1);
  arma::uword ptr = 0;
  arma::vec lp_death = Z_sorted * coef_death;

  for (arma::uword k = 0; k < n_D; k++) {
    double t_dk = td_sorted(k);
    while (ptr < n_pseudo && L_sorted(ptr) <= t_dk) {
      active[(int)I_sorted(ptr) - 1] = (int)ptr;
      ptr++;
    }
    for (arma::uword s = 0; s < n; s++) {
      idx_de(s, k) = (double)(active[s] + 1);
      if (active[s] >= 0)
        cum_integral(s) += lambda0_d(k) * std::exp(lp_death(active[s]));
      wt_de(s, k) = 1.0 / (1.0 + theta * cum_integral(s));
    }
  }

  // 6. Build recurrent index matrix: L < t (strict)
  arma::mat idx_re(n, n_R);
  std::vector<int> active_re(n, -1);
  arma::uword ptr_re = 0;

  for (arma::uword k = 0; k < n_R; k++) {
    double t_rk = tr_sorted(k);
    while (ptr_re < n_pseudo && L_sorted(ptr_re) < t_rk) {
      active_re[(int)I_sorted(ptr_re) - 1] = (int)ptr_re;
      ptr_re++;
    }
    for (arma::uword s = 0; s < n; s++)
      idx_re(s, k) = (double)(active_re[s] + 1);
  }

  // 7. Event pseudo-entry indices (0-based for C++ score)
  Rcpp::IntegerVector de_epi(n_D), re_epi(n_R);
  for (arma::uword k = 0; k < n_D; k++)
    de_epi[k] = (int)idx_de((int)td_id_sorted(k) - 1, k) - 1;
  for (arma::uword k = 0; k < n_R; k++)
    re_epi[k] = (int)idx_re((int)tr_id_sorted(k) - 1, k) - 1;

  // 8. Entry subject mapping (0-based)
  Rcpp::IntegerVector entry_subject(n_pseudo);
  for (arma::uword i = 0; i < n_pseudo; i++)
    entry_subject[i] = (int)I_sorted(i) - 1;

  return Rcpp::List::create(
    Rcpp::Named("Z_pseudo") = Z_sorted,
    Rcpp::Named("entry_times") = L_sorted,
    Rcpp::Named("exit_times") = exit_times,
    Rcpp::Named("is_last") = is_last_vec,
    Rcpp::Named("entry_subject") = entry_subject,
    Rcpp::Named("idx_de") = idx_de,
    Rcpp::Named("idx_re") = idx_re,
    Rcpp::Named("wt_de") = wt_de,
    Rcpp::Named("td_id") = Rcpp::NumericVector(td_id_sorted.begin(), td_id_sorted.end()),
    Rcpp::Named("tr_id") = Rcpp::NumericVector(tr_id_sorted.begin(), tr_id_sorted.end()),
    Rcpp::Named("td_sorted") = Rcpp::NumericVector(td_sorted.begin(), td_sorted.end()),
    Rcpp::Named("tr_sorted") = Rcpp::NumericVector(tr_sorted.begin(), tr_sorted.end()),
    Rcpp::Named("de_epi") = de_epi,
    Rcpp::Named("re_epi") = re_epi,
    Rcpp::Named("n") = (int)n,
    Rcpp::Named("n_D") = (int)n_D,
    Rcpp::Named("n_R") = (int)n_R
  );
}
