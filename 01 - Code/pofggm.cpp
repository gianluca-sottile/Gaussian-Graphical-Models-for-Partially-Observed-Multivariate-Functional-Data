#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

using namespace arma;
using namespace Rcpp;

// [[Rcpp::export]]
arma::uvec generate_jj(int p, int d) {
  arma::uvec jj(p * d);
  
  for (int i = 0; i < p; ++i) {
    jj.subvec(i * d, (i + 1) * d - 1).fill(i);  // R è 1-based
  }
  
  return jj;
}

// [[Rcpp::export]]
Rcpp::List compute_S_and_rho_max(const arma::cube& Xi) {
  int K = Xi.n_slices;
  int n = Xi.n_rows;
  int p = Xi.n_cols;
  
  arma::cube S(p, p, K, fill::zeros);
  arma::mat S_sum_sq(p, p, fill::zeros);
  for (int k = 0; k < K; ++k) {
    arma::mat Xk = Xi.slice(k);
    S.slice(k) = arma::cov(Xk) * ((n - 1.0) / n);  // varianza scalata
    // Somma dei quadrati delle matrici 
    S_sum_sq += arma::square(S.slice(k));
  }
  
  // Imposta la diagonale a zero
  S_sum_sq.diag().zeros();
  
  // Calcolo rho.max
  double rho_max = std::sqrt(S_sum_sq.max()) / K;
  
  // return rho_max;
  return List::create(
    Named("S") = S,
    Named("rho.max") = rho_max
  );
}

// [[Rcpp::export]]
arma::cube compute_Ximputed(const arma::mat& Phi,
                            const arma::cube& Xi) {
  int n = Xi.n_rows;
  int p = Xi.n_cols;
  // int K = Xi.n_slices;
  int d = Phi.n_rows;
  
  arma::cube Ximputed(d, n, p, arma::fill::zeros);
  
  for (int j = 0; j < p; ++j) {
    arma::mat Xi_j = Xi.col(j);            // dim: n × K
    Ximputed.slice(j) = Phi * Xi_j.t();    // risultato: (d x K) x (K x n) = d × n
  }
  
  return Ximputed;
}

// [[Rcpp::export]]
double integrate_scalar(const arma::vec& x, const arma::vec& fx) {
  // if (x.n_elem != fx.n_elem || x.n_elem < 2)
    //   Rcpp::stop("x e fx devono avere stessa lunghezza ≥ 2.");
  // if (!arma::is_sorted(x))
    //   Rcpp::stop("x deve essere ordinato in modo crescente.");
  // if (!arma::is_finite(x).all() || !arma::is_finite(fx).all())
    //   Rcpp::stop("x e fx devono contenere solo valori finiti.");
  
  arma::vec dx = x.subvec(1, x.n_elem - 1) - x.subvec(0, x.n_elem - 2);
  arma::vec f_sum = fx.subvec(1, fx.n_elem - 1) + fx.subvec(0, fx.n_elem - 2);
  
  return 0.5 * arma::dot(dx, f_sum);
}

// [[Rcpp::export]]
arma::mat integrate_mat(const arma::mat& A, const arma::mat& B, const arma::vec& x) {
  // if (A.n_rows != x.n_elem || B.n_rows != x.n_elem)
    //   Rcpp::stop("Le righe di A e B devono coincidere con la lunghezza di x.");
  // if (!arma::is_sorted(x))
    //   Rcpp::stop("x deve essere ordinato in modo crescente.");
  // if (!arma::is_finite(A).all() || !arma::is_finite(B).all() || !arma::is_finite(x).all())
    //   Rcpp::stop("A, B e x devono contenere solo valori finiti.");
  
  int n = A.n_cols;
  int K = B.n_cols;
  
  arma::mat C(n, K, arma::fill::zeros);
  
  // Pre-calcolo differenze e fattori trapezoidali
  arma::vec dx = x.subvec(1, x.n_elem - 1) - x.subvec(0, x.n_elem - 2);
  
  for (int i = 0; i < n; ++i) {
    arma::vec Ai = A.col(i);
    
    // Costruisci f_matrix: ogni colonna è Ai % B.col(k)
    arma::mat f_matrix(Ai.n_elem, K);
    for (int k = 0; k < K; ++k)
      f_matrix.col(k) = Ai % B.col(k);
    
    // Applica trapezio vettorizzato a ogni colonna
    arma::mat f_up = f_matrix.rows(1, x.n_elem - 1);
    arma::mat f_low = f_matrix.rows(0, x.n_elem - 2);
    arma::rowvec integrals = 0.5 * dx.t() * (f_up + f_low);
    
    C.row(i) = integrals;
  }
  
  return C;
}

// [[Rcpp::export]]
arma::mat integrate_mat_dx(const arma::mat& A, const arma::mat& B, double dx) {
  
  // if (A.n_rows < 2 || B.n_rows < 2)
    //   Rcpp::stop("A e B devono avere almeno 2 righe.");
  // if (A.n_rows != B.n_rows)
    //   Rcpp::stop("A e B devono avere lo stesso numero di righe.");
  // if (!arma::is_finite(A).all() || !arma::is_finite(B).all())
    //   Rcpp::stop("A e B devono contenere solo valori numerici finiti.");
  // if (dx <= 0)
    //   Rcpp::stop("dx deve essere positivo.");
  
  int n = A.n_cols;
  int K = B.n_cols;
  int d = A.n_rows;
  
  arma::mat C(n, K, arma::fill::zeros);
  
  for (int i = 0; i < n; ++i) {
    arma::vec Ai = A.col(i);
    arma::mat f_matrix(d, K);
    
    for (int k = 0; k < K; ++k)
      f_matrix.col(k) = Ai % B.col(k);
    
    arma::mat f_low = f_matrix.rows(0, d - 2);
    arma::mat f_up  = f_matrix.rows(1, d - 1);
    
    arma::rowvec integrals = 0.5 * dx * arma::sum(f_low + f_up, 0);  // somma su righe
    C.row(i) = integrals;
  }
  
  return C;
}

// [[Rcpp::export]]
arma::cube integrate_cube(const arma::cube& A, const arma::mat& B, const arma::vec& x) {
  // if (A.n_rows != x.n_elem || B.n_rows != x.n_elem)
    //   Rcpp::stop("Le righe di A e B devono coincidere con la lunghezza di x.");
  // if (!arma::is_sorted(x))
    //   Rcpp::stop("x deve essere ordinato.");
  // if (!arma::is_finite(A).all() || !arma::is_finite(B).all() || !arma::is_finite(x).all())
    //   Rcpp::stop("A, B e x devono contenere solo valori finiti.");
  
  const int n = A.n_cols;
  const int p = A.n_slices;
  const int K = B.n_cols;
  
  arma::cube C(n, p, K, arma::fill::zeros);
  
  for (int j = 0; j < p; ++j)
    C.col(j) = integrate_mat(A.slice(j), B, x);
  
  return C;
}

// [[Rcpp::export]]
mat kron_sum(const cube& Sigma, const mat& Phi) {
  int K = Phi.n_cols;
  int d = Sigma.n_rows;
  int m = Phi.n_rows;
  
  mat result(d * m, d * m, fill::zeros);
  
  for (int k = 0; k < K; ++k) {
    mat Sigma_k = Sigma.slice(k);               // d x d
    vec phi_k = Phi.col(k);                     // m
    mat outer_phi = phi_k * phi_k.t();          // m x m
    mat kron_prod = kron(Sigma_k, outer_phi); // (d*m) x (d*m)
    
    result += kron_prod;
  }
  
  return result;
}

// [[Rcpp::export]]
arma::mat compute_Ginv(const arma::vec& lambda, double alpha, const arma::mat& Q) {
  arma::vec scale = 1.0 / arma::sqrt(lambda + alpha);
  arma::mat Q_scaled = Q.each_row() % scale.t();
  return Q_scaled * Q_scaled.t();  // Matrice simmetrica
}

// [[Rcpp::export]]
double gcv_fun(double alpha,
               const vec& lambda,
               const mat& Q,
               const mat& Gom,
               const mat& ZidOi,
               const mat& ZidMi,
               const int n_obs) {
  // double df = sum(lambda / (lambda + alpha));
  double df = sum(lambda / (lambda + std::max(alpha, 1e-12)));
  mat Ginv = compute_Ginv(lambda, alpha, Q);
  mat A = Ginv * Gom;
  mat XhatMiss = ZidOi * A;
  mat diff = ZidMi - XhatMiss;
  double rss = accu(square(diff));
  double penalty = 1.0 - (df / n_obs);
  if (std::abs(penalty) < 1e-8) penalty = (penalty >= 0 ? 1e-8 : -1e-8);
  return rss / (penalty * penalty);
}

// Brent minimizer
// [[Rcpp::export]]
double Brent_fmin(double ax, double bx, double tol,
                  const vec& lambda,
                  const mat& Q,
                  const mat& Gom,
                  const mat& ZidOi,
                  const mat& ZidMi,
                  int n_obs) {
  
  const int ITMAX = 100;
  const double CGOLD = 0.3819660;
  const double ZEPS = 1e-10;
  
  double a = std::min(ax, bx);
  double b = std::max(ax, bx);
  double x = a + CGOLD * (b - a);
  double w = x, v = x;
  double fx = gcv_fun(x, lambda, Q, Gom, ZidOi, ZidMi, n_obs);
  double fw = fx, fv = fx;
  double d = 0.0, e = 0.0;
  
  for (int iter = 0; iter < ITMAX; iter++) {
    double xm = 0.5 * (a + b);
    double tol1 = tol * std::abs(x) + ZEPS;
    double tol2 = 2.0 * tol1;
    
    if (std::abs(x - xm) <= (tol2 - 0.5 * (b - a))) break;
    
    double p = 0.0, q = 0.0, r = 0.0;
    if (std::abs(e) > tol1) {
      r = (x - w) * (fx - fv);
      q = (x - v) * (fx - fw);
      p = (x - v) * q - (x - w) * r;
      q = 2.0 * (q - r);
      if (q > 0.0) p = -p;
      q = std::abs(q);
      double etemp = e;
      e = d;
      if ((std::abs(p) >= std::abs(0.5 * q * etemp)) ||
          (p <= q * (a - x)) || (p >= q * (b - x))) {
        e = (x >= xm ? a - x : b - x);
        d = CGOLD * e;
      } else {
        d = p / q;
        double u = x + d;
        if ((u - a < tol2) || (b - u < tol2))
          d = tol1 * ((xm - x >= 0) ? 1 : -1);
      }
    } else {
      e = (x >= xm ? a - x : b - x);
      d = CGOLD * e;
    }
    
    double u = x + (std::abs(d) >= tol1 ? d : tol1 * ((d > 0) ? 1 : -1));
    double fu = gcv_fun(u, lambda, Q, Gom, ZidOi, ZidMi, n_obs);
    
    if (fu <= fx) {
      if (u >= x) a = x; else b = x;
      v = w; fv = fw;
      w = x; fw = fx;
      x = u; fx = fu;
    } else {
      if (u < x) a = u; else b = u;
      if ((fu <= fw) || (w == x)) {
        v = w; fv = fw;
        w = u; fw = fu;
      } else if ((fu <= fv) || (v == x) || (v == w)) {
        v = u; fv = fu;
      }
    }
  }
  
  return x;
}

// Robust Brent minimizer with grid presearch
// [[Rcpp::export]]
double Brent_fmin_robust(double ax, double bx, double tol,
                         const arma::vec& lambda,
                         const arma::mat& Q,
                         const arma::mat& Gom,
                         const arma::mat& ZidOi,
                         const arma::mat& ZidMi,
                         int n_obs,
                         int grid_size = 25) {
  
  const double LARGE_VAL = 1e300; // fallback per NaN/Inf
  arma::vec grid(grid_size);
  
  // 1. Costruisci griglia log-scale
  for (int i = 0; i < grid_size; ++i) {
    double log_ax = std::log(ax);
    double log_bx = std::log(bx);
    grid[i] = std::exp(log_ax + i * (log_bx - log_ax) / (grid_size - 1));
  }
  
  // 2. Valuta la funzione sulla griglia
  arma::vec fvals(grid_size);
  for (int i = 0; i < grid_size; ++i) {
    double val = gcv_fun(grid[i], lambda, Q, Gom, ZidOi, ZidMi, n_obs);
    if (!arma::is_finite(val)) val = LARGE_VAL;
    fvals[i] = val;
  }
  
  // 3. Trova minimi locali sulla griglia
  std::vector<int> cand_idx;
  for (int i = 1; i < grid_size - 1; ++i) {
    if (fvals[i] < fvals[i-1] && fvals[i] < fvals[i+1]) {
      cand_idx.push_back(i);
    }
  }
  
  // 4. Applica Brent locale su ciascun candidato
  double best_x = grid[0];
  double best_f = fvals[0];
  
  for (size_t k = 0; k < cand_idx.size(); ++k) {
    int idx = cand_idx[k];
    int left = std::max(0, idx - 2);
    int right = std::min(grid_size - 1, idx + 2);
    
    double a_local = grid[left];
    double b_local = grid[right];
    
    double xstar = Brent_fmin(a_local, b_local, tol,
                              lambda, Q, Gom, ZidOi, ZidMi, n_obs);
    double fstar = gcv_fun(xstar, lambda, Q, Gom, ZidOi, ZidMi, n_obs);
    if (!arma::is_finite(fstar)) fstar = LARGE_VAL;
    
    if (fstar < best_f) {
      best_f = fstar;
      best_x = xstar;
    }
  }
  
  // 5. Fallback: se non ci sono minimi locali, prendi il minimo della griglia
  if (cand_idx.size() == 0) {
    arma::uword idx_min = fvals.index_min();
    best_x = grid[idx_min];
  }
  
  return best_x;
}

// [[Rcpp::export]]
arma::mat reconsX_internal(const arma::umat& O,
                           const arma::mat& Z,
                           const arma::uvec& id_obs,
                           const arma::mat& G,
                           arma::uword i) {
  
  // Estrazione Oi e Mi
  arma::uvec Oi = arma::find(O.row(i) == 1);
  arma::uvec Mi = arma::find(O.row(i) == 0);
  
  // Costruzione ZidOi e ZidMi
  arma::mat ZidOi = Z.submat(id_obs, Oi);
  arma::mat ZidMi = Z.submat(id_obs, Mi);
  
  // Estrazione Gom e Goo
  arma::mat Gom = G.submat(Oi, Mi);
  arma::mat Goo = G.submat(Oi, Oi);

  // Diagonalizzazione Goo
  arma::vec lambda;
  arma::mat Q;
  arma::eig_sym(lambda, Q, Goo);
  lambda = arma::clamp(lambda, 0.0, arma::datum::inf);

  // Ottimizzazione alpha
  // double alpha_opt = Brent_fmin(1e-12, 1.0, 1e-12,
  //                               lambda, Q, Gom, ZidOi, ZidMi, id_obs.n_elem);
  double alpha_opt = Brent_fmin_robust(2.22e-16, 100.0, 1e-6,
                                       lambda, Q, Gom, ZidOi, ZidMi, id_obs.n_elem);
  if(alpha_opt == 100.0) alpha_opt = 1e-12;
  // double alpha_opt = gradient_descent(.01, 1e-12, 10.0, lambda, Q, Gom, ZidOi, ZidMi,
  //                                     id_obs.n_elem, 1e-2, 1e-8, 1000);

  // Ricostruzione curva
  arma::mat Ginv = compute_Ginv(lambda, alpha_opt, Q);
  arma::mat A = Ginv * Gom;
  
  // Estrazione ZOi
  arma::uvec row_idx = { i };
  arma::mat Zi = Z.row(i);
  Zi.cols(Mi) = (Zi.cols(Oi) * A);

  return Zi;
}

// [[Rcpp::export]]
Rcpp::List dfs(int v, const arma::umat& Adj) {
  int p = Adj.n_rows;

  // Vettori di output
  arma::uvec Ck(p, arma::fill::zeros);
  int pk = 0;

  // Stack e visited
  std::vector<int> stack;
  std::vector<bool> visited(p, false);

  visited[v] = true;
  Ck(pk) = v;
  pk++;
  stack.push_back(v);

  while (!stack.empty()) {
    int i = stack.back();
    bool isolated = true;

    for (int j = 0; j < p; ++j) {
      if (!visited[j] && Adj(i, j)) {
        visited[j] = true;
        Ck(pk) = j;
        pk++;
        stack.push_back(j);
        isolated = false;
        break;  // simula `exit` in Fortran
      }
    }

    if (isolated) stack.pop_back();
  }

  // Riduci Ck al numero effettivo di nodi visitati
  arma::uvec Ck_trimmed = Ck.head(pk);

  return Rcpp::List::create(
    Rcpp::Named("Ck") = Ck_trimmed,
    Rcpp::Named("pk") = pk
  );
}

// [[Rcpp::export]]
Rcpp::List graph_adjacency(const arma::cube& S,
                           const arma::cube& W,
                           const arma::vec& fk,
                           double rho,
                           double alpha) {

  int p = S.n_rows;
  int N = S.n_slices;

  // Step 1: Calcolo normsoftTht
  arma::mat normsoftTht(p, p, arma::fill::zeros);
  double rho2_scalar = (1.0 - alpha) * rho;
  arma::mat rho2 = rho2_scalar * arma::ones(p, p);

  for (int k = 0; k < N; ++k) {
    arma::mat temp = fk(k) * arma::abs(S.slice(k)) - alpha * rho * W.slice(k);
    temp.for_each([](double& val) { if (val < 0) val = 0; });
    normsoftTht += arma::square(temp);
  }

  // Step 2: Costruzione matrice di adiacenza binaria
  arma::umat Adj(p, p, arma::fill::zeros);
  int ncount = 0;

  for (int j = 0; j < p; ++j) {
    Adj(j, j) = 1;
    for (int i = j + 1; i < p; ++i) {
      if (normsoftTht(i, j) > std::pow(rho2(i, j), 2)) {
        Adj(i, j) = 1;
        Adj(j, i) = 1;
        ncount++;
      }
    }
  }

  // Step 3: Casi speciali
  arma::uvec Ck;
  arma::ivec pk;

  int maxncount = (p - 1) * p / 2;
  if (ncount == 0) {
    Ck = arma::regspace<arma::uvec>(0, p - 1);
    pk = arma::ivec(p, arma::fill::ones);
    return Rcpp::List::create(
      Rcpp::Named("k") = p,
      Rcpp::Named("Ck") = Ck,
      Rcpp::Named("pk") = pk
    );
  }

  if (ncount == maxncount) {
    Ck = arma::regspace<arma::uvec>(0, p - 1);
    pk = arma::ivec(p, arma::fill::zeros);
    pk(0) = p;
    return Rcpp::List::create(
      Rcpp::Named("k") = 1,
      Rcpp::Named("Ck") = Ck,
      Rcpp::Named("pk") = pk
    );
  }

  // Step 4: DFS per componenti connesse
  arma::uvec connected(p, arma::fill::zeros);
  std::vector<int> pk_vec;
  std::vector<int> Ck_accum;
  int k = 0;

  for (int v = 0; v < p; ++v) {
    if (!connected[v]) {
      Rcpp::List dfs_result = dfs(v, Adj);
      arma::uvec component = dfs_result["Ck"];
      int pk_k = dfs_result["pk"];

      for (int i = 0; i < pk_k; ++i)
        connected[component[i]] = 1;

      pk_vec.push_back(pk_k);
      Ck_accum.insert(Ck_accum.end(), component.begin(), component.end());

      k++;
      if ((int)Ck_accum.size() == p) break;
    }
  }

  Ck = arma::conv_to<arma::uvec>::from(Ck_accum);
  pk = arma::ivec(pk_vec);

  return Rcpp::List::create(
    Rcpp::Named("k") = k,
    Rcpp::Named("Ck") = Ck,
    Rcpp::Named("pk") = pk
  );
}

inline double soft_group_cpp(double s, double tau, double norm) {
  if (norm == 0.0) return 0.0;
  double scale = 1.0 - tau / norm;
  return (scale > 0.0) ? scale * s : 0.0;
}

void softmatrix_tht_joint_cpp(arma::mat& S, const arma::mat& Tau, const arma::mat& Norm) {
  int p = S.n_rows;

  for (int i = 0; i < p; ++i) {
    for (int j = i; j < p; ++j) {
      double s_ij = soft_group_cpp(S(i, j), Tau(i, j), Norm(i, j));
      S(i, j) = s_ij;
      S(j, i) = s_ij;  // simmetrico
    }
  }
}

inline double softc_cpp(double s, double tau) {
  if ((s > 0 && s > tau))      return s - tau;
  else if ((s < 0 && -s > tau)) return s + tau;
  else                          return 0.0;
}

void softmatrix_tht_cpp(arma::mat& S, const arma::mat& Tau) {
  int p = S.n_rows;

  for (int i = 0; i < p; ++i) {
    for (int j = i; j < p; ++j) {
      double s_ij = softc_cpp(S(i, j), Tau(i, j));
      S(i, j) = s_ij;
      S(j, i) = s_ij;  // simmetrico
    }
  }
}

arma::mat ridgec_joint_cpp(double fk, const arma::mat& S, double lam) {
  // int p = S.n_rows;

  // Decomposizione spettrale
  arma::vec eigval;
  arma::mat eigvec;

  arma::eig_sym(eigval, eigvec, S);

  arma::vec w = arma::sqrt(fk * (-eigval + arma::sqrt(arma::square(eigval) + 4 * lam / fk)) / (2 * lam));
  arma::mat A = eigvec * arma::diagmat(w);

  // Omega = A * Aᵗ
  arma::mat Omega = A * A.t();
  return Omega;
}

void glasso_trace_2_2(int i, int k) {
  Rcpp::Rcout << "\tconnected component number " << i << "/" << k << "\n";
}

void admm_trace_2(int nit, double r, double eps1, double s, double eps2, double rho) {
  Rcpp::Rcout << "\titer num " << std::setw(6) << nit
  << "\trho: "   << std::fixed << std::setprecision(4) << rho
  << "\tr: "     << std::setprecision(6) << r
  << "\teps1: "  << eps1
  << "\ts: "     << s
  << "\teps2: "  << eps2 << "\n";
}

// [[Rcpp::export]]
Rcpp::List admm_tht_joint(int p, int K,
                          const arma::vec& fk,
                          arma::cube& S,
                          arma::cube& Theta,
                          const arma::cube& wTht,
                          double lam, double alpha,
                          int maxit, double thr,
                          double rho, int trace) {

  arma::cube Z(p, p, K, arma::fill::zeros);
  arma::cube U(p, p, K, arma::fill::zeros);
  arma::cube Z2(p, p, K);

  int conv = 0, nit = 0;
  int criterion = 1;
  int iter = 0;

  arma::mat normsoftTht(p, p, arma::fill::zeros);

  while (criterion == 1 && iter < maxit) {
    Z2 = Z;
    normsoftTht.zeros();
    iter++;

    // Step 1: Ridge + Soft threshold
    for (int k2 = 0; k2 < K; ++k2) {
      arma::mat input = S.slice(k2) - (rho * Z.slice(k2) - U.slice(k2)) / fk(k2);
      Theta.slice(k2) = ridgec_joint_cpp(fk(k2), input, rho);

      arma::mat Z_k = Theta.slice(k2) + U.slice(k2) / rho;
      arma::mat Tau_k = alpha * lam * wTht.slice(k2) / rho;

      softmatrix_tht_cpp(Z_k, Tau_k);
      Z.slice(k2) = Z_k;

      normsoftTht += arma::square(Z_k);
    }

    normsoftTht = arma::sqrt(normsoftTht);

    // Step 2: Soft matrix joint + update U
    double s2 = 0.0, r2 = 0.0, eps1 = 0.0, eps2 = 0.0;

    for (int k2 = 0; k2 < K; ++k2) {
      arma::mat Tau_k = (1 - alpha) * lam * wTht.slice(k2) / rho;

      softmatrix_tht_joint_cpp(Z.slice(k2), Tau_k, normsoftTht);
      U.slice(k2) += rho * (Theta.slice(k2) - Z.slice(k2));

      s2 += fk(k2) * arma::norm(rho * (Z2.slice(k2) - Z.slice(k2)), "fro");
      r2 += fk(k2) * arma::norm(Theta.slice(k2) - Z.slice(k2), "fro");

      eps1 += fk(k2) * (p * thr + thr * std::max(arma::norm(Theta.slice(k2), "fro"), arma::norm(Z.slice(k2), "fro")));
      eps2 += fk(k2) * (p * thr + thr * arma::norm(U.slice(k2), "fro"));
    }

    criterion = ((r2 >= eps1) || (s2 >= eps2)) ? 1 : 0;

    // Optional trace logging
    if (trace == 9 || (trace == 2 && criterion == 0)) {
      admm_trace_2(iter, r2, eps1, s2, eps2, rho);
    }

    if (r2 > 10.0 * s2) rho *= 2.0;
    else if (s2 > 10.0 * r2) rho /= 2.0;
  }

  if (iter == maxit) conv = 1;
  Theta = Z;
  nit = iter;

  // Calcolo df
  arma::ivec df(K, arma::fill::zeros);
  for (int k2 = 0; k2 < K; ++k2) {
    for (int i = 0; i < p - 1; ++i) {
      for (int j = i + 1; j < p; ++j) {
        if (Theta(i, j, k2) != 0.0)
          df(k2)++;
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("Tht") = Theta,
    Rcpp::Named("df")  = df,
    Rcpp::Named("nit") = nit,
    Rcpp::Named("conv") = conv
  );
}


// [[Rcpp::export]]
Rcpp::List admm_tht_sub(int p, int N,
                        arma::vec fk,
                        arma::cube S,
                        arma::cube wTht,
                        bool pendiag,
                        double rho, double alpha,
                        int maxit, double thr,
                        arma::cube Tht,
                        int trace) {

  arma::uvec Ck;
  arma::ivec pk;
  int nit = 0;
  int conv = 0;

  // 🔍 Trova componenti connesse
  Rcpp::List adj = graph_adjacency(S, wTht, fk, rho, alpha);
  int k = Rcpp::as<int>(adj["k"]);
  Ck = Rcpp::as<arma::uvec>(adj["Ck"]);
  pk = Rcpp::as<arma::ivec>(adj["pk"]);

  // return Rcpp::List::create(
    //   Rcpp::Named("k")  = k,
    //   Rcpp::Named("Ck")  = Ck,
    //   Rcpp::Named("pk")  = pk
    // );

  // 🔁 Ciclo su componenti connesse
  for (int i = 0; i < k; ++i) {
    int p_k1 = pk[0];
    int p_k2 = p - p_k1;

    // Controllo dei limiti
    arma::uvec idx1 = Ck.subvec(0, p_k1 - 1);
    arma::uvec idx2;
    if (p_k2 > 0)
      idx2 = Ck.subvec(p_k1, p - 1);

    if (p_k1 == 1) {
      for (int k2 = 0; k2 < N; ++k2) {
        if (pendiag) {
          S(idx1[0], idx1[0], k2) += alpha * rho * wTht(idx1[0], idx1[0], k2);
          Tht(idx1[0], idx1[0], k2) = 1.0 / S(idx1[0], idx1[0], k2);
        }
      }
    } else {
      if (trace == 2 || trace == 9)
        glasso_trace_2_2(i + 1, k);  // 1-based per output utente

      // 🚧 Estrai sottoblocchi
      arma::cube S_k(p_k1, p_k1, N);
      arma::cube wTht_k(p_k1, p_k1, N);
      arma::cube Tht_k(p_k1, p_k1, N);

      for (int k2 = 0; k2 < N; ++k2) {
        S_k.slice(k2)     = S.slice(k2).submat(idx1, idx1);
        wTht_k.slice(k2)  = wTht.slice(k2).submat(idx1, idx1);
        Tht_k.slice(k2)   = Tht.slice(k2).submat(idx1, idx1);

        if (pendiag) {
          for (int m = 0; m < p_k1; ++m)
            S_k(m, m, k2) += alpha * rho * wTht_k(m, m, k2);
        }
      }

      // 🚀 Stima modello glasso con ADMM
      double tau = 2.0;
      Rcpp::List fit = admm_tht_joint(p_k1, N, fk, S_k, Tht_k, wTht_k,
                                      rho, alpha, maxit, thr, tau, trace);

      if (Rcpp::as<int>(fit["conv"]) != 0)
        conv = -1;
        // return Rcpp::List::create(Rcpp::Named("conv") = -1);

      nit += Rcpp::as<int>(fit["nit"]);
      Tht_k = Rcpp::as<arma::cube>(fit["Tht"]);

      for (int k2 = 0; k2 < N; ++k2)
        Tht.slice(k2).submat(idx1, idx1) = Tht_k.slice(k2);
    }

    // ⚙️ Azzera blocchi fuori diagonale
    if (p_k2 > 0) {
      for (int k2 = 0; k2 < N; ++k2) {
        Tht.slice(k2).submat(idx1, idx2).zeros();
        Tht.slice(k2).submat(idx2, idx1).zeros();
      }

      // 🔁 Riordina Ck e pk per ciclo successivo
      Ck.head(p_k2)      = idx2;
      Ck.tail(p_k1)      = idx1;
      pk.head(k - 1)     = pk.tail(k - 1);
      pk[k - 1]          = p_k1;
    }
  }

  // 📊 Calcola gradi di libertà
  arma::ivec df(N, arma::fill::zeros);
  for (int k2 = 0; k2 < N; ++k2) {
    for (int i = 0; i < p - 1; ++i)
      for (int j = i + 1; j < p; ++j)
        if (Tht(i, j, k2) != 0.0)
          df[k2]++;
  }

  return Rcpp::List::create(
    Rcpp::Named("Tht")  = Tht,
    Rcpp::Named("df")   = df,
    Rcpp::Named("nit")  = nit,
    Rcpp::Named("k")  = k,
    Rcpp::Named("Ck")  = Ck,
    Rcpp::Named("pk")  = pk,
    Rcpp::Named("conv")  = conv
  );
}

// [[Rcpp::export]]
Rcpp::List reconsX(const umat& O,
                   const mat& Z,
                   const uvec& id_obs,
                   const uvec& id_pobs,
                   const mat& Phi,
                   const cube& Sgm,
                   const cube& X,
                   const cube& Xi,
                   const vec& tp,
                   const vec& fk,
                   const cube& wTht,
                   bool pendiag,
                   double rho,
                   double alpha,
                   int maxit,
                   double thr,
                   const cube& Tht,
                   const int n_threads = 1) {
  
  int d = X.n_rows;
  // int n = X.n_cols;
  int p = X.n_slices;
  int K = Phi.n_cols;
  int N = id_pobs.n_elem;
  
  arma::mat G = kron_sum(Sgm, Phi);
  
  // Contenitore thread-safe per i risultati
  std::vector<arma::mat> X_temp(N);
  
#ifdef _OPENMP
  omp_set_num_threads(n_threads);
#pragma omp parallel for schedule(dynamic)
#endif
  for (int idx = 0; idx < N; ++idx) {
    try {
      unsigned int i = id_pobs[idx];
      X_temp[idx] = reconsX_internal(O, Z, id_obs, G, i);
    } catch (std::exception& ex) {
      Rcpp::Rcout << "⚠️ Errore nel thread " << idx << ": " << ex.what() << std::endl;
      X_temp[idx] = arma::mat(); // fallback sicuro
    }
  }
  
  arma::uword N_uw = static_cast<arma::uword>(N);
  arma::uword p_uw = static_cast<arma::uword>(p);
  arma::cube Xhat = X;
  arma::cube Xihat = Xi;
  
  arma::uvec jj = generate_jj(p, d);
  
  for (arma::uword idx = 0; idx < N_uw; ++idx) {
    arma::uword i = id_pobs[idx];
    
    for (arma::uword j = 0; j < p_uw; ++j) {
      arma::uvec id_j = arma::find(jj == j); // colonne relative al blocco j
      Xhat.slice(j).col(i) = X_temp[idx].cols(id_j).t(); // trasposta per compatibilità d x 1
    }
  }
  
  Xihat = integrate_cube(Xhat, Phi, tp);
  
  Rcpp::List tmp = compute_S_and_rho_max(Xihat);
  tmp["X"] = Xhat;
  tmp["Xi"] = Xihat;
  
  arma::cube S = tmp["S"];
  int trace = 0;
  Rcpp::List tmp2 = admm_tht_sub(p, K, fk, S, wTht, pendiag, rho, alpha,
                                 maxit, thr, Tht, trace);
  arma::cube Tht_hat = tmp2["Tht"];
  arma::vec df = tmp2["df"];
  arma::cube Sgm_hat(p, p, K, fill::zeros);
  double err = 0.0;
  for (int k = 0; k < K; ++k) {
    Sgm_hat.slice(k) = arma::inv_sympd(Tht_hat.slice(k));
    double frob = arma::norm(Sgm_hat.slice(k) - Sgm.slice(k), "fro");
    err += frob / (df[k] + p);
  }
  tmp2["Sgm"] = Sgm_hat;
  tmp["Theta"] = tmp2;
  tmp["err"] = err / K;
  
  return tmp;
}

