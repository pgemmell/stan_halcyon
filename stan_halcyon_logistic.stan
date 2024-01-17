functions {
  real Joe_F(array[,] int X, vector r, int L, int N, int E, int R, int S,
    matrix Q,
    row_vector pi_distro,
    array[,] int edge,
    vector edge_length) {
    real log_p_X = 0.0;
    array[S] matrix[N,4] lp_Xup; // prob subtree at i if nuc(i) = j, indexed by site
    for (s in 1:S) {
      for (l in 1:L) {
        for (x in 1:4) {
          if (x == X[l,s]) { lp_Xup[s][l,x] = log(1); } else { lp_Xup[s][l,x] = log(0.0); }
        }
      }
      int e = 1;
      while (e <= E) {
        int i = edge[e  ,1];
        int j = edge[e  ,2];
        int k = edge[e+1,2];
        real di = r[e  ] * edge_length[e  ];
        real dj = r[e+1] * edge_length[e+1];
        matrix[4,4] p_ij = matrix_exp(Q * di);
        matrix[4,4] p_ik = matrix_exp(Q * dj);
        for (x in 1:4) {
          real p_subj = dot_product(p_ij[x], exp(lp_Xup[s][j])); // x, j, k index rows
          real p_subk = dot_product(p_ik[x], exp(lp_Xup[s][k]));
          lp_Xup[s][i,x] = log(p_subj) + log(p_subk);
        }
        e += 2;
      }
      log_p_X += log_sum_exp(lp_Xup[s][R] + log(pi_distro));
    }
    return log_p_X;
  }
  matrix VCV(array[,] int Path, int L, int E, vector scaled_edge_length) {
    matrix[L,L] Sig;
    for (j in 1:L) {
      for (k in 1:j) {
        array[E] int common_edges;
        for (ii in 1:E) {
          common_edges[ii] =  Path[j,ii] * Path[k,ii];
        }
        real s = 0.0;
        for (ii in 1:E) {
          s = s + common_edges[ii] * scaled_edge_length[ii];
        }
        Sig[j,k] = s;
        Sig[k,j] = s;
      }
    }
    return Sig;
  }
  vector sig_fn(real h, real a, real b, vector r) {
    return h ./ (1 + exp(-a - b * log(r)));
  }
}

data {
  int<lower=0> L;
  int<lower=0> N;
  int<lower=0> E;
  int<lower=0> R;
  int<lower=0> S;
  array[L,S] int X; // alignment atgc todo: add 5 = gap
  matrix[4,4] Q;
  row_vector[4] pi_distro;
  array[E,2] int edge;
  vector[E] edge_length;
  array[N,E] int path_to_leaf;
  row_vector[L] y;
}

//transformed data {
//}

parameters {
  vector<lower=0>[E] r;
  real a;
  real b;
  real<lower=0> h;
  real yR;
}

model {
  h ~ uniform(0, 4);
  a ~ normal(0, 2);
  b ~ normal(0, 2);
  yR ~ normal(0, 1);
  r ~ lognormal(0, 1);
  target += Joe_F(X, r, L, N, E, R, S, Q, pi_distro, edge, edge_length);
  matrix[L,L] Sig = VCV(path_to_leaf, L, E, edge_length .* sig_fn(h, a, b, r));
  y ~ multi_normal(rep_vector(yR, L), Sig);
}

