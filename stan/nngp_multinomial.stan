/* IIA multinomial response NNGP
   Three alternatives: 1 = car, 2 = transit, 3 = non-motorized.
   Non-motorized utility is fixed at 0. Car and transit contrasts each
   have an NNGP, independent and sharing (sigma, phi). Nugget tausq = 1.
*/

functions {
  real nngp_lpdf(vector Y, vector X_beta, real sigmasq, real tausq,
                 real phi, matrix NN_dist, matrix NN_distM, array[,] int NN_ind,
                 int N, int M) {

    vector[N] V;
    vector[N] YXb = Y - X_beta;
    vector[N] U = YXb;
    real kappa_p_1 = tausq / sigmasq + 1;
    int dim;
    int h;

    for (i in 2:N) {
      matrix[i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1) : M]
        iNNdistM;
      matrix[i < (M + 1) ? (i - 1) : M, i < (M + 1) ? (i - 1) : M]
        iNNCholL;
      vector[i < (M + 1) ? (i - 1) : M] iNNcorr;
      vector[i < (M + 1) ? (i - 1) : M] v;
      row_vector[i < (M + 1) ? (i - 1) : M] v2;
      dim = (i < (M + 1)) ? (i - 1) : M;

      if (dim == 1) {
        iNNdistM[1, 1] = kappa_p_1;
      } else {
        h = 0;
        for (j in 1:(dim - 1)) {
          for (k in (j + 1):dim) {
            h = h + 1;
            iNNdistM[j, k] = exp(-phi * NN_distM[(i - 1), h]);
            iNNdistM[k, j] = iNNdistM[j, k];
          }
        }
        for (j in 1:dim) {
          iNNdistM[j, j] = kappa_p_1;
        }
      }

      iNNCholL = cholesky_decompose(iNNdistM);
      iNNcorr = to_vector(exp(-phi * NN_dist[(i - 1), 1:dim]));

      v = mdivide_left_tri_low(iNNCholL, iNNcorr);

      V[i] = kappa_p_1 - dot_self(v);

      v2 = mdivide_right_tri_low(v', iNNCholL);

      U[i] = U[i] - v2 * YXb[NN_ind[(i - 1), 1:dim]];
    }
    V[1] = kappa_p_1;
    return -0.5 * (1 / sigmasq * dot_product(U, (U ./ V)) +
                   sum(log(V)) + N * log(sigmasq));
  }
}

data {
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> P;
  int<lower=1> S;
  array[N] int<lower=1, upper=3> y;
  array[N] int<lower=0, upper=1> mask;
  matrix[N, P] X_car;
  matrix[N, P] X_public;
  matrix[N, P] X_non;
  matrix[N, S] Socio;
  array[N - 1, M] int NN_ind;
  matrix[N - 1, M] NN_dist;
  matrix[N - 1, (M * (M - 1) / 2)] NN_distM;
  vector[P] uB;
  matrix[P, P] VB;
  vector[S] ugm;
  matrix[S, S] Vgm;
  real ss;
  real st;
  real ap;
  real bp;
}

transformed data {
  real tausq = 1;
  vector[N] zeros = rep_vector(0, N);
  cholesky_factor_cov[P] L_VB = cholesky_decompose(VB);
  cholesky_factor_cov[S] L_Vgm = cholesky_decompose(Vgm);
}

parameters {
  vector[P] beta;
  vector[S] gm;
  real asc_car;
  real asc_public;
  matrix[N, 2] u_raw;
  real<lower=0> sigma;
  real<lower=0> phi;
}

transformed parameters {
  real sigmasq = square(sigma);
  matrix[N, 3] u = append_col(u_raw, zeros);
}

model {
  beta ~ multi_normal_cholesky(uB, L_VB);
  gm ~ multi_normal_cholesky(ugm, L_Vgm);
  asc_car ~ normal(0, 1);
  asc_public ~ normal(0, 1);
  phi ~ gamma(ap, bp);
  sigma ~ normal(1, ss);
  u[:, 1] ~ nngp((X_car - X_non) * beta + Socio * gm + asc_car,
                 sigmasq, tausq, phi, NN_dist, NN_distM, NN_ind, N, M);
  u[:, 2] ~ nngp((X_public - X_non) * beta + asc_public,
                 sigmasq, tausq, phi, NN_dist, NN_distM, NN_ind, N, M);
  for (n in 1:N) {
    if (mask[n] == 1) {
      y[n] ~ categorical_logit(u[n]');
    }
  }
}

generated quantities {
  vector[N] log_lik;
  array[N] int y_est;

  for (n in 1:N) {
    log_lik[n] = categorical_logit_lpmf(y[n] | u[n]');
    y_est[n] = categorical_logit_rng(u[n]');
  }
}
