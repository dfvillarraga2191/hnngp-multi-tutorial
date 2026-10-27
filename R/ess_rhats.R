# Effective sample size (ESS) and Gelman–Rubin R-hat
# IIA multinomial HNNGP
#
# Requires: sns, asbio
# Assumes an rstan fit object `fit` is in the workspace, with:
#   chains = 8, iter = 6000, warmup = 1000  (adjust below if different)

library(sns)
library(asbio)

iterations <- 6000   # total iterations per chain (including warmup)
warmup     <- 1000
chains     <- 8
post_start <- warmup + 1

c1 <- fit@sim[["samples"]][[1]]
c2 <- fit@sim[["samples"]][[2]]
c3 <- fit@sim[["samples"]][[3]]
c4 <- fit@sim[["samples"]][[4]]
c5 <- fit@sim[["samples"]][[5]]
c6 <- fit@sim[["samples"]][[6]]
c7 <- fit@sim[["samples"]][[7]]
c8 <- fit@sim[["samples"]][[8]]

# Parameter names as stored in the Stan samples list
par_list <- c(
  "beta[1]", "beta[2]",
  "gm[1]", "gm[2]", "gm[3]", "gm[4]",
  "asc_car", "asc_public",
  "sigmasq", "phi"
)

par_labels <- c(
  "beta[1]"    = "beta_t (time)",
  "beta[2]"    = "beta_c (cost)",
  "gm[1]"      = "gamma_Income",
  "gm[2]"      = "gamma_Male",
  "gm[3]"      = "gamma_Vehicle",
  "gm[4]"      = "gamma_Manhattan",
  "asc_car"    = "ASC_car",
  "asc_public" = "ASC_transit",
  "sigmasq"    = "sigmasq",
  "phi"        = "phi"
)

ess_rhat <- data.frame(
  parameter = character(0),
  ESS = numeric(0),
  R_hat = numeric(0),
  stringsAsFactors = FALSE
)

for (par in par_list) {
  post1 <- c1[[par]][post_start:iterations]
  post2 <- c2[[par]][post_start:iterations]
  post3 <- c3[[par]][post_start:iterations]
  post4 <- c4[[par]][post_start:iterations]
  post5 <- c5[[par]][post_start:iterations]
  post6 <- c6[[par]][post_start:iterations]
  post7 <- c7[[par]][post_start:iterations]
  post8 <- c8[[par]][post_start:iterations]

  ess_total <- ess(post1, method = c("coda", "ise")) +
               ess(post2, method = c("coda", "ise")) +
               ess(post3, method = c("coda", "ise")) +
               ess(post4, method = c("coda", "ise")) +
               ess(post5, method = c("coda", "ise")) +
               ess(post6, method = c("coda", "ise")) +
               ess(post7, method = c("coda", "ise")) +
               ess(post8, method = c("coda", "ise")) 

  rhat <- R.hat(cbind(post1, post2, post3, post4, post5, post6, post7, post8))

  label <- if (!is.null(par_labels[[par]])) par_labels[[par]] else par
  cat(label, "\n")
  cat("  ESS  :", ess_total, "\n")
  cat("  R_hat:", rhat, "\n\n")

  ess_rhat <- rbind(
    ess_rhat,
    data.frame(parameter = label, ESS = ess_total, R_hat = rhat)
  )
}

print(ess_rhat)
# Optional: write diagnostics to disk
# write.csv(ess_rhat, "results/ess_rhat.csv", row.names = FALSE)
