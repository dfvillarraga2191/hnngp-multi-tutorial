# End-to-end IIA multinomial HNNGP estimation with rstan
# Mode choice: car (1), transit (2), non-motorized (3)
#
# Usage (from the hnngp-multi-tutorial/ directory):
#   source("scripts/fit_hnngp.R")
#
# R/ess_rhats.R assumes chains = 4, iter = 3000, warmup = 1000.
# Change those three values together if you edit the stan() call.

library(rstan)
library(HDInterval)

source("R/NNMatrix.R")

# ---------------------------------------------------------------------------
# 1. Data
# ---------------------------------------------------------------------------
my_data <- read.csv("data/mode_multi.csv", header = TRUE, sep = ",",
                    dec = ".", check.names = FALSE)
my_data <- my_data[complete.cases(my_data), ]

# 1 = car (PMODE_R3 == 1), 2 = transit (== 2), 3 = non-motorized (== 97)
y <- rep(NA_integer_, nrow(my_data))
y[my_data$PMODE_R3 == 1]  <- 1L
y[my_data$PMODE_R3 == 2]  <- 2L
y[my_data$PMODE_R3 == 97] <- 3L
if (anyNA(y)) stop("Unexpected PMODE_R3 code.")

# Travel times in the file are seconds. /60 converts car and transit to minutes.
# Non-motorized time is the average of walking and bicycling, in minutes.
X_car <- cbind(
  time = my_data$T_Driving / 60,
  cost = my_data$Cost_Driving
)
X_public <- cbind(
  time = my_data$T_Transit / 60,
  cost = my_data$Cost_Transit
)
X_non <- cbind(
  time = (my_data$T_Walking + my_data$T_Bicycling) / 120,
  cost = 0
)
Socio <- cbind(
  high_inc  = as.integer(my_data$INCOM_R == 4),
  male      = as.integer(my_data$GENDER == "Male"),
  no_veh    = as.integer(my_data$HHVEH_R == 0),
  manhattan = as.integer(my_data$GEO_GROUP1_D == 1)
)

# Trip direction (destination - origin), as in the binary HNNGP tutorial
s_o <- as.matrix(my_data[, c("INTPTLAT", "INTPTLONG")])
s_d <- as.matrix(my_data[, c("DLAT", "DLONG")])
s   <- s_d
N   <- nrow(s)

# Optional: subsample for a faster demo (do this before the mask and the NN graph)
# set.seed(42); keep <- sample(N, 400)
# y <- y[keep]
# X_car <- X_car[keep, , drop = FALSE]
# X_public <- X_public[keep, , drop = FALSE]
# X_non <- X_non[keep, , drop = FALSE]
# Socio <- Socio[keep, , drop = FALSE]
# s <- s[keep, , drop = FALSE]
# N <- nrow(s)

# Train / test mask (1 = training), in the original row order
train_prop <- 0.8
set.seed(42)
mask <- sample(c(1L, 0L), N, replace = TRUE, prob = c(train_prop, 1 - train_prop))

# ---------------------------------------------------------------------------
# 2. Nearest-neighbor matrix (NNGP)
# ---------------------------------------------------------------------------
M <- 20
NN.matrix <- NNMatrix(
  coords = s,
  n.neighbors = M,
  n.omp.threads = 2,
  ord = order(s[, 1])
)

# ---------------------------------------------------------------------------
# 3. Priors
# ---------------------------------------------------------------------------
P  <- ncol(X_car)          # generic time and cost coefficients
S  <- ncol(Socio)          # socioeconomic shifters of the car utility
uB <- rep(0, P)
VB <- diag(P)
ugm <- rep(0, S)
Vgm <- diag(S)
ss <- 0.5                  # sigma ~ Normal(1, ss), truncated at 0 in Stan
st <- 0.001                # unused: tausq is fixed in the Stan file
ap <- 3; bp <- 0.5         # phi ~ Gamma(shape = ap, rate = bp)

# ---------------------------------------------------------------------------
# 4. Stan data
#    y, the attribute matrices, Socio, AND mask must share NN.matrix$ord.
# ---------------------------------------------------------------------------
ord <- NN.matrix$ord
stan_data <- list(
  N = N, M = M, P = P, S = S,
  y = y[ord],
  X_car = X_car[ord, , drop = FALSE],
  X_public = X_public[ord, , drop = FALSE],
  X_non = X_non[ord, , drop = FALSE],
  Socio = Socio[ord, , drop = FALSE],
  NN_ind = NN.matrix$NN_ind,
  NN_dist = NN.matrix$NN_dist,
  NN_distM = NN.matrix$NN_distM,
  mask = mask[ord],
  uB = uB, VB = VB, ugm = ugm, Vgm = Vgm,
  ss = ss, st = st, ap = ap, bp = bp
)

# ---------------------------------------------------------------------------
# 5. MCMC
# ---------------------------------------------------------------------------
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

parameters <- c(
  "beta", "gm", "asc_car", "asc_public",
  "sigmasq", "phi", "log_lik", "y_est"
)

fit <- stan(
  file = "stan/nngp_multinomial.stan",
  data = stan_data,
  pars = parameters,
  chains = 8,
  iter = 6000,
  warmup = 1000,
  refresh = 50,
  verbose = FALSE
)

dir.create("results", showWarnings = FALSE)
saveRDS(fit, "results/hnngp_multinomial_fit.rds")

print(fit, pars = c("beta", "gm", "asc_car", "asc_public", "sigmasq", "phi"))

# ---------------------------------------------------------------------------
# 6. Posterior summaries of interest
# ---------------------------------------------------------------------------
post <- rstan::extract(fit, pars = c("beta", "sigmasq", "phi"))

# Time is in minutes and cost in USD, so (beta_t / beta_c) * 60 is USD/hour.
vot <- (post$beta[, 1] / post$beta[, 2]) * 60
cat("\nVOT posterior mean:", mean(vot),
    " median:", median(vot), "\n")
cat("VOT 95% HDI:", paste(round(hdi(vot, credMass = 0.95), 2), collapse = " to "), "\n")

# ---------------------------------------------------------------------------
# 7. ESS and R-hat (see R/ess_rhats.R)
# ---------------------------------------------------------------------------
source("R/ess_rhats.R")
