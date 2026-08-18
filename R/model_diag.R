
# Model diagnostics -------------------------------------------------------
# https://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html

library("bayesplot")
library("ggplot2")
library("rstanarm")

## 0. Table summary fit
write.csv(summary_fit, "summary_fit.csv", row.names = TRUE)

## 1.Convergence checks 
### R-hat and ESS (rstan)
summary_fit <- summary(fit)$summary
head(summary_fit)

summary_fit[, c("Rhat", "n_eff")]

### R-hat plot
rhat_vals <- summary_fit[, "Rhat"]
mcmc_rhat(rhat_vals)
### Trace plots
traceplot(fit, pars = c("beta", "gamma", "sigma_beta", "sigma_gamma"))
#traceplot(fit, pars = c("log_beta", "gamma", "sigma_beta", "sigma_gamma"))

### (cmdstan)
summary_fit <- fit$summary()
summary_fit[, c("rhat", "ess_bulk", "ess_tail")]

## R-hat plot
mcmc_rhat(summary_fit$rhat)

## Trace plots
mcmc_trace(fit$draws(c("beta", "gamma", "sigma_beta", "sigma_gamma")))
#mcmc_trace(fit$draws(c("log_beta", "log_gamma", "z[1]", "z_gamma[1]", "etasq", "rhosq")))

## 2.Divergence
sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)

n_divergent <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
n_divergent

n_treedepth <- sum(sapply(sampler_params, function(x) sum(x[, "treedepth__"] == 10)))
n_treedepth

energy <- sapply(sampler_params, function(x) x[, "energy__"])
bayesplot::bfmi(energy)

## 3.Posterior gemometry
### Pair plots
pairs(fit, pars = c("beta", "gamma"))

### Posterior dist
library("bayesplot")
library("ggplot2")
library("rstanarm")
# https://cran.r-project.org/web/packages/bayesplot/vignettes/plotting-mcmc-draws.html
mcmc_dens(as.array(fit), pars = c("beta", "gamma"))
mcmc_dens_overlay(as.array(fit), pars = c("beta", "gamma")) #chain overlaid

### (cmdstan)
mcmc_pairs(fit$draws(c("log_beta", "log_gamma")))
mcmc_dens_overlay(fit$draws(c("log_beta", "log_gamma")))

## 5.LOO-CV (Model evaluation)
library(loo)
log_lik <- extract_log_lik(fit)
loo_result <- loo(log_lik)

print(loo_result)

pdf("loo_diagnostic.pdf", width = 8, height = 6)
plot(loo_result)
dev.off()

#(cmdstan)
#log_lik <- fit$draws("log_lik")
#log_lik_matrix <- posterior::as_draws_matrix(log_lik)
#loo_result <- loo(log_lik_matrix)
#print(loo_result)

#pdf("loo_diagnostics2.pdf", width = 8, height = 6)
#plot(loo_result, diagnostic = "k")
#dev.off()

log_lik <- fit$draws("log_lik")
log_lik_matrix <- posterior::as_draws_matrix(log_lik)
loo_result <- loo(log_lik_matrix)

# Capture the printed loo output as text
loo_txt <- capture.output(print(loo_result))

# One PDF page with text on top and plot below
pdf("loo_2705604.pdf", width = 8, height = 11)

# Layout: top text panel, bottom plot panel
layout(matrix(c(1, 2), nrow = 2, byrow = TRUE), heights = c(2.5, 4.5))

par(mar = c(0, 0, 0, 0))

# Top panel: printed loo results
plot.new()
text(
  x = 0, y = 1,
  labels = paste(loo_txt, collapse = "\n"),
  adj = c(0, 1),
  cex = 0.8,
  family = "mono"
)

# Bottom panel: diagnostic plot
par(mar = c(5, 5, 2, 1))

plot(loo_result, diagnostic = "k")


dev.off()
