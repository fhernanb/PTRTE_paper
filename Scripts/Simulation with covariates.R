
library(gamlss)
library(DiscreteDists)

# In the PTRTE distribution
# mu > 0, so we"ll use the log link function
# 0 < sigma < 1, so we"ll use the logit link function

# The parameters ----------------------------------------------------------
# The next values correspond to the true betas in the regression model

true_b0_mu <- 0.8078  # intercept for mu
true_b1_mu <- 0.50  # slope for mu
true_g0_si <- 1.7424   # intercept for sigma
true_g1_si <- 0.5  # slope for sigma

# Useful functions to the simulation study --------------------------------

gendat <- function(n) {
  x1 <- runif(n)
  x2 <- runif(n)
  mu    <- exp(true_b0_mu + true_b1_mu * x1)
  sigma <- logit_inv(true_g0_si + true_g1_si* x2)
  y <- rPTRTE(n=n, mu=mu, sigma=sigma)
  data.frame(y=y, x1=x1, x2=x2)
}

logit_inv <- function(x) exp(x)/(1+exp(x))

# Function to obtain beta_hat for a fixed value of n
simul_one <- function(size) {
  dat <- gendat(size)
  mod <- NULL
  mod <- try(gamlss(y~x1, sigma.fo=~x2, family="PTRTE", data=dat,
                    control=gamlss.control(n.cyc=2500, trace=FALSE)),
             silent=TRUE)
  if (class(mod)[1] == "try-error")
    res <- rep(NA, 4)
  else
    res <- c(coef(mod, what="mu"), coef(mod, what="sigma"))
  res
}

# To perform the simulation -----------------------------------------------

library(parSim)

# Instruction to simulate 
parSim(
  ### SIMULATION CONDITIONS
  n = seq(from=1000, to=10000, by=1000),
  
  replications = 1000,                              # Repetitions
  write = TRUE,                           # Writing to a file
  name = "Simuls/case2_13", # Name of the file
  nCores = 1,                             # Number of cores to use
  
  expression = {
    res <- simul_one(size=n)
    
    # Results list:
    Results <- list(
      b0_mu_hat = res[1],
      b1_mu_hat = res[2],
      g0_si_hat = res[3],
      g1_si_hat = res[4]
    )
    
    # Return:
    Results
  }
)

# To load the results -----------------------------------------------------
archivos <- list.files(pattern = "^case2.*\\.txt$", 
                       path="Simuls",
                       full.names = TRUE)

archivos

lista_datos <- lapply(archivos, read.table, header = TRUE, 
                      sep = "", stringsAsFactors = FALSE)
datos <- do.call(rbind, lista_datos)

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

# To obtain the metrics mean and MSE
trim <- 0.10

res <- datos %>% 
  group_by(n) %>% 
  summarise(mean_b0_mu=mean(b0_mu_hat, trim=trim),
            mean_b1_mu=mean(b1_mu_hat, trim=trim),
            mean_g0_si=mean(g0_si_hat, trim=trim),
            mean_g1_si=mean(g1_si_hat, trim=trim),
            mse_b0_mu=mean((true_b0_mu - b0_mu_hat)^2, trim=trim), 
            mse_b1_mu=mean((true_b1_mu - b1_mu_hat)^2, trim=trim),
            mse_g0_si=mean((true_g0_si - g0_si_hat)^2, trim=trim), 
            mse_g1_si=mean((true_g1_si - g1_si_hat)^2, trim=trim),
            nobs=n())

res

mean(res$nobs)
min(res$nobs)

# Mean -----------------------------------------------------
k <- 0.35

p1 <- ggplot(data=res, aes(x=n, y=mean_b0_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(beta)[0])) +
  #ylim(true_b0_mu-k, true_b0_mu+k) +
  geom_line(y=true_b0_mu, col="red", lty="dashed")

p2 <- ggplot(data=res, aes(x=n, y=mean_b1_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(beta)[1])) +
  ylim(true_b1_mu-0.05, true_b1_mu+0.05) +
  geom_line(y=true_b1_mu, col="red", lty="dashed")

p3 <- ggplot(data=res, aes(x=n, y=mean_g0_si)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(gamma)[0])) +
  #ylim(true_g0_si-k, true_g0_si+k) +
  geom_line(y=true_g0_si, col="red", lty="dashed")

p4 <- ggplot(data=res, aes(x=n, y=mean_g1_si)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(gamma)[1])) +
  ylim(true_g1_si-k, true_g1_si+k) +
  geom_line(y=true_g1_si, col="red", lty="dashed")

mean2 <- grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
mean2
ggsave(filename="Figs/case2_mean.pdf", 
       plot=mean2, 
       width=10, height=8)

# MSE -----------------------------------------------------
p1 <- ggplot(data=res, aes(x=n, y=mse_b0_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE(hat(beta)[0])))

p2 <- ggplot(data=res, aes(x=n, y=mse_b1_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE(hat(beta)[1])))

p3 <- ggplot(data=res, aes(x=n, y=mse_g0_si)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE(hat(gamma)[0])))

p4 <- ggplot(data=res, aes(x=n, y=mse_g1_si)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE(hat(gamma)[1])))

mse2 <- grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
mse2
ggsave(filename="Figs/case2_mse.pdf", 
       plot=mse2, 
       width=10, height=8)



