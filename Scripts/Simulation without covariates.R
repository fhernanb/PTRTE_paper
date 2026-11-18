library(gamlss)
library(DiscreteDists)

# In the PTRTE distribution
# mu > 0, so we'll use the log link function
# 0 < sigma < 1, so we'll use the logit link function

# Useful functions to the simulation study --------------------------------

# Function to obtain mu_hat and sigma_hat for a fixed value of n
simul_one <- function(size, true_mu, true_sigma) {
  y <- rPTRTE(n=size, mu=true_mu, sigma=true_sigma)
  mod <- NULL
  mod <- try(gamlss(y~1, sigma.fo=~1, family="PTRTE",
                    control=gamlss.control(n.cyc=2500, trace=FALSE)))
  if (class(mod)[1] == "try-error") {
    res <- c(NA, NA)
  }
  else {
    res <- c(exp(coef(mod, what="mu")),
             logit_inv(coef(mod, what="sigma")))
  }
  res
}

logit_inv <- function(x) exp(x)/(1+exp(x))

# To perform the simulation -----------------------------------------------

library(parSim)

# Instruction to simulate
parSim(
  ### SIMULATION CONDITIONS
  mu = c(2.88),
  sigma = c(0.88),
  
  n = seq(from=2000, to=10000, by=1000),
  
  replications = 1000,                                 # Repetitions
  write = TRUE,                              # Writing to a file
  name = "Simuls/case1_6", # Name of the file
  nCores = 1,                                # Number of cores to use
  
  expression = {
    res <- simul_one(size=n, true_mu=mu, true_sigma=sigma)
    
    # Results list:
    Results <- list(
      mu_hat = res[1],
      sigma_hat = res[2]
    )
    
    # Return:
    Results
  }
)

# To load the results -----------------------------------------------------
archivos <- list.files(pattern = "^case1.*\\.txt$", 
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

# Number of observations for each n
num <- datos %>% group_by(n, mu, sigma) %>% count()
num
mean(num$nn)
min(num$nn)

# To obtain the metrics mean and MSE
trim <- 0.05
res <- datos %>% 
  group_by(n, mu, sigma) %>% 
  summarise(mean_mu=mean(mu_hat, trim=trim, na.rm=TRUE),
            mean_sigma=mean(sigma_hat, trim=trim, na.rm=TRUE),
            mse_mu=mean((mu - mu_hat)^2, trim=trim, na.rm=TRUE),
            mse_sigma=mean((sigma - sigma_hat)^2, trim=trim, na.rm=TRUE),
            nobs=n()
  )

res

# Mean -----------------------------------------------------

# Case 1

true_mu <- 2.88
true_sigma <- 0.88

p1 <- ggplot(data=res, aes(x=n, y=mean_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(mu))) +
  geom_line(y=true_mu, col="red", lty="dashed")

p1

p2 <- ggplot(data=res, aes(x=n, y=mean_sigma)) + 
  geom_line() + 
  labs(x="n", y=expression(hat(sigma))) +
  geom_line(y=true_sigma, col="red", lty="dashed")

p2

the_mean <- grid.arrange(p1, p2,
                         nrow = 1, ncol=2)
the_mean

ggsave(filename=paste0("Figs/case1_mean.pdf"), 
       plot=the_mean, 
       width=10, height=4)


# MSE -----------------------------------------------------

p1 <- ggplot(data=res, aes(x=n, y=mse_mu)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE(hat(mu))))

p1

p2 <- ggplot(data=res, aes(x=n, y=mse_sigma)) + 
  geom_line() + 
  labs(x="n", y=expression(MSE(hat(sigma))))

p2

the_mse <- grid.arrange(p1, p2,
                         nrow = 1, ncol=2)
the_mse

ggsave(filename=paste0("Figs/case1_mse.pdf"), 
       plot=the_mse, 
       width=10, height=4)





