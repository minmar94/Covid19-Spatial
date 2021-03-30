# Compute metrics (RelMSE, RootMSE, Coverage, etc.)
do_metric <- function(yhat, ytrue, lb, ub){
  
  cover <- mean(ytrue>=lb & ytrue<=ub)
  MSE <- mean((yhat-ytrue)^2)
  rMSE <- sqrt(MSE)
  RelMSE <- MSE/var(ytrue)
  PIW <- mean(ub-lb)
  
  return(round(c("Coverage" = cover, 
                 "MSE" = MSE, 
                 "Root_MSE" = rMSE,
                 "Rel_MSE" = RelMSE,
                 "PIW" = PIW), 
               5))
  
}

# extract post simulations for Y from Stan object
extract_postY <- function(fit){
  post <- fit@sim$samples[[1]]
  ypreds <- post[which(grepl("^Y_pred", names(post)))] %>% reduce(rbind)
  # burned
  ypredsQ <- rbind(apply(ypreds[,(M/2 + 1):M], 1, quantile, probs = c(0.025, 0.975)), rowMeans(ypreds[,(M/2 + 1):M])) %>% t
  return(list(post_samp = post, ypreds = ypreds, ypredsQ = ypredsQ))
}

# matplot of the last simulations
last_sim <- function(object, last = 2000, M = M, dat = dat){
  matplot(dat$t, object[,(M-2000):M], 
          type = "l", lwd = 0.1, col = "#00000008", lty = 1, 
          main = "Negative Binomial", 
          ylab = "Daily positives", 
          xlab = "Date")
  lines(dat$t, dat$Y, type = "b", pch = 19, col = "red")
}

# summary stats for model parameters
stats_pars <- function(pars){
  pars %>% 
    map_dfr(function(x){
      round(c(quantile(x, c(0.005, 0.025, 0.975, 0.995)), Mean = mean(x))[c(1,2,5,3,4)], 6)
    }) %>% mutate(Param = names(pars))
}