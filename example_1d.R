## generate data:
set.seed(32198734)

x <- matrix(seq(-2,2,length.out = 40),ncol=1)
x_pred <- matrix(seq(-2,2,length.out = 100),ncol = 1)
# 1 jump
f <- function(xx,sd=0.5){
  nn <- length(xx)
  ifelse(xx >= 0,
         5,0) +
    rnorm(n = nn,sd = sd)
  
}  


df1 <- data.frame()# cbind(x,y0)
p <- 1 
nrows_train <- nrow(x)
train_ind <- as.data.frame((1:nrows_train)-1)
nrows_test <- nrow(x_pred)


y_0 <- as.matrix(f(x),ncol = 1)
y_sim <- as.matrix(f(x_pred),ncol=1)
n_init <- nrow(df1)

## add an intercept:
x <- cbind(1,x)
x_pred <- cbind(1,x_pred)


if(!require('mvtnorm')){
  install.packages('mvtnorm')
}
source('function_several_layers.R')
s1 <- stable_KP(x,c(y_0),x_pred)

## use median because mean might not be well-defined
## burn in of 1000 simulations
preds <- apply(s1$y_new[1001:3000,],2,FUN = median)
# quantiles at the 5% and 95% for a credible interval
p05 <- apply(s1$y_new[1001:3000,],2,FUN = quantile,probs=0.05)
p95 <- apply(s1$y_new[1001:3000,],2,FUN = quantile,probs=0.95)

# display fit with out-of-sample observations
plot(x_pred[,2],preds,type = 'l',ylim = range(p05,p95),
     xlab = 'x',ylab='predictions')
lines(x_pred[,2],p05,col='blue')
lines(x_pred[,2],p95,col='blue')
points(x_pred[,2],y_sim)

# out-of-sample error:
c('RMSE'=sqrt(mean((y_sim-preds)^2)),
  'MAE' = mean(abs(y_sim-preds)))

resids_matrix <- cbind(p05-y_sim, # we want this negative 
                       p95-y_sim  # we want this positive
                       )
covered_i <- rep(NA,100)
for(i in 1:100){
  covered_i[i] <- (0>= resids_matrix[i,1]) & (0<=resids_matrix[i,2])
}
## coverage in this sample
mean(covered_i) 



