## generate data:
set.seed(32198734)

x <- as.matrix(expand.grid(seq(-1,1,length.out = 7),
                           seq(-1,1,length.out = 7)))
# matrix(runif(2*100,min=-1,max = 1),ncol = 2)
x_pred <- as.matrix(expand.grid(seq(-1,1,length.out = 9),
                                seq(-1,1,length.out = 9)))
f <- function(x, sd = 0.5){
  nn <- nrow(x)
  ifelse(x[,1] > 0, 5,0) + 
    ifelse(x[,2] > 0, 5,0) + 
    rnorm(n=nn,sd=sd)
}


df1 <- data.frame()# cbind(x,y0)
p <- 1 
nrows_train <- nrow(x)
train_ind <- as.data.frame((1:nrows_train)-1)
nrows_test <- nrow(x_pred)


y_0 <- as.matrix(f(x),ncol = 1)
y_sim <- as.matrix(f(x_pred),ncol=1)
n_init <- nrow(df1)

if(!require('mvtnorm')){
  install.packages('mvtnorm')
}
source('function_several_layers.R')
s1 <- stable_KP(x,c(y_0),x_pred)

## use median because mean might not be well-defined
preds <- apply(s1$y_new[1001:3000,],2,FUN = median)
p05 <- apply(s1$y_new[1001:3000,],2,FUN = quantile,probs=0.05)
p95 <- apply(s1$y_new[1001:3000,],2,FUN = quantile,probs=0.95)

c('RMSE'= sqrt(mean((y_sim-preds)^2)),
  'MAE' = mean(abs(y_sim-preds)))

resids_matrix <- cbind(p05-y_sim,
                         p95-y_sim)
covered_i <- rep(NA,81)
for(i in 1:81){
  covered_i[i] <- 0>= resids_matrix[i,1] & 0<=resids_matrix[i,2]
}


## coverage in this sample
mean(covered_i) 

if(!require('ggplot2')){
  install.packages('ggplot2')
}
require(ggplo2)

data.frame(x1 = x_pred[,1],
           x2 = x_pred[,2],
           preds = preds) |>
  ggplot2::ggplot(ggplot2::aes(x=x1,y=x2,fill = preds)) + 
  ggplot2::geom_tile() + 
  ggplot2::scale_fill_viridis_c(option = 'B') + 
  ggplot2::theme_minimal() + 
  ggplot2::labs(fill = 'Predictions')

data.frame(x1 = x_pred[,1],
           x2 = x_pred[,2],
           Residual = y_sim - preds,
           Prediction= preds,
           p95_resid = p95-preds,
          p05_resid = p05-preds
           ) |>
  ggplot2::ggplot(ggplot2::aes(x=Prediction,y=Residual)) + 
  ggplot2::geom_point() +
  ggplot2::geom_errorbar(ggplot2::aes(x=Prediction,ymin = p05_resid,
                                      ymax=p95_resid),color='blue') + 
  ggplot2::theme_minimal()





