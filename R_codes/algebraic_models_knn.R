## Calculate the variance based sensitivity indices in presence of correlation with the sobolshap_knn
## method from the sensitivity package in CRAN.
# https://cran.r-project.org/web/packages/sensitivity/index.html

# U=1 -> main eff
# U=0 -> total eff
# U=NULL and return_shap=TRUE -> shapley effects

library(sensitivity)
library (ggplot2)
library(mvtnorm)
library(whitening)
library(reshape2)

### define functions and select the model -------------------------------------------------------------------------------------

model_alg <- function(n_mod, X){
  
  if(n_mod==1){
    Y <- X[,1] + X[,2] + X[,2]*X[,3]
  }else if (n_mod==2){
    Y <- X[,1] + X[,2] + X[,1]*X[,3]
  }else{
    Y <- X[,1] + X[,2] + X[,3] + X[,4]
  }
  
  return(Y)
}

# to change the model assign different values of model_select (1,2,3)
model_select <- 1

### one evaluation of the simple algebraic model with fixed correlation coefficient ---------------------------------------------------------
# to change the correlation assign different values to ro
# HP: 
#    - all the factors were considered normal with mean=0 and sd=1
#    - correlation was considered between X1 and X4
#    - X4 does not appear in models 1 and 2
#    - here we have not considered the latent variable


ro <- 0.9 # change the correlation coefficient

n <- 50000
d <- 4
mu1 <- rep(0, d)
sig1 <- c(1,1,1,1)
Cormat1 <- matrix(c(1,0,0,ro,0,1,0,0,0,0,1,0,ro,0,0,1), d, d)
Covmat1 <- ( sig1 %*% t(sig1) ) * Cormat1
Xall_ind <- function(n) mvtnorm::rmvnorm(n,mu1,Covmat1)
X <- Xall_ind(n)

Y <- model_alg(model_select, X)

x_main <- sobolshap_knn(model = NULL, X = X, U = 1, method = "knn", n.knn = 5)
tell(x_main, Y)

x_total <- sobolshap_knn(model = NULL, X = X, U = 0, method = "knn", n.knn = 5)
tell(x_total, Y)


### variable correlation coefficient ------------------------------------------------------------------

ro_vect <- seq(from=-0.9, to=0.9, by=0.1)
lrv <- length(ro_vect)
S_main <- matrix(rep(0, lrv*d), nrow = lrv)
S_total <- matrix(rep(0, lrv*d), nrow = lrv)

for(i in 1:lrv){
  
  ro_i  <- ro_vect[i]
  
  Cormat1 <- matrix(c(1,0,0,ro_i,0,1,0,0,0,0,1,0,ro_i,0,0,1), d, d)
  Covmat1 <- ( sig1 %*% t(sig1) ) * Cormat1
  Xall_ind <- function(n) mvtnorm::rmvnorm(n,mu1,Covmat1)
  X <- Xall_ind(n)
  
  Y <- model_alg(model_select, X)
  
  x_main <- sobolshap_knn(model = NULL, X = X, U = 1, method = "knn", n.knn = 5)
  tell(x_main, Y)
  
  x_total <- sobolshap_knn(model = NULL, X = X, U = 0, method = "knn", n.knn = 5)
  tell(x_total, Y)
  
  S_main[i,] <- x_main$S
  S_total[i,] <- x_total$S
  
}


# draw some plots
colnames(S_main) <- c('X1_main','X2_main','X3_main','X4_main')
colnames(S_total) <- c('X1_total','X2_total','X3_total','X4_total')

df_main <- as.data.frame(S_main[,c(1,2,3,4)])
df_total <- as.data.frame(S_total[,c(1,2,3,4)])
df1 <- cbind(df_main, df_total)
df1$ro <- ro_vect
df2 <- melt(df1, id.vars="ro")

cols <- c("X1 main"="blue","X1 total"="blue","X2 main"="red","X2 total"="red", "X3 main")


ggplot(data=df2, aes(x=ro, y=value, group=variable)) +
  geom_line(aes(linetype=variable, color=variable), size=0.7) + 
  geom_point(aes(shape=variable, color=variable), size=2) + 
  scale_color_manual(values=c('dodgerblue3','brown3','goldenrod1','darkorchid3','dodgerblue3','brown3','goldenrod1','darkorchid3')) + 
  scale_linetype_manual(values = c('dashed','dashed','dashed','dashed','solid','solid','solid','solid')) + 
  scale_shape_manual(values = c(3,3,3,3,5,5,5,5))



