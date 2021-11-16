library(MASS)
library(glmnet)

tol <- 1e-5

ln <- function(coefs, x, a) {
    X <- cbind(1, x)
    eta <- as.numeric((X %*% coefs))
    sum(
        -a*eta + log(1 + exp(eta))
    )
}

dln <- function(coefs, x, a, intercept=F) {
    if(intercept) {
        X <- cbind(1, x)
    } else {
        X <- x
    }
    
    deta <- t(X)
    eta <- as.numeric((X %*% coefs))
    mu <- plogis(eta)
    crossprod(X, (a - mu))
}

test_glmnet_fit <- function(
  X, A, intercept, pen, lambda1
) {
  logit_oal <- glmnet(X, A,
    family = binomial(), intercept = intercept,
    standardize = F,
    penalty.factor = pen
  )
  
  coefs_all <- coef(
    logit_oal, s = mean(pen) * lambda1,
    exact = T,
    x = X, y = A,
    family = binomial(), intercept = intercept,
    standardize = F,
    penalty.factor = pen
  )
}

# coefs_all should be unedited from glmnet
kkt_active <- function(coefs_all, X, A, intercept, pen, lambda1) {
  if(intercept){
    pen <- c(0, pen)
  } else {
    coefs_all <- coefs_all[-1]
  }
  sel <- which(abs(coefs_all) > 1e-8)
  sgn <- sign(coefs_all)[sel]
  
  dln(coefs_all, X, A, intercept)[sel]/nrow(X) -
    (pen[sel] * lambda1 * sgn)
}



n <- 200
p <- 10
rho <- 0
sig_x <- 0.25
lambda1 <- 0.02

### simulate data
set.seed(2021)
Sigma_x <- matrix(rho * sig_x^2, nrow = p, ncol = p)
diag(Sigma_x) <- sig_x^2
Mean_x <- rep(0, p)
X <- as.matrix(mvrnorm(n = n, mu = Mean_x, Sigma = Sigma_x, empirical = FALSE))

gA_x <- as.numeric(rowMeans(X[,1:4]))
pA <- plogis(gA_x)
summary(pA)

A <- as.numeric(runif(n = length(pA)) < pA) # simulate A
X <- scale(X, T, T)

#### Test 0: Verify that dln truly captures the score info
coef <- coef(glm(A ~ X[,1:4], family = binomial()))
stopifnot(max(abs(
  dln(coef, X[,1:4], A, T)
)) < tol)

##### Test 1: Verify lambda scaling without penalty.factor
pen <- rep(1,p)
lambda1 <- 0.02
coefs_all <- test_glmnet_fit(
  X, A, intercept = T, pen = pen, lambda1 = lambda1
)
stopifnot(max(abs(
kkt_active(coefs_all, X, A, intercept=T, pen=pen, lambda1=lambda1)
)) < tol)


##### Test 2: Verify lambda scaling with penalty.factor
pen <- rep(1,p)
pen[5] <- 10
pen[6] <- 200
coefs_all <- test_glmnet_fit(
  X, A, intercept = T, pen = pen, lambda1 = lambda1
)
stopifnot(max(abs(
kkt_active(coefs_all, X, A, intercept=T, pen=pen, lambda1=lambda1)
)) < tol)


##### Test 3: Verify lambda scaling with penalty.factor, no intercept
pen <- rep(1,p)
pen[5] <- 10
pen[6] <- 200
coefs_all <- test_glmnet_fit(
  X, A, intercept = F, pen = pen, lambda1 = lambda1
)
stopifnot(max(abs(
kkt_active(coefs_all, X, A, intercept=F, pen=pen, lambda1=lambda1)
)) < tol)

pen <- abs(betaXY)^{grid_min$gamma}
lambda1 <- grid_min$lambda1

stopifnot(max(abs(
kkt_active(coefs_all, X, A, intercept=T, pen=pen, lambda1=lambda1)
)) < tol)
