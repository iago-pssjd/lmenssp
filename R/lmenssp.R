lmenssp <-
function(formula, data = NULL, id, process = "bm", timeVar, init = NULL, tol = 1e-5, maxiter = 100, silent = TRUE){ 
 
#library(nlme)
#library(MASS)

mf <- model.frame(formula = formula, data = data)
y  <- as.matrix(model.extract(mf, "response"))
x  <- as.matrix(model.matrix(attr(mf, "terms"), data = mf)[, -1])

nsubj  <- length(unique(id))
idlist <- unique(id)

Time  <- tapply(timeVar, id, function(x) x)
Intercept <- rep(1, nrow(data))
data2 <- data.frame(cbind(id, Intercept, x))
DM    <- split(data2[, -1], data2$id)
YM    <- tapply(y, id, function(x) x)
nobs  <- tapply(id, id, function(x) length(x))

#########################################################
################### BROWNIAN MOTION #####################
#########################################################

if(process == "bm"){

if (length(init) == 0){
data$timeVar <- timeVar
init <- as.numeric(VarCorr(lme(formula, random = ~ timeVar|id, method = "ML", data = data))[,1])
}

theta.new <- init
theta.old <- init * 20
tol       <- tol
iter      <- 1
Niter     <- maxiter

while (sqrt((theta.old - theta.new) %*% (theta.old - theta.new)) > tol & iter <= Niter){

theta.old <- theta.new

a <- theta.old[1]  ## omegasq
b <- theta.old[2]  ## sigmasq
c <- theta.old[3]  ## nu


### betahat

sum.left.beta <- sum.right.beta <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv <- chol2inv(chol(Vi))
Vi.inv  <- solve(Vi)
Xi.transp <- t(Xi)

sum.left.beta  <- sum.left.beta + Xi.transp %*% Vi.inv %*% Xi  
sum.right.beta <- sum.right.beta + Xi.transp %*% Vi.inv %*% Yi

}#for(i in idlist)

#b.hat     <- chol2inv(chol(sum.left.beta)) %*% sum.right.beta
b.hat      <- solve(sum.left.beta) %*% sum.right.beta

####### score for theta

### a = omegasq

sum.a <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv   <- chol2inv(chol(Vi))
Vi.inv <- solve(Vi)
ri       <- Yi - Xi %*% b.hat 

derVa <- Ji

sum.a <- sum.a + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVa))

}#for (i in idlist)


### b = sigmasq

sum.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv   <- chol2inv(chol(Vi))
Vi.inv <- solve(Vi)
ri       <- Yi - Xi %*% b.hat 

derVb <- Ri

sum.b <- sum.b + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVb))

}#for (i in idlist)


# c = tausq

sum.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv   <- chol2inv(chol(Vi))
Vi.inv <- solve(Vi)
ri       <- Yi - Xi %*% b.hat 

#derVc <- Ii
#sum.c  <- sum.c + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVc))
sum.c  <- sum.c + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv))

}#for (i in idlist)

score.theta <- 0.5 * matrix(c(sum.a, sum.b, sum.c))


######### INFORMATION MATRIX

#### 1) a = omegasq, 2) a = omegasq

sum.a.a <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv   <- chol2inv(chol(Vi))
Vi.inv <- solve(Vi)

derVa <- Ji

sum.a.a  <- sum.a.a + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVa))

}#for (i in idlist)


#### 1) a = omegasq, 2) b = sigmasq

sum.a.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii
#Vi.inv   <- chol2inv(chol(Vi))

Vi.inv <- solve(Vi)

derVa <- Ji
derVb <- Ri

sum.a.b  <- sum.a.b + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVb))

}#for (i in idlist)


#### 1) a = omegasq, 2) c = tausq

sum.a.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv   <- chol2inv(chol(Vi))

Vi.inv <- solve(Vi)

derVa <- Ji
#derVc <- Ii
#sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVc))

sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv))

}#for (i in idlist)


#### 1) b = sigmasq, 2) b = sigmasq

sum.b.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv   <- chol2inv(chol(Vi))

Vi.inv <- solve(Vi)

derVb <- Ri

sum.b.b  <- sum.b.b + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVb))

}#for (i in idlist)


#### 1) b = sigmasq, 2) c = tausq

sum.b.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv   <- chol2inv(chol(Vi))

Vi.inv <- solve(Vi)

derVb <- Ri
#derVc <- Ii
#sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVc))

sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv))

}#for (i in idlist)


#### 1) c = tausq, 2) c = tausq

sum.c.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv   <- chol2inv(chol(Vi))

Vi.inv <- solve(Vi)

#derVc <- Ii
#sum.c.c  <- sum.c.c + sum(diag(Vi.inv %*% derVc %*% Vi.inv %*% derVc))
sum.c.c  <- sum.c.c + sum(Vi.inv * t(Vi.inv))

}#for (i in idlist)

expected.hessian <- -0.5 * matrix(c(sum.a.a, sum.a.b, sum.a.c,
                                    sum.a.b, sum.b.b, sum.b.c,
                                    sum.a.c, sum.b.c, sum.c.c), 
                                  ncol = 3, byrow = T)

theta.new <- as.numeric(theta.old - ginv(expected.hessian) %*% score.theta)

if(silent == FALSE){

cat("iteration = ", iter, "\n")
cat("theta.old = ", theta.old, "\n")
cat("beta.hat  = ", b.hat, "\n")
cat("theta.new = ", theta.new, "\n")
cat("sqrt.diff=", sqrt((theta.old - theta.new) %*% (theta.old - theta.new)), "\n")      
print("-----------------------------")

}

iter <- iter + 1

}#while


theta.old <- theta.new

a <- theta.old[1]  ## omegasq
b <- theta.old[2]  ## sigmasq
c <- theta.old[3]  ## tausq

### betahat

sum.left.beta <- sum.right.beta <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv    <- chol2inv(chol(Vi))
Vi.inv <- solve(Vi)

Xi.transp <- t(Xi)

sum.left.beta  <- sum.left.beta + Xi.transp %*% Vi.inv %*% Xi  
sum.right.beta <- sum.right.beta + Xi.transp %*% Vi.inv %*% Yi

}#for(i in idlist)

b.hat    <- solve(sum.left.beta) %*% sum.right.beta
b.varcov <- solve(sum.left.beta)


######### INFORMATION MATRIX

#### 1) a = omegasq, 2) a = omegasq

sum.a.a <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv   <- chol2inv(chol(Vi))
Vi.inv <- solve(Vi)

derVa <- Ji

sum.a.a  <- sum.a.a + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVa))

}#for (i in idlist)


#### 1) a = omegasq, 2) b = sigmasq

sum.a.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii
#Vi.inv   <- chol2inv(chol(Vi))
Vi.inv <- solve(Vi)

derVa <- Ji
derVb <- Ri

sum.a.b  <- sum.a.b + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVb))

}#for (i in idlist)


#### 1) a = omegasq, 2) c = tausq

sum.a.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y)) 
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv   <- chol2inv(chol(Vi))
Vi.inv <- solve(Vi)

derVa <- Ji
#derVc <- Ii
#sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVc))

sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv))

}#for (i in idlist)


#### 1) b = sigmasq, 2) b = sigmasq

sum.b.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv   <- chol2inv(chol(Vi))
Vi.inv <- solve(Vi)

derVb <- Ri

sum.b.b  <- sum.b.b + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVb))

}#for (i in idlist)


#### 1) b = sigmasq, 2) c = tausq

sum.b.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv   <- chol2inv(chol(Vi))
Vi.inv <- solve(Vi)

derVb <- Ri
#derVc <- Ii
#sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVc))

sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv))

}#for (i in idlist)


#### 1) c = tausq, 2) c = tausq

sum.c.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y)) 
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv   <- chol2inv(chol(Vi))
Vi.inv <- solve(Vi)

#derVc <- Ii
#sum.c.c  <- sum.c.c + sum(diag(Vi.inv %*% derVc %*% Vi.inv %*% derVc))

sum.c.c  <- sum.c.c + sum(diag(Vi.inv %*% Vi.inv))

}#for (i in idlist)

expected.hessian <- -0.5 * matrix(c(sum.a.a, sum.a.b, sum.a.c,
                                    sum.a.b, sum.b.b, sum.b.c,
                                    sum.a.c, sum.b.c, sum.c.c), 
                                  ncol = 3, byrow = T)


sd.theta <- sqrt(diag(ginv(-expected.hessian)))

## loglik

sum.loglik <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) pmin(x, y))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

#Vi.inv <- chol2inv(chol(Vi))
Vi.inv <- solve(Vi)

ri     <- Yi - Xi %*% b.hat

#sum.loglik  <- sum.loglik + log(det(Vi)) + t(ri) %*% Vi.inv %*% ri

sum.loglik  <- sum.loglik + as.numeric(determinant(Vi, logarithm=TRUE)$modulus) +  t(ri) %*% Vi.inv %*% ri

}#for (i in idlist)

max.loglik <- as.numeric(- 0.5 * nrow(YM) * log(2 * pi) - 0.5 * sum.loglik)

result <- cbind(c(c(b.hat), theta.new) , c(sqrt(diag(b.varcov)), sd.theta),
                c(c(b.hat)/sqrt(diag(b.varcov)), NA, NA, NA))
result <- cbind(result, 2 * (1 - pnorm(abs(result[, 3]))))
colnames(result) <- c("Estimate", "Standard error", "Z-estimate", "p-value")
rownames(result)[(nrow(result) - 2) : nrow(result)] <- c("omegasq", "sigmasq", "tausq")

score.theta <- matrix(score.theta, nrow = 1)
colnames(score.theta) <- c("omegasq", "sigmasq", "tausq")

output           <- list()
output$title     <- "Mixed effects model with random intercept and Brownian motion"
output$date      <- date()
output$estimates <- result
output$maxloglik <- max.loglik 
output$score     <- score.theta
output
}

##########################################################
################### INTEGRATED BROWNIAN MOTION ###########
##########################################################

if(process == "ibm"){

if (length(init) == 0){
data$timeVar <- timeVar
init <- as.numeric(VarCorr(lme(formula, random = ~ timeVar|id, method = "ML", data = data))[,1])
}

theta.new <- init
theta.old <- init * 20
tol       <- tol
iter      <- 1
Niter     <- maxiter


while (sqrt((theta.old - theta.new) %*% (theta.old - theta.new)) > tol & iter <= Niter){

theta.old <- theta.new

a <- theta.old[1]  ## omegasq
b <- theta.old[2]  ## sigmasq
c <- theta.old[3]  ## nu


### betahat

sum.left.beta <- sum.right.beta <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv  <- solve(Vi)
Xi.transp <- t(Xi)

sum.left.beta  <- sum.left.beta + Xi.transp %*% Vi.inv %*% Xi  
sum.right.beta <- sum.right.beta + Xi.transp %*% Vi.inv %*% Yi

}#for(i in idlist)

b.hat <- solve(sum.left.beta) %*% sum.right.beta

####### score for theta

### a = omegasq

sum.a <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)
ri       <- Yi - Xi %*% b.hat 

derVa <- Ji

sum.a <- sum.a + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVa))

}#for (i in idlist)


### b = sigmasq

sum.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)
ri       <- Yi - Xi %*% b.hat 

derVb <- Ri

sum.b <- sum.b + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVb))

}#for (i in idlist)


# c = tausq

sum.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)
ri       <- Yi - Xi %*% b.hat 

derVc <- Ii

sum.c  <- sum.c + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVc))

}#for (i in idlist)

score.theta <- 0.5 * matrix(c(sum.a, sum.b, sum.c))


######### INFORMATION MATRIX

#### 1) a = omegasq, 2) a = omegasq

sum.a.a <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)

derVa <- Ji

sum.a.a  <- sum.a.a + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVa))

}#for (i in idlist)


#### 1) a = omegasq, 2) b = sigmasq

sum.a.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)

derVa <- Ji
derVb <- Ri

sum.a.b  <- sum.a.b + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVb))

}#for (i in idlist)


#### 1) a = omegasq, 2) c = tausq

sum.a.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)

derVa <- Ji
#derVc <- Ii
#sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVc))

sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv))

}#for (i in idlist)


#### 1) b = sigmasq, 2) b = sigmasq

sum.b.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)

derVb <- Ri

sum.b.b  <- sum.b.b + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVb))

}#for (i in idlist)


#### 1) b = sigmasq, 2) c = tausq

sum.b.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)

derVb <- Ri
#derVc <- Ii
#sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVc))

sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv))

}#for (i in idlist)


#### 1) c = tausq, 2) c = tausq

sum.c.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)

#derVc <- Ii
#sum.c.c  <- sum.c.c + sum(diag(Vi.inv %*% derVc %*% Vi.inv %*% derVc))
sum.c.c  <- sum.c.c + sum(Vi.inv * t(Vi.inv))

}#for (i in idlist)

expected.hessian <- -0.5 * matrix(c(sum.a.a, sum.a.b, sum.a.c,
                                    sum.a.b, sum.b.b, sum.b.c,
                                    sum.a.c, sum.b.c, sum.c.c), 
                                  ncol = 3, byrow = T)

theta.new <- as.numeric(theta.old - ginv(expected.hessian) %*% score.theta)

if(silent == FALSE){

cat("iteration = ", iter, "\n")
cat("theta.old = ", theta.old, "\n")
cat("beta.hat  = ", b.hat, "\n")
cat("theta.new = ", theta.new, "\n")
cat("sqrt.diff=", sqrt((theta.old - theta.new) %*% (theta.old - theta.new)), "\n")
cat("score=", score.theta, "\n")      
print("-----------------------------")

}

iter <- iter + 1

}#while


theta.old <- theta.new

a <- theta.old[1]  ## omegasq
b <- theta.old[2]  ## sigmasq
c <- theta.old[3]  ## tausq

### betahat

sum.left.beta <- sum.right.beta <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)

Xi.transp <- t(Xi)

sum.left.beta  <- sum.left.beta + Xi.transp %*% Vi.inv %*% Xi  
sum.right.beta <- sum.right.beta + Xi.transp %*% Vi.inv %*% Yi

}#for(i in idlist)

b.hat    <- solve(sum.left.beta) %*% sum.right.beta
b.varcov <- solve(sum.left.beta)

####### score for theta

### a = omegasq

sum.a <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)
ri       <- Yi - Xi %*% b.hat 

derVa <- Ji

sum.a <- sum.a + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVa))

}#for (i in idlist)


### b = sigmasq

sum.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)
ri       <- Yi - Xi %*% b.hat 

derVb <- Ri

sum.b <- sum.b + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVb))

}#for (i in idlist)


# c = tausq

sum.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)
ri       <- Yi - Xi %*% b.hat 

derVc <- Ii

sum.c  <- sum.c + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVc))

}#for (i in idlist)

score.theta <- 0.5 * matrix(c(sum.a, sum.b, sum.c))

######### INFORMATION MATRIX

#### 1) a = omegasq, 2) a = omegasq

sum.a.a <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)

derVa <- Ji

sum.a.a  <- sum.a.a + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVa))

}#for (i in idlist)


#### 1) a = omegasq, 2) b = sigmasq

sum.a.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)

derVa <- Ji
derVb <- Ri

sum.a.b  <- sum.a.b + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVb))

}#for (i in idlist)


#### 1) a = omegasq, 2) c = tausq

sum.a.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)

derVa <- Ji
#derVc <- Ii
#sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVc))

sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv))

}#for (i in idlist)


#### 1) b = sigmasq, 2) b = sigmasq

sum.b.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)

derVb <- Ri

sum.b.b  <- sum.b.b + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVb))

}#for (i in idlist)


#### 1) b = sigmasq, 2) c = tausq

sum.b.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)

derVb <- Ri
#derVc <- Ii
#sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVc))

sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv))

}#for (i in idlist)


#### 1) c = tausq, 2) c = tausq

sum.c.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)

#derVc <- Ii
#sum.c.c  <- sum.c.c + sum(diag(Vi.inv %*% derVc %*% Vi.inv %*% derVc))

sum.c.c  <- sum.c.c + sum(diag(Vi.inv %*% Vi.inv))

}#for (i in idlist)

expected.hessian <- -0.5 * matrix(c(sum.a.a, sum.a.b, sum.a.c,
                                    sum.a.b, sum.b.b, sum.b.c,
                                    sum.a.c, sum.b.c, sum.c.c), 
                                  ncol = 3, byrow = T)

sd.theta <- sqrt(diag(ginv(-expected.hessian)))

## loglik

sum.loglik <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 0.5 * pmin(x, y)^2 * (pmax(x, y) - pmin(x, y)/3))  
Ii <- diag(ni)

Vi <- a * Ji + b * Ri + c * Ii

Vi.inv <- solve(Vi)

ri     <- Yi - Xi %*% b.hat

sum.loglik  <- sum.loglik + as.numeric(determinant(Vi, logarithm=TRUE)$modulus) +  t(ri) %*% Vi.inv %*% ri

}#for (i in idlist)

max.loglik <- as.numeric(- 0.5 * nrow(YM) * log(2 * pi) - 0.5 * sum.loglik)

result <- cbind(c(c(b.hat), theta.new) , c(sqrt(diag(b.varcov)), sd.theta),
                c(c(b.hat)/sqrt(diag(b.varcov)), NA, NA, NA))
result <- cbind(result, 2 * (1 - pnorm(abs(result[, 3]))))
colnames(result) <- c("Estimate", "Standard error", "Z-estimate", "p-value")
rownames(result)[(nrow(result) - 2) : nrow(result)] <- c("omegasq", "sigmasq", "tausq")

score.theta <- matrix(score.theta, nrow = 1)
colnames(score.theta) <- c("omegasq", "sigmasq", "tausq")

output           <- list()
output$title     <- "Mixed effects model with random intercept and integrated Brownian motion"
output$date      <- date()
output$estimates <- result
output$maxloglik <- max.loglik 
output$score     <- score.theta
output

}


#####################################################################
################### INTEGRATED ORNSTEIN-UHLENBECK ###################
#####################################################################

if(process == "iou"){

if (length(init) == 0){
data$timeVar <- timeVar
lme.est <- as.numeric(VarCorr(lme(formula, random = ~ timeVar|id, method = "ML", data = data))[,1])
init    <- c(lme.est[1], lme.est[2]/2, lme.est[2]/2, lme.est[3]) 
}

theta.new <- init
theta.old <- init * 20
tol       <- tol
iter      <- 1
Niter     <- maxiter

while (sqrt((theta.old - theta.new) %*% (theta.old - theta.new)) > tol & iter <= Niter){

theta.old <- theta.new

a <- theta.old[1]  ## omegasq
b <- theta.old[2]  ## sigmasq
c <- theta.old[3]  ## alpha or nu
d <- theta.old[4]  ## tausq


### betahat

sum.left.beta <- sum.right.beta <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv    <- solve(Vi)
Xi.transp <- t(Xi)

sum.left.beta  <- sum.left.beta + Xi.transp %*% Vi.inv %*% Xi  
sum.right.beta <- sum.right.beta + Xi.transp %*% Vi.inv %*% Yi

}#for(i in idlist)

b.hat    <- solve(sum.left.beta) %*% sum.right.beta


####### score for theta

### a = omegasq

sum.a <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)
ri       <- Yi - Xi %*% b.hat 

derVa <- Ji

sum.a  <- sum.a + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVa))

}#for (i in idlist)


### b = sigmasq

sum.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv    <- solve(Vi)
ri       <- Yi - Xi %*% b.hat 

derVb <- outer(Timei, Timei, function(x,y) 
0.5*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)

sum.b <- sum.b + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVb))

}#for (i in idlist)


### c = alpha

sum.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv    <- solve(Vi)
ri       <- Yi - Xi %*% b.hat 

derVc <- outer(Timei, Timei, function(x,y) b * (
0.5*c^(-3)*(-4*pmin(x,y)-pmax(x,y)*exp(-c*pmax(x,y))-pmin(x,y)*exp(-c*pmin(x,y))+abs(pmax(x,y)-pmin(x,y))*exp(-c*abs(pmax(x,y)-pmin(x,y))))
+
1.5*c^(-4)*(-exp(-c*pmax(x,y))-exp(-c*pmin(x,y))+1+exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
)

sum.c  <- sum.c + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVc))

}#for (i in idlist)


# d = tausq

sum.d <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)
ri       <- Yi - Xi %*% b.hat 

#derVd <- Ii
#sum.d  <- sum.d + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv %*% derVd))
sum.d  <- sum.d + sum(diag(Vi.inv %*% (ri %*% t(ri) - Vi) %*% Vi.inv))

}#for

score.theta <- 0.5 * matrix(c(sum.a, sum.b, sum.c, sum.d))


######### INFORMATION MATRIX

#### 1) a = omegasq, 2) a = omegasq

sum.a.a <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVa <- Ji

sum.a.a  <- sum.a.a + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVa))

}#for (i in idlist)


#### 1) a = omegasq, 2) b = sigmasq

sum.a.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVa <- Ji
derVb <- outer(Timei, Timei, function(x,y) 
0.5*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)

sum.a.b  <- sum.a.b + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVb))

}#for (i in idlist)


#### 1) a = omegasq, 2) c = alpha

sum.a.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVa <- Ji
derVc <- outer(Timei, Timei, function(x,y) b * (
0.5*c^(-3)*(-4*pmin(x,y)-pmax(x,y)*exp(-c*pmax(x,y))-pmin(x,y)*exp(-c*pmin(x,y))+abs(pmax(x,y)-pmin(x,y))*exp(-c*abs(pmax(x,y)-pmin(x,y))))
+
1.5*c^(-4)*(-exp(-c*pmax(x,y))-exp(-c*pmin(x,y))+1+exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
)

sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVc))

}#for (i in idlist)


#### 1) a = omegasq, 2) d = tausq

sum.a.d <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVa <- Ji
#derVd <- Ii
#sum.a.d  <- sum.a.d + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVd))
sum.a.d  <- sum.a.d + sum(diag(Vi.inv %*% derVa %*% Vi.inv))

}#for (i in idlist)


#### 1) b = sigmasq, 2) b = sigmasq

sum.b.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVb <- outer(Timei, Timei, function(x,y) 
0.5*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)

sum.b.b  <- sum.b.b + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVb))

}#for (i in idlist)


#### 1) b = sigmasq, 2) c = alpha

sum.b.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVb <- outer(Timei, Timei, function(x,y) 
0.5*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
derVc <- outer(Timei, Timei, function(x,y) b * (
0.5*c^(-3)*(-4*pmin(x,y)-pmax(x,y)*exp(-c*pmax(x,y))-pmin(x,y)*exp(-c*pmin(x,y))+abs(pmax(x,y)-pmin(x,y))*exp(-c*abs(pmax(x,y)-pmin(x,y))))
+
1.5*c^(-4)*(-exp(-c*pmax(x,y))-exp(-c*pmin(x,y))+1+exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
)

sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVc))

}#for (i in idlist)


#### 1) b = sigmasq, 2) d = tausq

sum.b.d<- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVb <- outer(Timei, Timei, function(x,y) 
0.5*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
#derVd <- Ii
#sum.b.d  <- sum.b.d + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVd))
sum.b.d  <- sum.b.d + sum(diag(Vi.inv %*% derVb %*% Vi.inv))

}#for (i in idlist)


#### 1) c = alpha, 2) c = alpha

sum.c.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVc <- outer(Timei, Timei, function(x,y) b * (
0.5*c^(-3)*(-4*pmin(x,y)-pmax(x,y)*exp(-c*pmax(x,y))-pmin(x,y)*exp(-c*pmin(x,y))+abs(pmax(x,y)-pmin(x,y))*exp(-c*abs(pmax(x,y)-pmin(x,y))))
+
1.5*c^(-4)*(-exp(-c*pmax(x,y))-exp(-c*pmin(x,y))+1+exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
)

sum.c.c  <- sum.c.c + sum(diag(Vi.inv %*% derVc %*% Vi.inv %*% derVc))

}#for (i in idlist)


#### 1) c = alpha, 2) d = tausq

sum.c.d <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVc <- outer(Timei, Timei, function(x,y) b * (
0.5*c^(-3)*(-4*pmin(x,y)-pmax(x,y)*exp(-c*pmax(x,y))-pmin(x,y)*exp(-c*pmin(x,y))+abs(pmax(x,y)-pmin(x,y))*exp(-c*abs(pmax(x,y)-pmin(x,y))))
+
1.5*c^(-4)*(-exp(-c*pmax(x,y))-exp(-c*pmin(x,y))+1+exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
)
#derVd <- Ii
#sum.c.d  <- sum.c.d + sum(diag(Vi.inv %*% derVc %*% Vi.inv %*% derVd))
sum.c.d  <- sum.c.d + sum(diag(Vi.inv %*% derVc %*% Vi.inv))

}#for (i in idlist)


#### 1) d = log(tausq), 2) d = log(tausq)

sum.d.d <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

#derVd <- Ii
#sum.d.d  <- sum.d.d + sum(diag(Vi.inv %*% derVd %*% Vi.inv %*% derVd))
sum.d.d  <- sum.d.d + sum(diag(Vi.inv %*% Vi.inv))

}#for (i in idlist)

expected.hessian <- -0.5 * matrix(c(sum.a.a, sum.a.b, sum.a.c, sum.a.d,
                                    sum.a.b, sum.b.b, sum.b.c, sum.b.d,
                                    sum.a.c, sum.b.c, sum.c.c, sum.c.d,
                                    sum.a.d, sum.b.d, sum.c.d, sum.d.d), 
                                  ncol = 4, byrow = T)

theta.new <- as.numeric(theta.old - ginv(expected.hessian) %*% score.theta)

if(silent == FALSE){

cat("iteration = ", iter, "\n")
cat("theta.old = ", theta.old, "\n")
cat("beta.hat  = ", b.hat, "\n")
cat("theta.new = ", theta.new, "\n")
cat("sqrt.diff=", sqrt((theta.old - theta.new) %*% (theta.old - theta.new)), "\n")      
print("-----------------------------")

}

iter <- iter + 1

}#while

#proc.time() - ptm

theta.old <- theta.new

a <- theta.old[1]  ## omegasq
b <- theta.old[2]  ## sigmasq
c <- theta.old[3]  ## alpha
d <- theta.old[4]  ## tausq

### betahat

sum.left.beta <- sum.right.beta <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv    <- solve(Vi)
Xi.transp <- t(Xi)

sum.left.beta  <- sum.left.beta + Xi.transp %*% Vi.inv %*% Xi  
sum.right.beta <- sum.right.beta + Xi.transp %*% Vi.inv %*% Yi

}#for(i in idlist)

b.hat    <- solve(sum.left.beta) %*% sum.right.beta
b.varcov <- solve(sum.left.beta)

######### INFORMATION MATRIX

#### 1) a = omegasq, 2) a = omegasq

sum.a.a <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVa <- Ji

sum.a.a  <- sum.a.a + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVa))

}#for (i in idlist)


#### 1) a = omegasq, 2) b = sigmasq

sum.a.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVa <- Ji
derVb <- outer(Timei, Timei, function(x,y) 
0.5*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)

sum.a.b  <- sum.a.b + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVb))

}#for (i in idlist)


#### 1) a = omegasq, 2) c = alpha

sum.a.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVa <- Ji
derVc <- outer(Timei, Timei, function(x,y) b * (
0.5*c^(-3)*(-4*pmin(x,y)-pmax(x,y)*exp(-c*pmax(x,y))-pmin(x,y)*exp(-c*pmin(x,y))+abs(pmax(x,y)-pmin(x,y))*exp(-c*abs(pmax(x,y)-pmin(x,y))))
+
1.5*c^(-4)*(-exp(-c*pmax(x,y))-exp(-c*pmin(x,y))+1+exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
)

sum.a.c  <- sum.a.c + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVc))

}#for (i in idlist)


#### 1) a = omegasq, 2) d = tausq

sum.a.d <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVa <- Ji
#derVd <- Ii
#sum.a.d  <- sum.a.d + sum(diag(Vi.inv %*% derVa %*% Vi.inv %*% derVd))
sum.a.d  <- sum.a.d + sum(diag(Vi.inv %*% derVa %*% Vi.inv))

}#for (i in idlist)


#### 1) b = sigmasq, 2) b = sigmasq

sum.b.b <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVb <- outer(Timei, Timei, function(x,y) 
0.5*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)

sum.b.b  <- sum.b.b + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVb))

}#for (i in idlist)


#### 1) b = sigmasq, 2) c = alpha

sum.b.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVb <- outer(Timei, Timei, function(x,y) 
0.5*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
derVc <- outer(Timei, Timei, function(x,y) b * (
0.5*c^(-3)*(-4*pmin(x,y)-pmax(x,y)*exp(-c*pmax(x,y))-pmin(x,y)*exp(-c*pmin(x,y))+abs(pmax(x,y)-pmin(x,y))*exp(-c*abs(pmax(x,y)-pmin(x,y))))
+
1.5*c^(-4)*(-exp(-c*pmax(x,y))-exp(-c*pmin(x,y))+1+exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
)

sum.b.c  <- sum.b.c + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVc))

}#for (i in idlist)


#### 1) b = sigmasq, 2) d = tausq

sum.b.d<- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVb <- outer(Timei, Timei, function(x,y) 
0.5*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
#derVd <- Ii
#sum.b.d  <- sum.b.d + sum(diag(Vi.inv %*% derVb %*% Vi.inv %*% derVd))
sum.b.d  <- sum.b.d + sum(diag(Vi.inv %*% derVb %*% Vi.inv))

}#for (i in idlist)


#### 1) c = alpha, 2) c = alpha

sum.c.c <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVc <- outer(Timei, Timei, function(x,y) b * (
0.5*c^(-3)*(-4*pmin(x,y)-pmax(x,y)*exp(-c*pmax(x,y))-pmin(x,y)*exp(-c*pmin(x,y))+abs(pmax(x,y)-pmin(x,y))*exp(-c*abs(pmax(x,y)-pmin(x,y))))
+
1.5*c^(-4)*(-exp(-c*pmax(x,y))-exp(-c*pmin(x,y))+1+exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
)

sum.c.c  <- sum.c.c + sum(diag(Vi.inv %*% derVc %*% Vi.inv %*% derVc))

}#for (i in idlist)


#### 1) c = alpha, 2) d = tausq

sum.c.d <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

derVc <- outer(Timei, Timei, function(x,y) b * (
0.5*c^(-3)*(-4*pmin(x,y)-pmax(x,y)*exp(-c*pmax(x,y))-pmin(x,y)*exp(-c*pmin(x,y))+abs(pmax(x,y)-pmin(x,y))*exp(-c*abs(pmax(x,y)-pmin(x,y))))
+
1.5*c^(-4)*(-exp(-c*pmax(x,y))-exp(-c*pmin(x,y))+1+exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
)
#derVd <- Ii
#sum.c.d  <- sum.c.d + sum(diag(Vi.inv %*% derVc %*% Vi.inv %*% derVd))
sum.c.d  <- sum.c.d + sum(diag(Vi.inv %*% derVc %*% Vi.inv))

}#for (i in idlist)


#### 1) d = log(tausq), 2) d = log(tausq)

sum.d.d <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)

#derVd <- Ii
#sum.d.d  <- sum.d.d + sum(diag(Vi.inv %*% derVd %*% Vi.inv %*% derVd))
sum.d.d  <- sum.d.d + sum(diag(Vi.inv %*% Vi.inv))

}#for (i in idlist)

expected.hessian <- -0.5 * matrix(c(sum.a.a, sum.a.b, sum.a.c, sum.a.d,
                                    sum.a.b, sum.b.b, sum.b.c, sum.b.d,
                                    sum.a.c, sum.b.c, sum.c.c, sum.c.d,
                                    sum.a.d, sum.b.d, sum.c.d, sum.d.d), 
                                  ncol = 4, byrow = T)

sd.theta <- sqrt(diag(ginv(-expected.hessian)))

## loglik

sum.loglik <- 0

for (i in 1 : length(idlist)){

Timei <- Time[[i]]
ni    <- nobs[i]
Xi    <- as.matrix(DM[[i]], nrow = ni)
Yi    <- as.matrix(YM[[i]], nrow = ni)

Ji <- matrix(1, ni, ni)
Ri <- outer(Timei, Timei, function(x,y) 
0.5*b*c^(-3)*(2*c*pmin(x,y)+exp(-c*pmax(x,y))+exp(-c*pmin(x,y))-1-exp(-c*abs(pmax(x,y)-pmin(x,y))))
)
Ii <- diag(ni)

Vi <- a * Ji + Ri + d * Ii

Vi.inv   <- solve(Vi)
ri     <- Yi - Xi %*% b.hat

#sum.loglik  <- sum.loglik + log(det(Vi)) + t(ri) %*% Vi.inv %*% ri
sum.loglik  <- sum.loglik + as.numeric(determinant(Vi, logarithm=TRUE)$modulus) +  t(ri) %*% Vi.inv %*% ri

}#for (i in idlist)

max.loglik <- as.numeric(- 0.5 * nrow(YM) * log(2 * pi) - 0.5 * sum.loglik)


result <- cbind(c(c(b.hat), theta.new) , c(sqrt(diag(b.varcov)), sd.theta),
                c(c(b.hat)/sqrt(diag(b.varcov)), NA, NA, NA, NA))
result <- cbind(result, 2 * (1 - pnorm(abs(result[, 3]))))
colnames(result) <- c("Estimate", "Standard error", "Z-estimate", "p-value")
rownames(result)[(nrow(result) - 3) : nrow(result)] <- c("omegasq", "kappasq", "nu", "tausq")

score.theta <- matrix(score.theta, nrow = 1)
colnames(score.theta) <- c("omegasq", "kappasq", "nu", "tausq")

output           <- list()
output$title     <- "Mixed effects model with random intercept and integrated Ornstein-Uhlenbeck Process"
output$date      <- date()
output$estimates <- result
output$maxloglik <- max.loglik 
output$score     <- score.theta
output
}

output

}
