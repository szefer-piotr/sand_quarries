results <- list()
n <- 100
for(i in 1:n) results[[i]] <- glm(y ~ x1 + x2, family="poisson")
x <- do.call(cbind, lapply(results, predict, type="response"))
DD(x)

DD <- function(x, ref.t=1, zero.rm=FALSE){
  lmb <- apply(x, 1, sum)
  D2 <- log(lmb/lmb[ref.t])
  x.p <- t(apply(x, 1, function(z)z/sum(z)))
  Pt <- x.p[ref.t,]
  if(zero.rm==FALSE){
    D1 <- -t(apply(x.p, 1, function(z)ifelse(Pt==0, 0, log(Pt/z)))) %*% Pt}else{D1 <- -t(apply(x.p, 1, function(z)ifelse(Pt==0|z==0, 0, log(Pt/z)))) %*% Pt}
  D <- D1 + D2
  data.frame(D, D1, D2)
}
