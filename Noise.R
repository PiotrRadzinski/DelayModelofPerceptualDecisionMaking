sig <<- 0.3
tau <<- 1.4 # assume tau is less than 2
tmax <<- 12 # maximal simulation time
N <<- 1000 # number of points in 1 sec
stst <<- 0.411466 # steady state (starting)
b <<- 1/20

hill <- function(l)  return(l^2/(1 + l^2))

#indicator (for 0.5sec)
ind <- function(l){
     if(l < 2.5 * N) return(1)
     else return(0)
}

deterministic <- function(){
     r <<- numeric(tmax * N) #r1
     s <<- numeric(tmax * N) #r2
     for(i in 1:(2*N)){
          r[i] <<- stst
          s[i] <<- stst
     }
     for(j in (2*N + 1):(N * tmax)){
          r[j] <<- r[j-1] + (0.4 + sig * ind(j) - r[j - tau * N] + s[j-1] * hill(s[j-1]*r[j-1]))/N
          s[j] <<- s[j-1] + (0.4 - s[j - tau * N] + r[j-1] * hill(s[j-1]*r[j-1]))/N
     }
}
deterministic()

R <- numeric(tmax * N) #r1 noised
S <- numeric(tmax * N) #r2 noised
for(i in 1:(2*N)){
     R[i] <- stst
     S[i] <- stst
}

stochastic <- function(){
     for(j in (2*N + 1):(N * tmax)){
          R[j] <<- R[j-1] + (0.4 + sig * ind(j) - R[j - tau * N] + S[j-1] * hill(S[j-1]*R[j-1]))/N + b*rnorm(1,0,(1/N)^(1/2))
          S[j] <<- S[j-1] + (0.4 - S[j - tau * N] + R[j-1] * hill(S[j-1]*R[j-1]))/N + b*rnorm(1,0,(1/N)^(1/2))
     }
}

stochastic()
plot(seq(0, tmax-2, length.out = (N*tmax)-(2*N)+1), R[(2*N):(N*tmax)], type="l", col="blue", ylim = c(0,1))
lines(seq(0, tmax-2, length.out = (N*tmax)-(2*N)+1), S[(2*N):(N*tmax)], type="l", col="orange", ylim = c(0,1))
lines(seq(0, tmax-2, length.out = (N*tmax)-(2*N)+1), r[(2*N):(N*tmax)], type="l", col="blue", ylim = c(0,1))
lines(seq(0, tmax-2, length.out = (N*tmax)-(2*N)+1), s[(2*N):(N*tmax)], type="l", col="orange", ylim = c(0,1))

# Calculating difference between stochastic and deterministic solutions
Rr <- numeric(N * tmax) # abs of difference r1
Sr <- numeric(N * tmax) # abs of difference r2
for(k in 1:(N * tmax)){
     Rr[k] <- abs(R[k] - r[k])
     Sr[k] <- abs(S[k] - s[k])
#     exR <- exR + abs(R[k] - r[k])
#     exS <- exS + abs(S[k] - s[k])
}

plot(seq(0, tmax-2, length.out = (N*tmax)-(2*N)+1), Rr[(2*N):(N*tmax)], type="l", col="blue", ylim = c(0,1))
lines(seq(0, tmax-2, length.out = (N*tmax)-(2*N)+1), Sr[(2*N):(N*tmax)], type="l", col="orange", ylim = c(0,1))
# exR/(N * tmax)
# exS/(N * tmax)


# MEAN DIFFERENCE
deterministic()
nb <- 1 # number of simulation for single b
ns <- 120 # number of diffrent b 
x <- numeric(ns)
y <- numeric(ns)
tx <- numeric(nb) # counts means for r1 
ty <- numeric(nb) # counts means for r2
exR <- 0 # mean difference r1
exS <- 0 # mean difference r2
for(i in 1:ns){
     b <<- i/N
     for(j in 1:nb){
          stochastic()
          exR <- 0 
          exS <- 0 
          for(k in (2*N+1):(N * tmax)){
               exR <- exR + abs(R[k] - r[k])
               exS <- exS + abs(S[k] - s[k])
          }
          tx[j] <- exR/(N * tmax)    
          ty[j] <- exS/(N * tmax)
     }
     x[i] <- sum(tx)/nb
     y[i] <- sum(ty)/nb
     print(i)
}

z <- (x + y)/2
plot(seq(0, ns/N, length.out = ns), z[1:ns], type="l", col="blue", ylim = c(0,0.5), xlab = "b", ylab = "mean difference")

# SUCCESS RATE
integral <- numeric(tmax * N) # integral of difference between r1 and r2
stochasticWI <- function(){ # calculating SDDE with integral of difference
     integral[2*N] <<- 1/2
     for(j in (2*N + 1):(N * tmax)){
          R[j] <<- R[j-1] + (0.4 + sig * ind(j) - R[j - tau * N] + S[j-1] * hill(S[j-1]*R[j-1]))/N + b*rnorm(1,0,(1/N)^(1/2))
          S[j] <<- S[j-1] + (0.4 - S[j - tau * N] + R[j-1] * hill(S[j-1]*R[j-1]))/N + b*rnorm(1,0,(1/N)^(1/2))
          integral[j] <<- integral[j-1] + R[j-1] - S[j-1]
     }
}

nb <- 10000 # number of simulation for single b
ns <- 120 # number of diffrent b 
sc <- numeric(nb) # counts successes for single b
msc <- numeric(ns) # counts mean success rate for each b
mm <- 0 # moment when we check which decision would be maked
for(i in 1:ns){
     b <<- i/N
     sc <- numeric(nb)
     for(j in 1:nb){
          stochasticWI()
          mm <- round(runif(1,7,9), digits = 3) # for N = 1000. For diffrent N it is needed to change digits
          if ((1/(1+exp(-integral[mm*N]))) > 1/2) sc[j] <- 1
          print(i + j/nb)
     }
     msc[i] <- mean(sc)
}

# write.table(msc[1:3], file="srb_1_3_1000a.csv")
# p1 <- unlist(read.table("srb_1_3_10000a.csv"))
# scf <- c(p1,...)
# plot(seq(0, ns/N, length.out = ns), scf[1:ns], type="l", col="blue", ylim = c(0.5,1), xlab = "b", ylab = "success rate", lwd=3)
