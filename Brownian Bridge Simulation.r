#Brownian Bridge Simulation
#set.seed(2)
BB <- function(a, b, t, T, m){
  z <- rnorm(m)
  Wt <- c(a)
  dt <- sort(c(t, runif((m-2), min = t, max = T), T))
  Wt_m = b
  #Wt_mean <- c()
  for(j in 2 : (m-1)){
    Wt_j = Wt[j-1] + (Wt_m - Wt[j-1]) *(dt[j] - dt[j-1])/(T - dt[j-1]) + sqrt((T - dt[j])*(dt[j] - dt[j-1])/(T - dt[j-1])) * z[j]
    Wt <- c(Wt, Wt_j)
  }
  Wt = c(Wt, Wt_m)
  plot(x = dt, y = Wt, type = "b")
  return (Wt)
}

t = 0
T = 1
m = 100000
a = rnorm(1) * sqrt(t)
b = rnorm(1) * sqrt(T)
bb_sim = BB(a, b, t, T, m)
