X = cac# you can change your data
n <- length(X)
##############################################################################################################################
##############################################################################################################################
#2. fit an ARMA - GARCH model to the simulated data
##############################################################################################################################
##############################################################################################################################
spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
                        mean.model = list(armaOrder = c(1, 0), include.mean = TRUE), 
                        distribution.model = "std")
fit <-ugarchfit(spec = spec, data = X)
fit
mu <- fitted(fit)# fitted hat{mu}_t (= hat{X}_t)
sig <- sigma(fit)# fitted hat{sigma}_t
## Data X_t and fitted hat{X}_t
plot(X, type = "l", xlab = 't', ylab = expression("Data"~X[t]~"and fitted values"~hat(mu)[t]))
lines(as.numeric(mu), col = adjustcolor("blue", alpha.f = 0.5))
legend("bottomright", bty = "n", lty = c(1,1),
       col = c("black", adjustcolor("blue", alpha.f = 0.5)),
       legend = c(expression(X[t]), expression(hat(mu)[t])))

resi <- as.numeric(residuals(fit))
plot(resi, type = "l", xlab = "t", ylab = expression(epsilon[t]))
##plot Zt
#plot(fit)

##############################################################################################################################
##############################################################################################################################
#3. Calculate VaR
##############################################################################################################################
##############################################################################################################################
alpha <- 0.99
VaR <- as.numeric(quantile(fit, probs = alpha))
nu. <- fit@fit$coef["shape"]
VaR. <- as.numeric(mu + sig * sqrt((nu. - 2)/nu.)*qt(alpha, df = nu.))
stopifnot(all.equal(VaR., VaR))
## => quantile(<rugarch object>, probs = alpha) provides VaR_alpha = hat{mu}_t + hat{sigma}_t * q_Z(alpha)

##############################################################################################################################
##############################################################################################################################
#4. Backtesting Via Randomness check
##############################################################################################################################
##############################################################################################################################
btest <- VaRTest(alpha, actual = X, VaR = VaR, conf.level = 0.99)
btest$expected.exceed
btest$actual.exceed
btest$uc.Decision

#independence test
len = length(X)
Ut = pt(X,df = nu.)
Nt = sort(rt(Ut, df = nu.))
J = seq(1, len, by = 1)
Nt.star = qt((J - 0.5)/length(X), df = nu.)
cor(Nt, Nt.star)
lmode = summary(lm(Nt~Nt.star))
beta = lmode$coefficients[2,1]
std.beta = lmode$coefficients[2,2]
t.value = (beta - 1)/std.beta
if(abs(t.value) > qt(0.975, df = (len - 2))){
  paste("Failed to reject H0");
}else{
  paste("Reject H0")
}

#dstn test
Ut.dstn = pt(X,df = nu.)
Nt.dstn = rt(Ut.dstn, df = nu.)
cor.max = max(acf(Nt.dstn, lag.max = 5)$acf[2:6])
cor.max

##############################################################################################################################
##############################################################################################################################
#5. predict VaR based on fitted model
##############################################################################################################################
##############################################################################################################################
m <- ceiling(n / 10) # number of steps to forecast; => roll m-1 times with frequency 1
fspec <- getspec(fit) # specification of the fitted process
setfixed(fspec) <- as.list(coef(fit))
pred <- ugarchforecast(fspec, data = X, n.ahead = 1, n.roll = m - 1, out.sample = m - 1) # predict from the fitted process

mu.predict <- fitted(pred)
sig.predict <- sigma(pred)
VaR.predict <- as.numeric(quantile(pred, probs = alpha))

## Sanity checks
stopifnot(all.equal(mu.predict, pred@forecast$seriesFor, check.attributes = FALSE),
          all.equal(sig.predict, pred@forecast$sigmaFor, check.attributes = FALSE)) # sanity check
nu. <- pred@model$fixed.pars$shape # extract (fitted) d.o.f. nu
VaR. <- as.numeric(mu.predict + sig.predict * sqrt((nu.-2)/nu.) *
                     qt(alpha, df = nu.)) # VaR_alpha computed manually
stopifnot(all.equal(VaR., VaR.predict))
##############################################################################################################################
##############################################################################################################################
#6. Simulate future trajectories of (Xt) and compute corresponding VaRs
##############################################################################################################################
##############################################################################################################################
B <- 1000
X.boot <- ugarchpath(fspec, n.sim = m, m.sim = B, rseed = 124)
## Bootstrap VaR
## Note: Each series is now an (m, B) matrix (each row is one path)
X.t.boot <- fitted(X.boot)
sig.t.boot <- sigma(X.boot)
eps.t.boot <- X.boot@path$residSim
VaR.boot <- (X.t.boot - eps.t.boot) + sig.t.boot * sqrt((nu. - 2)/nu.)*qt(alpha,df = nu.)# (m, B) matrix
## => Bootstrapped VaR_alpha computed manually
## Compute bootstrapped two-sided 99%-confidence intervals for VaR
VaR.CI <- apply(VaR.boot, 1, function(x) quantile(x, probs = c(0.005, 0.995)))

##############################################################################################################################
##############################################################################################################################
# 7. Display the results
##############################################################################################################################
##############################################################################################################################

yran <- range(X, mu, VaR, mu.predict, VaR.predict, VaR.CI)
myran <- max(abs(yran))
yran <- c(-myran, myran)
xran <- c(1, length(X) + m)

## Simulated data (X_t) and estimated conditional mean and VaR_alpha
plot(X, type = "l", xlim = xran, ylim = yran, xlab = "Time t", ylab = "",
     main = "Data, fitted AR-GARCH process, VaR, VaR predictions and VaR CIs",
     sub = paste0("Expected exceedances: ", btest$expected.exceed, "   Actual exceedances: ",
                  btest$actual.exceed, "   Test decision: ", btest$uc.Decision))
lines(as.numeric(mu), col = adjustcolor("darkblue", alpha.f = 0.5)) # hat{\mu}_t
lines(VaR, col = "darkred") # estimated VaR_alpha

## Predictions
t. <- length(X) + seq_len(m) # future time points
lines(t., mu.predict, col = "blue") # predicted process X_t (or mu_t)
lines(t., VaR.predict, col = "red") # predicted VaR_alpha
lines(t., VaR.CI[1,], col = "orange") # lower 99%-CI for VaR_alpha
lines(t., VaR.CI[2,], col = "orange") # upper 99%-CI for VaR_alpha
legend("bottom", bty = "o", lty = rep(1, 6), lwd = 1.2,
       col = c("black", adjustcolor("darkblue", alpha.f = 0.5), "blue",
               "darkred", "red", "orange"),
       cex = 0.5,
       y.intersp = 0.5,
       legend = c(expression(X[t]), expression(hat(mu)[t]),
                  expression("Predicted"~mu[t]~"(or"~X[t]*")"),
                  substitute(widehat(VaR)[a], list(a = alpha)),
                  substitute("Predicted"~VaR[a], list(a = alpha)),
                  substitute("99%-CI for"~VaR[a], list(a = alpha))))
