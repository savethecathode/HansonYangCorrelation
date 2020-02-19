# in R session: source("HansonYang.R) 

require(stats) # for lowess, rpois, rnorm

# Set path to R-packages:
.libPaths("/Users/LazaridisLab/Documents/Rstuff/Rpackageslibrary")
# Load library: Practical Math, erf()
library(pracma)

###### Autocorrelation #######
###### Hanson, Yang, JCP, 2008
###### Discrete Fourier Transform (DFT) method

# Shift Function for generating the correlation function
# NOTE: lag>0 shifts left
shift <- function(x, lag) {
    n <- length(x)
    xnew <- rep(NA, n)
    if (lag < 0) {
        xnew[1:abs(lag)] <- x[(n+lag+1):n]
        xnew[(abs(lag)+1):n] <- x[1:(n+lag)]
    } else if (lag > 0) {
        xnew[1:(n-lag)] <- x[(lag+1):n]
        xnew[(n-lag+1):n] <- x[1:lag]
    } else {
        xnew <- x
    }
    return(xnew)
}

# Empirical Correlation function: individual time-lag
# Pass time series, and lag index
CxxFT <- function(x, lag) {
    f <- sum(( x*shift(x,lag) ))/length(x) - mean(x)^2
    return(f)
}

# Autocorrelation function: Fourier Transform method
# Pass sequence, and number lags
Auto <- function(x, lags) {
    f = rep(NA, abs(lags)+1)	# initialize sequence, NA's place-holders
    f[1] = CxxFT(x,0)	# (note: term appears omitted in fig.2b)
    #f[1] = 0	# Set 1st val 0 (to include origin in plot-lines)

    for (i in seq(lags)) {
       f[i+1] <- CxxFT(x,i)
    }
    return(f)
}

# Calc Variance of time series
varCxxFT <- function(x) {
    N  = length(x)
    N3i= 1/N/N/N
    E  = mean(x)
    E2 = mean(x^2)
    E3 = mean(x^3)
    E4 = mean(x^4)

    # Keep formula on same line!
    y <- N3i*E4 - 4*N3i*E3*E + N3i*(N-3)*(N+1)*E2*E2 - 2*N3i*(N*N-2*N-6)*E2*E*E + N3i*(N*N-2*N-6)*E*E*E*E
    return(y)
}

# Poisson variance, pass: number points, lambda
varCxxFTpoi <- function(N, lambda) {
    lambda*lambda/N - 2*lambda*lambda/N/N + lambda/N/N/N
}

# Gaussian variance: pass mean(mu), S.D.(sigma)
varCxxFTgaus <- function(N, mu, sigma) {
    sigma^4/N - 2*sigma^4/N/N + 8*mu^2*sigma^2/N/N/N
}

# Calc Statistical Test-function for Autocorrelation: Eq.14
# NOTE: Sum omits variance (s=0 in Eq.4) index=1 in Auto()
Z_auto <- function(x, lags) {
    f <- abs(sum( Auto(x,lags)[2:length(Auto(x,lags))] )/lags)
    return(f)
}

# Critical Region: Autocorrelation
c_auto <- function(x, lags, alpha) {
    c <- sqrt(2)*erfcinv(alpha)*sqrt( varCxxFT(x)/lags )
    return(c)
}

# Statistical Test for Autocorrelation: Hanson & Yang, JCP, 2008
# alpha=0.05, 95% confidence interval: 1.96*barSIGMAxxFT
# alpha=0.10, 90% confidence interval: 1.64*barSIGMAxxFT
# alpha=0.31, 62% confidence interval: barSIGMAxxFT
TESTauto <- function(x, lags) {
    Z <- Z_auto(x, lags) 
    #c <- sqrt( varCxxFT(x)/lags )
    #c <- 1.64*sqrt( varCxxFT(x)/lags )
    c <- 1.96*sqrt( varCxxFT(x)/lags )
    test <- Z > c
    return(test)
    #return(list(test,Z,c))	# (testing: show all vars)
}

# TESTING Statistical Test for Autocorrelation: Poisson distribution
# alpha=0.05, 95% confidence interval: 1.96*barSIGMAxxFT
# alpha=0.10, 90% confidence interval: 1.64*barSIGMAxxFT
# alpha=0.31, 62% confidence interval: barSIGMAxxFT
TESTpoi <- function(x, lags, lambda) {
    N = length(x)
    Z <- Z_auto(x, lags) 
    #c <- sqrt( varCxxFTpoi(N,lambda)/lags )	# CI: alpha=0.31
    #c <- 1.64*sqrt( varCxxFTpoi(N,lambda)/lags )	# CI: alpha=0.10
    c <- 1.96*sqrt( varCxxFTpoi(N,lambda)/lags )	# CI: alpha=0.05
    #test <- Z_auto(x, lags) > 1.96*sqrt( varCxxFTpoi(lambda,N)/lags )
    test <- Z > c
    #return(list(test,Z,c))	# (testing: show all vars)
    return(test)	# TRUE=false positive
}

# Plot Autocorrelation function
plotCf <- function(f, lags) {
    y <- Auto(f,lags)
    x <- seq(0,lags)
    plot(x,y,xlim=c(0,lags),type="l")
}

###############################
###############################
###### Cross-Correlation ######
# Ref. Hanson, Yang. JPCB 2008

# Cross Correlation, (Eq.1): i-th time-lag
# Assm X, Y periodicity.
Cxy <- function(x, y, lag) {
    N = length(x)
    C = sum( ( x-mean(x) )*( shift(y,lag)-mean(y) ) )/N
    return(C)
}

# Cross Correlation function (Eq.1) in functional form
# Pass: sequences, and number lags
Cross <- function(x, y, lags) {
    f = rep(NA, abs(lags)+1)	# initialize sequence, NA's place-holders
    f[1] = Cxy(x,y,0)

    for (i in seq(lags)) {
       f[i+1] = Cxy(x,y,i)
    }
    return(f)
}

# Test Function for Cross Correlation (Eq.4)
# when X and Y both have no autocorrelation
Z_cross <- function(x, y, lags) {
    f = abs(sum( Cross(x,y,lags)[2:length(Cross(x,y,lags))] )/lags)
    return(f)
}

# Variance of cross correlation when X, Y have no autocorrelation
# assm X, Y independent.
# assm elements of {x_i}, {y_i} mutually independent.
# X and Y must be the same length
varCxy <- function(x,y) {
    N = length(x)
    var = (N-1)/N/N*var(x)*var(y)
    return(var)
}

# Critical Region for Cross Correlation (Eq.5)
# when X and Y both have no autocorrelation.
c_cross <- function(x, y, lags, alpha) {
    c = sqrt(2)*erfcinv(alpha)*sqrt( varCxy(x,y)/lags )
    return(c)
}

# Run Instasnce of Statistical Test for Cross Correlation (Eq.4)
# when X, and Y both have no autocorrelation
TESTcross <- function(x, y, lags) {
    Z = Z_cross(x,y,lags) 
    #c = sqrt( varCxy(x,y)/lags )
    #c = 1.64*sqrt( varCxy(x,y)/lags )
    c = 1.96*sqrt( varCxy(x,y)/lags )

    test = Z > c
    return(test)
    #return(list(test,Z,c))	# (testing: show all vars)
}

# Variance of Cross Correlation for time series with autocorrelation
# scale variance by longest autocorrelation decay
varCscaled <- function(x, y) {
    #m1=10		# (testing)
    #m2=10		# (testing)
    m1 = 1/exp(coef(fitAuto(Auto(x,25)))[3])
    m2 = 1/exp(coef(fitAuto(Auto(y,25)))[3])
    m  = max(m1,m2)	# find max
    N  = length(x)/m	# def Neff
    var = (N-1)/N/N*var(x)*var(y)
    return(var)
}

TESTscaled <- function(x, y, lags) {
    Z = Z_cross(x,y,lags) 
    #c = sqrt( varCscaled(x,y)/lags )
    #c = 1.64*sqrt( varCscaled(x,y)/lags )
    c = 1.96*sqrt( varCscaled(x,y)/lags )

    test = Z > c
    return(test)
    #return(list(test,Z,c))	# (testing: show all vars)
}

# 2-State Switching, pass: N, lambda (Poisson)
# Features:
# swtich If unif dist randon number is less than kf or kb.
# kf = kb = 1/20 (forward and reverse rates).
# signal-to-background ratio = 2: (chose: shift v. scalar)
# background has ave. 20 photon count: lambda.
# actual number counts generated w Poisson distribution.
# donor and acceptor traces generated separately.
donor <- function(N, lambda) {
    k = 0.05		# transition rate
    x = rpois(N,lambda) 	# initialize time series
    S = 2*mean(x)-mean(x)	# (adjust so signal-to-background ratio=2)
    s = -1		# initialize switch (set -1)
    f = 0		# initialize flip
    #set.seed(158299)

    for (i in 1:length(x)) {
      x[i] = x[i]+f*s*S		# Start: (flip=0)(switch=-1)=0
      r = runif(n=1,min=-0.0001,max=1.0001)		# run unif rand n. generator

      if (r < k) {		# 1st Transition: (flip=1)(switch=1)=1
	s=-1*s			# 2nd Trnasition: (flip=0)(switch=-1)=0
        f=f+s			# 3nd Trnasition: (flip=1)(switch=1)=1
      }
    }
    return(x)
}

acceptor <- function(N, lambda) {
    k = 0.1		# transition rate
    x = rpois(N,lambda) 	# initialize time series (Poisson)
    #x = rnorm(N,10,10) 	# initialize time series (Poisson)
    S = 2*mean(x)-mean(x)	# (adjust so signal-to-background ratio=2)
    x = x+S 		# initialize time series
    s = -1		# initialize switch (set -1)
    f = 0		# initialize flip
    #set.seed(845721)

    for (i in 1:length(x)) {
      x[i] = x[i]-f*s*S		# Start: (flip=0)(switch=-1)=0
      r = runif(n=1,min=-1.0001,max=1.0001)		# run unif rand n. generator

      if (abs(r) <= k) {		# 1st Transition: (flip=1)(switch=1)=1
	s=-1*s			# 2nd Trnasition: (flip=0)(switch=-1)=0
        f=f+s			# 3nd Trnasition: (flip=1)(switch=1)=1
      }
    }
    return(x)
}

# Calc Decay of Autocorrelation, pass: Autocorrelation function
# SSasymp(): Asym+(R0-Asym)*exp(-exp(lrc)*input)
# lrc: log of rate const
fitAuto <- function(C) {
    dat=list(x=seq(0,length(C)-1),y=C)

    # Convert to dataframe, name components in declaration
    df <- data.frame(X=dat$x, Y=dat$y)

    # Construct fit: nonlinear least-squares
    # use self-start(SS) model
    fit <- nls(Y ~ SSasymp(X, Asym, R0, lrc), data = df)

    return(fit)

    # Access coefficients in separate function...
    #coefs <- coef(fit)	# access coefs via [indices]...

    ###predict(fit, newdata=data.frame(X=x))

    #plot(df$X, df$Y, type="l", xlab="", ylab="", main=NULL)
}

fitLine <- function(C) {
    dat=list(x=seq(0,length(C)-1),y=C)

    # Convert to dataframe, name components in declaration
    df <- data.frame(X=dat$x, Y=dat$y)

    fit <- lm(Y ~ X, data=df)

    return(fit)
}
