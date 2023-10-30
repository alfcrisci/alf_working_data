
# from https://cran.r-project.org/web/packages/moments/index.html
# Authors:	Lukasz Komsta, Frederick Novomestky
# Licences


"agostino.test" <-
function (x, alternative=c("two.sided","less","greater"))
{
     DNAME <- deparse(substitute(x))
     x <- sort(x[complete.cases(x)])
     n <- length(x)
s <- match.arg(alternative)
alter <- switch(s, two.sided=0, less=1, greater=2)
     if ((n < 8 || n > 46340))
         stop("sample size must be between 8 and 46340")
s3 <- (sum((x-mean(x))^3)/n)/(sum((x-mean(x))^2)/n)^(3/2)
y <- s3*sqrt((n+1)*(n+3)/(6*(n-2)))
b2 <- 3*(n*n+27*n-70)*(n+1)*(n+3)/((n-2)*(n+5)*(n+7)*(n+9))
w <- sqrt(-1+sqrt(2*(b2-1)));
d <- 1/sqrt(log(w));
a <- sqrt(2/(w*w-1));
z <- d*log(y/a+sqrt((y/a)^2+1));
     pval <- pnorm(z, lower.tail = FALSE)
if (alter == 0) {
pval <- 2*pval
if (pval > 1) pval<-2-pval
alt <- "data have a skewness"
}
else if (alter == 1)
{
alt <- "data have positive skewness"
}
else
{
pval <- 1-pval
alt <- "data have negative skewness"
}
     RVAL <- list(statistic = c(skew = s3, z = z), p.value = pval,
alternative = alt, method = "D'Agostino skewness test",
         data.name = DNAME)
     class(RVAL) <- "htest"
     return(RVAL)
}

all.cumulants <- function(mu.raw)
{
    mu.central <- raw2central( mu.raw )
    if ( is.vector( mu.raw ) ) {
        order.max <- length( mu.raw ) - 1
        kappa <- rep( 0, order.max + 1 )
        kappa[2] <- mu.central[2]
        kappa[3] <- mu.central[3]
        n <- 3
        while( n <= order.max ) {
            np1 <- n + 1
            nm1 <- n - 1
            total <- mu.raw[np1]
            for ( k in 1:nm1 ) {
                km1 <- k - 1
                total <- total - choose(nm1,km1 ) * kappa[k+1] * mu.raw[n-k+1]
            }
            kappa[np1] <- total
            n <- n + 1
        }
        return( kappa )
    }
    else if ( is.matrix( mu.raw ) ) {
        order.max <- nrow( mu.raw ) - 1
        nvar <- ncol( mu.raw )
        kappa <- matrix( nrow=order.max+1,ncol=nvar )
        kappa[1,] <- rep( 0, nvar )
        kappa[2,] <- mu.central[2,]
        kappa[3,] <- mu.central[3,]
        for ( j in 1:nvar ) {
            n <- 3
            while ( n <= order.max ) {
                np1 <- n + 1
                nm1 <- n - 1
                total <- mu.raw[np1,j]
                for ( k in 1:nm1 ) {
                    km1 <- k - 1
                    total <- total -  choose(nm1,km1) * kappa[k+1,j] * mu.raw[n-k+1,j]
                }
                kappa[np1,j] <- total
                n <- n + 1
            }
        }
        return( kappa )
    }
    else if ( is.data.frame( mu.raw ) ) {
        order.max <- nrow( mu.raw ) - 1
        nvar <- ncol( mu.raw )
        kappa <- data.frame( matrix( nrow=order.max+1,ncol=nvar ) )
        kappa[1,] <- rep( 0, nvar )
        kappa[2,] <- mu.central[2,]
        kappa[3,] <- mu.central[3,]
        for ( j in 1:nvar ) {
            n <- 3
            while ( n <= order.max ) {
                np1 <- n + 1
                nm1 <- n - 1
                total <- mu.raw[np1,j]
                for ( k in 1:nm1 ) {
                    km1 <- k - 1
                    total <- total -  choose(nm1,km1) * kappa[k+1,j] * mu.raw[n-k+1,j]
                }
                kappa[np1,j] <- total
                n <- n + 1
            }
        }
        return( kappa )
    }
    else
        stop( "argument mu.raw is not a vector, matrix or data frame" )
    return( NULL )
}
all.moments <- function( x, order.max=2, central=FALSE, absolute=FALSE, na.rm=FALSE )
{
    if( order.max < 2 )
        stop( "maximum order should be at least 2" )
    if ( is.matrix( x ) ) {
        n <- ncol( x )
        mu <- matrix( nrow=order.max+1, ncol=n )
        for ( order in 0:order.max )
            mu[order+1,] <- moment( x, order, central, absolute, na.rm )
        return( mu )
    }
    else if ( is.vector( x ) ) {
        mu <- rep( 0, order.max+1 )
        for ( order in 0:order.max )
            mu[order+1] <- moment( x, order, central, absolute, na.rm )
        return( mu )
    }
    else if ( is.data.frame( x ) ) {
        n <- ncol( x )
        mu <- matrix( nrow=order.max+1, ncol=n )
        for ( order in 0:order.max )
            mu[order+1,] <- moment( x, order, central, absolute, na.rm )
        return( mu )
    }
    else
        return( all.moments( as.vector(x), order.max, central, absolute, na.rm ) )
}
"anscombe.test" <-
function (x, alternative=c("two.sided","less","greater"))
{
     DNAME <- deparse(substitute(x))
     x <- sort(x[complete.cases(x)])
     n <- length(x)
s <- match.arg(alternative)
alter <- switch(s, two.sided=0, less=1, greater=2)
b <- n*sum( (x-mean(x))^4 )/(sum( (x-mean(x))^2 )^2);
eb2 <- 3*(n-1)/(n+1);
vb2 <- 24*n*(n-2)*(n-3)/ ((n+1)^2*(n+3)*(n+5));
m3 <- (6*(n^2-5*n+2)/((n+7)*(n+9)))*sqrt((6*(n+3)*(n+5))/(n*(n-2)*(n-3)));
a <- 6 + (8/m3) * (2/m3 + sqrt(1 + 4/m3^2));
xx <- (b-eb2)/sqrt(vb2);
z <- ( 1-2/(9*a)-( (1-2/a) / (1+xx*sqrt(2/(a-4))) )^(1/3))/ sqrt(2/(9*a));
     pval <- pnorm(z, lower.tail = FALSE)
if (alter == 0) {
pval <- 2*pval
if (pval > 1) pval<-2-pval
alt <- "kurtosis is not equal to 3"
}
else if (alter == 1)
{
alt <- "kurtosis is greater than 3"
}
else
{
pval <- 1-pval
alt <- "kurtosis is lower than 3"
}
     RVAL <- list(statistic = c(kurt = b, z = z), p.value = pval,
alternative = alt, method = "Anscombe-Glynn kurtosis test",
         data.name = DNAME)
     class(RVAL) <- "htest"
     return(RVAL)
}

"bonett.test" <-
function (x, alternative=c("two.sided","less","greater"))
{
     DNAME <- deparse(substitute(x))
     x <- sort(x[complete.cases(x)])
     n <- length(x)
s <- match.arg(alternative)
alter <- switch(s, two.sided=0, less=1, greater=2)
rho <- sqrt(sum((x-mean(x))^2)/n);
tau <- sum(abs(x-mean(x)))/n;
omega <- 13.29*(log(rho)-log(tau));
z <- sqrt(n+2)*(omega-3)/3.54;
     pval <- pnorm(z, lower.tail = FALSE)
if (alter == 0) {
pval <- 2*pval
if (pval > 1) pval<-2-pval
alt <- "kurtosis is not equal to sqrt(2/pi)"
}
else if (alter == 1)
{
alt <- "kurtosis is greater than sqrt(2/pi)"
}
else
{
pval <- 1-pval
alt <- "kurtosis is lower than sqrt(2/pi)"
}
     RVAL <- list(statistic = c(tau = tau, z = z), alternative = alt,
p.value = pval, method = "Bonett-Seier test for Geary kurtosis",
         data.name = DNAME)
     class(RVAL) <- "htest"
     return(RVAL)
}

central2raw <- function( mu.central, eta )
{
    if ( is.vector( mu.central ) ) {
        if ( length( eta ) > 1 )
            stop( "argument eta has too many values" )
        np1 <- length( mu.central )
        n <- np1 - 1
        mu.raw <- rep( 0, np1 )
        mu.raw[1] <- 1
        j <- 2
        for ( k in 1:n ) {
            total <- 0
            for ( i in 0:k ) {
                total <- total + choose(k,i) * ( ( eta )^(k-i) ) * mu.central[i+1]
            }
            mu.raw[j] <- total
            j <- j + 1
        }
        return( mu.raw )
    }
    mu.raw <- NULL
    if ( is.matrix( mu.central ) ) {
        np1 <- nrow( mu.central )
        nrv <- ncol( mu.central )
        mu.raw <- matrix( nrow=np1, ncol=nrv )
    }
    if ( is.data.frame( mu.central ) ) {
        np1 <- nrow( mu.central )
        nrv <- ncol( mu.central )
        mu.raw <- data.frame( matrix( nrow=np1, ncol=nrv ) )
    }
    if ( is.null( mu.raw ) )
        stop( "argument mu.central is not a vector, matrix or data frame" )
    if ( length( eta ) != nrv )
        stop( "argument eta does not match the argument mu.central" )
    n <- np1 - 1
    for ( m in 1:nrv ) {
        mu.raw[1,m] <- 1
        eta.m <- eta[m]
        j <- 2
        for ( k in 1:n ) {
            total <- 0
            for ( i in 0:k ) {
                total <- total + choose(k,i) * ( ( eta.m )^(k-i) ) * mu.central[i+1,m]
            }
            mu.raw[j,m] <- total
            j <- j + 1
        }
    }
    return( mu.raw )

}

"geary" <-
function (x, na.rm = FALSE)
{
    if (is.matrix(x))
        apply(x, 2, geary, na.rm = na.rm)
    else if (is.vector(x)) {
        if (na.rm) x <- x[!is.na(x)]
        n <- length(x)
        rho <- sqrt(sum((x-mean(x))^2)/n);
        tau <- sum(abs(x-mean(x)))/n;
        tau/rho
        }
    else if (is.data.frame(x))
        sapply(x, geary, na.rm = na.rm)
    else geary(as.vector(x), na.rm = na.rm)
}

jarque.test <- function(x)
{
    if ( !is.vector( x ) )
        stop( "argument x is not a vector" )
    if ( !is.numeric( x ) )
        stop( "argument x is not numeric" )
    DNAME <- deparse( substitute( x ) )
    n <- length(x)
    ALTERNATIVE <- "greater"
    METHOD <- "Jarque-Bera Normality Test"
    K <- kurtosis( x )
    S <- skewness( x )
    JB  <- ( n / 6 ) * ( S^2 + 0.25 * ( ( K - 3 )^2 ) )
    pval <- 1 - pchisq( JB, df=2 )
    JBVAL <- list( statistic=c(JB=JB), p.value=pval, alternative=ALTERNATIVE, method=METHOD,
                   data.name=DNAME )
    class( JBVAL ) <- "htest"
    return( JBVAL )
}

"kurtosis" <-
function (x, na.rm = FALSE)
{
    if (is.matrix(x))
        apply(x, 2, kurtosis, na.rm = na.rm)
    else if (is.vector(x)) {
        if (na.rm) x <- x[!is.na(x)]
        n <- length(x)
        n*sum( (x-mean(x))^4 )/(sum( (x-mean(x))^2 )^2)
        }
    else if (is.data.frame(x))
        sapply(x, kurtosis, na.rm = na.rm)
    else kurtosis(as.vector(x), na.rm = na.rm)
}

"moment" <-
function(x, order = 1, central = FALSE, absolute = FALSE, na.rm = FALSE)
{
    if (is.matrix(x))
        apply(x, 2, moment, order = order, central = central, absolute = absolute, na.rm = na.rm)
    else if (is.vector(x)) {
          if (na.rm) x = x[!is.na(x)] ;
          if (central) x = x - mean(x)
          if (absolute) x = abs(x)
          sum(x^order)/length(x)
                }
    else if (is.data.frame(x))
        sapply(x, moment, order = order, central = central, absolute = absolute, na.rm = na.rm)
    else moment(as.vector(x), order = order, central = central, absolute = absolute, na.rm = na.rm)
}

raw2central <- function( mu.raw )
{
    if ( is.matrix( mu.raw ) )
        return( apply( mu.raw, 2, raw2central ) )
    else if ( is.vector( mu.raw ) ) {
        np1 <- length( mu.raw )
        n <- np1 - 1
        mu.central <- rep( 0, np1 )
        mu.central[1] <- 1
        eta <- mu.raw[2]
        j <- 2
        for ( k in 1:n ) {
            total <- 0
            for ( i in 0:k ) {
                total <- total + choose(k,i) * ( (- eta )^(k-i) ) * mu.raw[i+1]
            }
            mu.central[j] <- total
            j <- j + 1
        }
        return( mu.central )
    }
    else if ( is.data.frame( mu.raw ) )
        return( sapply( mu.raw, raw2central ) )
    else
        return( raw2central( as.vector( mu.raw ) ) )
}

"skewness" <-
function (x, na.rm = FALSE)
{
    if (is.matrix(x))
        apply(x, 2, skewness, na.rm = na.rm)
    else if (is.vector(x)) {
        if (na.rm) x <- x[!is.na(x)]
        n <- length(x)
     (sum((x-mean(x))^3)/n)/(sum((x-mean(x))^2)/n)^(3/2)
        }
    else if (is.data.frame(x))
        sapply(x, skewness, na.rm = na.rm)
    else skewness(as.vector(x), na.rm = na.rm)
}
