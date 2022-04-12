#library(DescTools)

KendallTauB_T<-function (x, y = NULL, conf.level = NA, ...) 
{
    if (!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
    }
    else {
        dname <- deparse(substitute(x))
    }
    if (!is.null(y)) 
        tab <- table(x, y, ...)
    else tab <- as.table(x)
    x <- ConDisPairs(tab)
    n <- sum(tab)
    n0 <- n * (n - 1)/2
    ti <- rowSums(tab)
    uj <- colSums(tab)
    n1 <- sum(ti * (ti - 1)/2)
    n2 <- sum(uj * (uj - 1)/2)
    taub <- (x$C - x$D)/sqrt((n0 - n1) * (n0 - n2))
    pi <- tab/sum(tab)
    pdiff <- (x$pi.c - x$pi.d)/sum(tab)
    Pdiff <- 2 * (x$C - x$D)/sum(tab)^2
    rowsum <- rowSums(pi)
    colsum <- colSums(pi)
    rowmat <- matrix(rep(rowsum, dim(tab)[2]), ncol = dim(tab)[2])
    colmat <- matrix(rep(colsum, dim(tab)[1]), nrow = dim(tab)[1], 
        byrow = TRUE)
    delta1 <- sqrt(1 - sum(rowsum^2))
    delta2 <- sqrt(1 - sum(colsum^2))
    tauphi <- (2 * pdiff + Pdiff * colmat) * delta2 * delta1 + 
        (Pdiff * rowmat * delta2)/delta1
    sigma2 <- ((sum(pi * tauphi^2) - sum(pi * tauphi)^2)/(delta1 * 
        delta2)^4)/n
    if (is.na(conf.level)) {
        result <- taub
    }
    else {
        pr2 <- 1 - (1 - conf.level)/2
        ci <- qt(pr2,length(x)-2) * sqrt(sigma2) * c(-1, 1) + taub
        result <- c(tau_b = taub, lwr.ci = max(ci[1], -1), upr.ci = min(ci[2], 
            1))
    }
    return(result)
}

#KendallTauB_T(A,B, conf.level=0.95)
