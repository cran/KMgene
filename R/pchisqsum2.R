pchisqsum2<-function (Q, lambda, delta = rep(0, length(lambda)), method = c("saddlepoint", 
    "integration", "liu"), acc = 1e-07) 
{
    method <- match.arg(method)
    delta <- delta[lambda > 0]
    lambda <- lambda[lambda > 0]
    if (method == "saddlepoint") {
        p = saddle(Q, lambda, delta)
        if (is.na(p)) {
            method <- "integration"
        }
        else {
            return(list(p = p, errflag = 0))
        }
    }
    if (method == "integration") {
        tmp <- davies(q = Q, lambda = lambda, delta = delta, 
            acc = acc)
        if (tmp$ifault > 0) {
            lambda <- zapsmall(lambda, digits = 2)
            delta <- delta[lambda > 0]
            lambda <- lambda[lambda > 0]
            tmp <- farebrother(q = Q, lambda = lambda, delta = delta)
            tmp$Qq <- tmp$res
        }
        return(list(p = tmp$Qq, errflag = 0))
    }
    if (method == "liu") {
        tmp <- liu(Q, lambda = lambda, delta = delta)
        return(list(p = tmp, errflag = 0))
    }
}
