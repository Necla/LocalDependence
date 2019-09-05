AY_x <- function (X,Y,x) {
        
        n <- length(Y)
        mu <- mean(X)
        s_x <- sqrt(var(X))
        pay = 0; payda = 0;
        h <- 1.06*s_x*n^{-1/5}
        U <- x-X
        U <- 2*(U-min(U))/(max(U)-min(U))-1
        
        for (i in 1:n){
                
                kernel <- 1-abs(U[i])
                a <- Y[i]*kernel
                pay = pay + a
                payda = payda + kernel
        }
        return (pay/payda)
} 