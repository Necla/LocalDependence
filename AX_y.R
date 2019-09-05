AX_y <- function (X,Y,y) {
        n <- length(X)
        s_y <- sqrt(var(Y))
        mu <- mean(Y)
        pay = 0; payda = 0;
        h <- 1.06*s_y*n^{-1/5}
        U <- y-Y
        U <- 2*(U-min(U))/(max(U)-min(U))-1
        
        for (i in 1:n){
                
                kernel <- 1-abs(U[i]) 
                a <- X[i]*kernel
                pay = pay + a
                payda = payda + kernel
        }
        return (pay/payda)
} 
