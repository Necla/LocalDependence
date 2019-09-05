H <- function(X,Y,x, y) {
        
        rho <- cor(X,Y)
        x_bar <- mean(X)
        y_bar <- mean(Y)
        s_x <- sqrt(var(X))
        s_y <- sqrt(var(Y))
        H <- (rho + ((x_bar-AX_y(X,Y,y))*(y_bar-AY_x(X,Y,x)))/(s_x*s_y)) / (sqrt(1+(x_bar-AX_y(X,Y,y))^2/(s_x^2))*sqrt(1+(y_bar-AY_x(X,Y,x))^2/(s_y^2))) 
        H = s_x*s_y*H
        return(H)
        
}