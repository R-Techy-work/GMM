
X <- iris[,1:4]
y <- iris[,5]
y_col <- c('#7DB0DD', '#86B875', '#E495A5')
pdf('dat_iris.pdf')
pairs(X, lower.panel = NULL, col = y_col[y])
par(xpd = T)
legend(x = 0.1, y = 0.4, legend = as.character(levels(y)), fill = y_col)
dev.off()

library(reshape)   #cast()
wd <- getwd()
# Uses EM algorithm with multivariate normal
# distribution to estimate cluster probability
mvnorm.cov.inv <- function(Sigma) {
  # Eigendecomposition of covariance matrix
  E <- eigen(Sigma)
  Lambda.inv <- diag(E$values^-1)   # diagonal matrix with inverse of eigenvalues
  Q <- E$vectors                    # eigenvectors
  return(Q %*% Lambda.inv %*% t(Q))
}
#multivariate Gaussian pdf
mvn.pdf.i <- function(xi, mu, Sigma)
  1/sqrt( (2*pi)^length(xi) * det(Sigma) ) * 
  exp(-(1/2) * t(xi - mu) %*% mvnorm.cov.inv(Sigma) %*% (xi - mu)  )
mvn.pdf <- function(X, mu, Sigma)
  apply(X, 1, function(xi) mvn.pdf.i(as.numeric(xi), mu, Sigma))
gmm.fromscratch <- function(X, k, plot = F){
  p <- ncol(X)  # number of parameters
  n <- nrow(X)  # number of observations
  Delta <- 1; iter <- 0; itermax <- 30
  class_col <- c('#7DB0DD', '#86B875', '#E495A5')
  while(Delta > 1e-4 && iter <= itermax){
    # initiation
    if(iter == 0){
      km.init <- km.fromscratch(X, k)
      mu <- km.init$centroid; mu_mem <- mu
      w <- sapply(1:k, function(i) length(which(km.init$cluster == i)))
      w <- w/sum(w)
      cov <- array(dim = c(p, p, k))
      for(i in 1:p) for(j in 1:p) for(c in 1:k) cov[i, j, c] <- 
        1/n * sum((X[km.init$cluster == c, i] - mu[c, i]) *
                    (X[km.init$cluster == c, j] - mu[c, j]))
    }
    
    # E-step
    mvn.c <- sapply(1:k, function(c) mvn.pdf(X, mu[c,], cov[,, c]))
    r_ic <- t(w*t(mvn.c)) / rowSums(t(w*t(mvn.c)))
    # M-step
    n_c <- colSums(r_ic)
    w <- n_c/sum(n_c)
    mu <- t(sapply(1:k, function(c) 1/n_c[c] * colSums(r_ic[, c] *
                                                         X)))
    for(i in 1:p) for(j in 1:p) for(c in 1:k) cov[i, j, c] <- 
      1/n_c[c] * sum(r_ic[, c] * (X[, i] - mu[c, i]) * r_ic[, c] *
                       (X[, j] - mu[c, j]))
    cluster <- apply(r_ic, 1, which.max)
    if(plot){
      i <- 1
      idx <- matrix(rep(seq(p), p), ncol = p, nrow = p)
      idx <- idx[lower.tri(idx)]
      idy <- matrix(rep(seq(p), each=p), ncol = p, nrow = p)
      idy <- idy[lower.tri(idy)]
      
      if(iter < 10) iter4plot <- paste0('0', iter) else iter4plot <- iter
      
      png(paste0(wd, '/figs_gmm/iter', iter4plot, '.png'))
      pairs(rbind(X, mu), lower.panel = NULL, asp = 1, 
            col = c(class_col[cluster], rep('black', k)), main =
              paste0('iter=',iter), panel=function(x, y, ...) {
                points(x, y, col = c(class_col[cluster], rep('black', k)))
                xi <- seq(min(X[, idx[i]])-1, max(X[, idx[i]])+1, 0.1)
                yi <- seq(min(X[, idy[i]])-1, max(X[, idy[i]])+1, 0.1)
                grid <- expand.grid(xi = xi, yi = yi)
                grid['z'] <- mvn.pdf(grid, mu[1,c(idx[i],idy[i])],
                                     cov[c(idx[i],idy[i]),c(idx[i],idy[i]), 1])
                z <- cast(grid, xi ~ yi)
                contour(xi, yi, as.matrix(z[,-1]), 
                        levels = c(.1, .5, .9), col = class_col[1], 
                        add = T, lty = 'solid', labels = '')
                grid <- expand.grid(xi = xi, yi = yi)
                grid['z'] <- mvn.pdf(grid, mu[2,c(idx[i],idy[i])], 
                                     cov[c(idx[i],idy[i]),c(idx[i],idy[i]), 2])
                z <- cast(grid, xi ~ yi)
                contour(xi, yi, as.matrix(z[,-1]), 
                        levels = c(.1, .5, .9), col = class_col[2], 
                        add = T, lty = 'solid', labels = '')
                grid <- expand.grid(xi = xi, yi = yi)
                grid['z'] <- mvn.pdf(grid, mu[3,c(idx[i],idy[i])], 
                                     cov[c(idx[i],idy[i]),c(idx[i],idy[i]), 3])
                z <- cast(grid, xi ~ yi)
                contour(xi, yi, as.matrix(z[,-1]), 
                        levels = c(.1, .5, .9), col = class_col[3], 
                        add = T, lty = 'solid', labels = '')
                i <<- i+1
              })
      dev.off()
    }
    
    Delta <- sum((mu - mu_mem)^2)
    iter <- iter + 1; mu_mem <- mu
  }
  return(list(softcluster = r_ic, cluster = cluster))
}
gmm <- gmm.fromscratch(X, 3, plot = T)
table(y, gmm$cluster)
library(magick)
list.files(path = paste0(wd, "/figs_gmm/"), pattern = "*.png", full.names = T) %>% 
  image_read() %>%
  image_join() %>%
  image_animate(fps=1) %>%
  image_write("fig_gmm_anim.gif")