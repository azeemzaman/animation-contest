library(MASS)
library(ggplot2)
library(animation)
require(gridExtra)
mu1 <- c(1, -1)
mu2 <- c(-1, 1)
Sigma <- matrix(c(1, .4, .4, .5), byrow = T, nrow = 2)
x1 <- mvrnorm(n = 1000, mu1, Sigma)
x2 <- mvrnorm(n = 1000, mu2, Sigma)
x1 <- cbind(x1,1)
x2 <- cbind(x2,-1)
all.obs <- as.data.frame(rbind(x1,x2))
colnames(all.obs) <- c("X", "Y", "label")
all.obs$label <- as.factor(all.obs$label)

ggplot(all.obs) + geom_point(aes(x = X, y = Y, col = label))

mu.bar <- (mu1 + mu2)/2
mu.diff <- mu1 - mu2
alpha <- solve(Sigma)%*%(mu.diff)
# alpha1*x + alpha2*y = alpha*mu.bar
# so y = (alpha*mu.bar - alpha1*x)/alpha2
# y = alpha*mu.bar/alpha2 - alpha1/alpha2*x
# so slope = -alpha1/alpha2
# intercept = alpha*mu.bar/alpha2

lin.slope <- -alpha[1]/alpha[2]
intercept <- t(alpha)%*%mu.bar

ggplot(all.obs) + geom_point(aes(x = X, y = Y, col = label)) +
  geom_abline(slope = lin.slope, intercept = intercept)

# now figure out parameters for projection
# to do this, rotate so that boundary is vertical
# what is the angle formed by the line?
# it is inverse tangent of the slope
lin.angle <- atan2(1, lin.slope)
perp.angle <- lin.angle

# this function takes an angle, and the true parameters
# and returns a plot
anglePlot <- function(perp.angle, mu1, mu2, Sigma, opt.slope){
  # calculate rotation matrix for this angle
  ort.angle <- perp.angle - pi/2
  rot.mat <- matrix(c(cos(perp.angle), sin(perp.angle),
                      -sin(perp.angle), cos(perp.angle)), nrow = 2)
  ort.mat <- matrix(c(cos(ort.angle), sin(ort.angle),
                      -sin(ort.angle), cos(ort.angle)), nrow = 2)
  # using the rotation matrix we can get the parameters of the projected distribution
  # the means become
  mu1.rot <- rot.mat%*%mu1
  mu2.rot <- rot.mat%*%mu2
  # now get the rotated covariance matrix
  cov.rot <- rot.mat%*%Sigma%*%t(rot.mat)
  # get the marginal distributions
  mu1.marg <- mu1.rot[1,]
  mu2.marg <- mu2.rot[1,]
  # get marginal variance
  var.marg <- cov.rot[1,1]
  
  # get values for densities
  # set x range
  Xs <- seq(from = -3, to = 3, by = 0.01)
  d1 <- dnorm(Xs, mean = mu1.marg, sd = sqrt(var.marg))
  d2 <- dnorm(Xs, mean = mu2.marg, sd = sqrt(var.marg))
  # then the angle of the orthongal space should be theta - pi/2
  d1.proj <- t(t(rot.mat)%*%t(cbind(Xs,d1)))
  d2.proj <- t(t(rot.mat)%*%t(cbind(Xs,d2)))
  d1.df <- data.frame(X = d1.proj[,1],
                      Y = d1.proj[,2])
  d2.df <- data.frame(X = d2.proj[,1],
                      Y = d2.proj[,2])
  # need to calculate slope for separating and perp spaces
  lin.slope <- 1/tan(perp.angle)
  p <- ggplot(all.obs) + geom_point(aes(x = X, y = Y, color = label), alpha = 0.1) +
    geom_abline(slope = lin.slope, intercept = intercept) +
    geom_point(data = d1.df, aes(x = X, y = Y), 
               color = "blue", size = .5, pch = 21) + 
    geom_point(data = d2.df, aes(x = X, y = Y), 
               col = "darkseagreen", size = 0.5, pch = 21) +
    geom_abline(slope = -1/lin.slope, intercept = intercept, color = "red") +
    guides(color = FALSE) + 
    scale_color_manual(values = c("darkseagreen", "blue")) +
    coord_fixed() +
    theme_minimal() +
    geom_abline(slope = opt.slope, intercept = 0, linetype = "dashed") +
    geom_text(label = "LDA", x = 2.5, y = 2.4)
  
  # make the second plot
  q <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))
  q <- q + stat_function(fun = function(x) {dnorm(x,mean = mu1.marg,
                         sd = sqrt(var.marg))}, col = "blue") + xlim(-4,4) +
    ylim(0, .85) +
    stat_function(fun = function(x) {dnorm(x,mean = mu2.marg,
                                           sd = sqrt(var.marg))},
                  col = "darkseagreen") +
    theme_minimal()
  return(list(p1 = p, p2 = q))
}
    
angles <- seq(from = pi/2, to = -pi/2, length.out = 20)
angles <- c(angles, rev(angles))
# Finish the code inside saveGIF
saveGIF({
  
  # Loop through all time points
  for (ang in angles) {
  
    
    # Finish the ggplot command
    plots <- anglePlot(ang, mu1, mu2, Sigma, lin.slope)
    
    grid.arrange(plots$p1, plots$p2, ncol = 2, 
                 widths = c(2,1))
    
    
  }
  
}, movie.name = "AzeemZaman_artifact.gif", interval = 1/2)
 
  
  
