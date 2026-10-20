
# Lotka-Volterra Predator-Prey Phase Plot

# Model parameters
alpha <- 1.1
beta <- 0.4
delta <- 0.1
gamma <- 0.4

# Lotka-Volterra equations
lotka_volterra <- function(x, y) {
  dx <- alpha*x - beta*x*y
  dy <- delta*x*y - gamma*y
  return(c(dx, dy))
}

# Time settings
dt <- 0.01
steps <- 10000

# Starting populations
initial_conditions <- list(
  c(4, 1.3),
  c(5, 1.2),
  c(6, 1.1),
  c(7, 1.0),
  c(8, 0.9)
)

# Create phase plot
plot(NA,
     xlim = c(0, 12),
     ylim = c(0, 9),
     xlab = "Prey Population",
     ylab = "Predator Population",
     main = "Lotka-Volterra Phase Portrait")

# Solve and plot each trajectory
for (initial in initial_conditions) {
  
  x <- numeric(steps + 1)
  y <- numeric(steps + 1)
  
  x[1] <- initial[1]
  y[1] <- initial[2]
  
  for (i in 1:steps) {
    
    k1 <- lotka_volterra(x[i], y[i])
    
    k2 <- lotka_volterra(
      x[i] + dt*k1[1]/2,
      y[i] + dt*k1[2]/2
    )
    
    k3 <- lotka_volterra(
      x[i] + dt*k2[1]/2,
      y[i] + dt*k2[2]/2
    )
    
    k4 <- lotka_volterra(
      x[i] + dt*k3[1],
      y[i] + dt*k3[2]
    )
    
    x[i+1] <- x[i] + dt/6 *
      (k1[1] + 2*k2[1] + 2*k3[1] + k4[1])
    
    y[i+1] <- y[i] + dt/6 *
      (k1[2] + 2*k2[2] + 2*k3[2] + k4[2])
  }
  
  lines(x, y, col = "blue", lwd = 2)
}

# Equilibrium point
x_eq <- gamma / delta
y_eq <- alpha / beta

points(x_eq, y_eq,
       col = "red",
       pch = 19,
       cex = 1.5)

legend("topright",
       legend = c("Trajectories", "Equilibrium"),
       col = c("blue", "red"),
       lty = c(1, NA),
       pch = c(NA, 19))

# Display equilibrium coordinates
print(c(x_eq, y_eq))

