
library(deSolve)

simple <- function(t, state, p) {
  T <- state["T"]
  E <- state["E"]
  with(as.list(p), {
    dT <- a*T*(1-b*T) - n*E*T
    dE <- s + p_rate*E*T/(g+T) - m*E*T - d*E
    list(c(T = dT, E = dE))
  })
}

extended <- function(t, state, p) {
  T <- state["T"]
  E <- state["E"]
  with(as.list(p), {
    dT <- (alpha_T*(1-T/kappa)-theta*E)*T
    dE <- (theta*rho*T/(T+eta) -
             theta*mu*T-delta)*E+sigma
    list(c(T = dT, E = dE))
  })
}

p1 <- c(
  s = 1.3e4,
  p_rate = 0.1245,
  g = 2.019e7,
  m = 3.422e-10,
  n = 1.101e-7,
  d = 0.0412,
  a = 0.18,
  b = 2e-9
)

p2 <- c(
  alpha_T = 1.636-0.4908,
  theta = 0.7,
  kappa = 500,
  rho = 1.131,
  eta = 20.19,
  mu = 0.00311,
  delta = 0.3743,
  sigma = 0.1181
)

find_roots <- function(f, lower, upper) {
  xx <- seq(lower, upper, length.out = 3001)
  yy <- sapply(xx, f)
  roots <- numeric(0)
  
  for (i in 1:(length(xx)-1)) {
    if (is.finite(yy[i]) && is.finite(yy[i+1]) &&
        yy[i]*yy[i+1] < 0) {
      roots <- c(roots,
                 uniroot(f, c(xx[i], xx[i+1]))$root)
    }
  }
  unique(roots)
}

E1 <- function(T) {
  p1["a"]*(1-p1["b"]*T)/p1["n"]
}

F1 <- function(T) {
  E <- E1(T)
  p1["s"]+p1["p_rate"]*E*T/(p1["g"]+T)-
    p1["m"]*E*T-p1["d"]*E
}

r1 <- find_roots(F1, 1, 1/p1["b"])

eq1 <- rbind(
  c(0, p1["s"]/p1["d"]),
  cbind(r1, sapply(r1, E1))
)

E2 <- function(T) {
  p2["alpha_T"]*(1-T/p2["kappa"])/p2["theta"]
}

F2 <- function(T) {
  E <- E2(T)
  (p2["theta"]*p2["rho"]*T/(T+p2["eta"])-
      p2["theta"]*p2["mu"]*T-
      p2["delta"])*E+p2["sigma"]
}

r2 <- find_roots(F2, 0.001, p2["kappa"])

eq2 <- rbind(
  c(0, p2["sigma"]/p2["delta"]),
  cbind(r2, sapply(r2, E2))
)

phase_plot <- function(model, pars, initials,
                       eq, title, xlim, ylim,
                       scaleT, scaleE, tmax) {
  
  plot(NA,
       xlim = xlim,
       ylim = ylim,
       xlab = "Tumour Population",
       ylab = "Effector Population",
       main = title,
       cex.main = 1.2,
       cex.lab = 1.1,
       cex.axis = 1,
       las = 1)
  
  grid(col = "gray92")
  
  for (u0 in initials) {
    sol <- ode(
      y = c(T = u0[1], E = u0[2]),
      times = seq(0, tmax, length.out = 2000),
      func = model,
      parms = pars,
      method = "lsoda"
    )
    
    lines(sol[, "T"]*scaleT,
          sol[, "E"]*scaleE,
          col = "blue",
          lwd = 2)
  }
  
  points(eq[,1]*scaleT,
         eq[,2]*scaleE,
         col = "red",
         pch = 19,
         cex = 1.4)
  
  legend("topright",
         legend = c("Trajectories", "Equilibrium"),
         col = c("blue", "red"),
         lty = c(1, NA),
         pch = c(NA, 19),
         bty = "n")
}

draw_simple <- function() {
  par(mar = c(5, 5, 4, 2))
  
  phase_plot(
    simple, p1,
    list(
      c(2e7, 2e5),
      c(8e7, 4e5),
      c(1.5e8, 8e5),
      c(3e8, 5e5)
    ),
    eq1,
    "Simple Kuznetsov Phase Portrait",
    c(0, 5), c(0, 20),
    1e-8, 1e-5, 200
  )
}

draw_extended <- function() {
  par(mar = c(5, 5, 4, 2))
  
  phase_plot(
    extended, p2,
    list(
      c(20, 0.4),
      c(80, 0.8),
      c(180, 1.3),
      c(350, 0.5)
    ),
    eq2,
    "Extended Kuznetsov Phase Portrait",
    c(0, 500), c(0, 2),
    1, 1, 100
  )
}

png("simple_kuznetsov_phase.png",
    width = 1000, height = 800, res = 130)
draw_simple()
dev.off()

png("extended_kuznetsov_phase.png",
    width = 1000, height = 800, res = 130)
draw_extended()
dev.off()

draw_simple()
draw_extended()

print(eq1)
print(eq2)
