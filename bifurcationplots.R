
alpha <- 1.636
alpha1 <- 1.636 - 0.4908
alpha2 <- 1.636
eta <- 20.19
kappa <- 500
delta <- 0.3743
rho <- 1.131
sigma <- 0.1181
mu <- 0.00311

find_roots <- function(f, lower, upper, n = 1000) {
  x <- seq(lower, upper, length.out = n)
  y <- vapply(x, f, numeric(1))
  roots <- numeric(0)
  
  for (i in seq_len(n - 1)) {
    if (is.finite(y[i]) && abs(y[i]) < 1e-9) {
      roots <- c(roots, x[i])
    }
    
    if (is.finite(y[i]) && is.finite(y[i + 1]) &&
        y[i] * y[i + 1] < 0) {
      roots <- c(
        roots,
        uniroot(f, c(x[i], x[i + 1]))$root
      )
    }
  }
  
  unique(round(roots, 8))
}

get_stability <- function(J) {
  values <- eigen(J, only.values = TRUE)$values
  
  if (all(Re(values) < -1e-8)) {
    "Stable"
  } else {
    "Unstable"
  }
}

kuznetsov_eq <- function(source) {
  E_fun <- function(T) {
    alpha * (1 - T / kappa)
  }
  
  f <- function(T) {
    E <- E_fun(T)
    
    E * (
      rho * T / (T + eta) -
        mu * T - delta
    ) + source
  }
  
  roots <- find_roots(f, 1e-6, kappa - 1e-6)
  
  out <- data.frame()
  
  for (T in c(0, roots)) {
    E <- if (T == 0) source / delta else E_fun(T)
    
    J <- matrix(c(
      alpha * (1 - 2 * T / kappa) - E,
      E * (rho * eta / (T + eta)^2 - mu),
      -T,
      rho * T / (T + eta) - mu * T - delta
    ), nrow = 2)
    
    out <- rbind(out, data.frame(
      parameter = source,
      tumour = T,
      effector = E,
      stability = get_stability(J)
    ))
  }
  
  out
}

extended_eq <- function(theta, mechanism) {
  E_fun <- function(T) {
    alpha1 / theta * (1 - T / kappa)
  }
  
  f <- function(T) {
    E <- E_fun(T)
    
    if (mechanism == "RA") {
      activation <- rho * theta * T / (T + eta)
      exhaustion <- mu * theta * T
    } else {
      activation <- rho * T / (T + eta)
      exhaustion <- mu * (2 - theta) * T
    }
    
    E * (activation - exhaustion - delta) + sigma
  }
  
  roots <- find_roots(f, 1e-6, kappa - 1e-6)
  
  out <- data.frame()
  
  for (T in roots) {
    E <- E_fun(T)
    
    J11 <- -alpha1 * T / kappa
    J12 <- -alpha1 * T / kappa
    J13 <- -theta * T
    
    J21 <- 0
    J22 <- alpha2 * (1 - T / kappa) - E
    J23 <- 0
    
    if (mechanism == "RA") {
      J31 <- E * (
        rho * theta * eta / (T + eta)^2 -
          mu * theta
      )
      
      J32 <- E * (
        rho * (1 - theta) * T * eta / (T + eta)^2 -
          mu
      )
    } else {
      J31 <- E * (
        rho * eta / (T + eta)^2 -
          mu * (2 - theta)
      )
      
      J32 <- E * (
        rho * eta / (T + eta)^2 - mu
      )
    }
    
    J33 <- -sigma / E
    
    J <- matrix(c(
      J11, J21, J31,
      J12, J22, J32,
      J13, J23, J33
    ), nrow = 3)
    
    out <- rbind(out, data.frame(
      parameter = theta,
      tumour = T,
      effector = E,
      stability = get_stability(J)
    ))
  }
  
  out
}

sigma_values <- seq(0.01, 0.8, length.out = 300)
theta_values <- seq(0.001, 1, length.out = 500)

simple_data <- do.call(
  rbind,
  lapply(sigma_values, kuznetsov_eq)
)

ra_data <- do.call(
  rbind,
  lapply(theta_values, function(x) extended_eq(x, "RA"))
)

ici_data <- do.call(
  rbind,
  lapply(theta_values, function(x) extended_eq(x, "ICI"))
)

draw_bifurcation <- function(data, title, xlab) {
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  par(
    mfrow = c(1, 2),
    mar = c(5, 5, 4, 1.5),
    oma = c(2, 1, 3, 1),
    mgp = c(2.7, 0.8, 0)
  )
  
  colors <- c(
    Stable = "#187B65",
    Unstable = "#D46A96"
  )
  
  for (variable in c("tumour", "effector")) {
    ylabel <- if (variable == "tumour") {
      "Tumour Cell Population"
    } else {
      "Effector Cell Population"
    }
    
    values <- data[[variable]]
    
    plot(
      NA,
      xlim = range(data$parameter),
      ylim = range(c(0, values), finite = TRUE),
      xlab = xlab,
      ylab = ylabel,
      main = ylabel,
      cex.main = 1.05,
      cex.lab = 0.95,
      cex.axis = 0.85,
      las = 1
    )
    
    grid(col = "gray92")
    
    for (status in c("Unstable", "Stable")) {
      selected <- data$stability == status
      
      points(
        data$parameter[selected],
        values[selected],
        pch = 16,
        cex = 0.42,
        col = colors[status]
      )
    }
    
    legend(
      "topright",
      legend = c("Stable", "Unstable"),
      col = colors,
      pch = 16,
      bty = "n",
      cex = 0.8
    )
  }
  
  mtext(title, outer = TRUE, side = 3,
        line = 0.6, cex = 1.25, font = 2)
}

draw_bifurcation(
  simple_data,
  "Original Kuznetsov Model",
  expression("Effector source rate, " * sigma)
)

draw_bifurcation(
  ra_data,
  "Extended Kuznetsov - Reduced Antigen Presentation",
  expression("Kill-rate reduction, " * theta[1])
)

draw_bifurcation(
  ici_data,
  "Extended Kuznetsov - Immune Checkpoint Evasion",
  expression("Kill-rate reduction, " * theta[1])
)

