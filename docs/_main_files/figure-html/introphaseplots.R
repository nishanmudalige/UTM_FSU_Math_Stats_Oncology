
phase_plot <- function(A, title, initials, tmax) {
  
  par(mar = c(5, 5, 4, 2))
  
  plot(
    NA,
    xlim = c(-4, 4),
    ylim = c(-4, 4),
    xlab = "State Variable x",
    ylab = "State Variable y",
    main = title,
    asp = 1,
    cex.main = 1.3,
    cex.lab = 1.2,
    cex.axis = 1
  )
  
  grid(col = "gray90")
  
  dt <- 0.01
  steps <- round(tmax / dt)
  
  for (initial in initials) {
    
    trajectory <- matrix(
      NA_real_,
      nrow = steps + 1,
      ncol = 2
    )
    
    trajectory[1, ] <- initial
    
    for (i in 1:steps) {
      
      state <- trajectory[i, ]
      
      k1 <- as.vector(A %*% state)
      k2 <- as.vector(A %*% (state + dt * k1 / 2))
      k3 <- as.vector(A %*% (state + dt * k2 / 2))
      k4 <- as.vector(A %*% (state + dt * k3))
      
      trajectory[i + 1, ] <- state +
        dt * (k1 + 2*k2 + 2*k3 + k4) / 6
    }
    
    lines(
      trajectory[, 1],
      trajectory[, 2],
      col = "#2563EB",
      lwd = 2
    )
    
    visible <- which(
      trajectory[, 1] > -3.8 &
        trajectory[, 1] < 3.8 &
        trajectory[, 2] > -3.8 &
        trajectory[, 2] < 3.8
    )
    
    if (length(visible) > 20) {
      
      k <- visible[round(length(visible) * 0.45)]
      j <- min(k + 8, nrow(trajectory))
      
      arrows(
        trajectory[k, 1],
        trajectory[k, 2],
        trajectory[j, 1],
        trajectory[j, 2],
        length = 0.09,
        col = "#2563EB",
        lwd = 1.5
      )
    }
  }
  
  points(
    0, 0,
    pch = 19,
    col = "#DC2626",
    cex = 1.6
  )
  
  legend(
    "topright",
    legend = c("Trajectories", "Equilibrium"),
    col = c("#2563EB", "#DC2626"),
    lty = c(1, NA),
    pch = c(NA, 19),
    lwd = c(2, NA),
    bty = "n",
    cex = 0.9
  )
}

initials <- list(
  c(-3, -3),
  c(-3, 2),
  c(3, 3),
  c(3, -2),
  c(-1, 3),
  c(1, -3),
  c(-3, 0.5),
  c(3, 0.5)
)

A_node <- matrix(
  c(-0.7, 0,
    0, -1.5),
  nrow = 2,
  byrow = TRUE
)

A_spiral <- matrix(
  c(-0.25, -1.4,
    1.4, -0.25),
  nrow = 2,
  byrow = TRUE
)

A_saddle <- matrix(
  c(0.7, 0,
    0, -0.7),
  nrow = 2,
  byrow = TRUE
)

phase_plot(
  A_node,
  "Stable Node Phase Portrait",
  initials,
  7
)

phase_plot(
  A_spiral,
  "Stable Spiral Phase Portrait",
  initials,
  12
)

phase_plot(
  A_saddle,
  "Saddle Point Phase Portrait",
  initials,
  4
)
