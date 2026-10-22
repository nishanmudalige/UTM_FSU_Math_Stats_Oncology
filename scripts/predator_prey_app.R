# =============================================================================
# Predator-prey (sheep - wolf - grass) agent-based model
# R / Shiny / plotly port of the Agents.jl v5.6 example:
#   https://juliadynamics.github.io/Agents.jl/v5.6/examples/predator_prey/
#
# Run with:  shiny::runApp("app.R")   (needs packages: shiny, plotly)
# =============================================================================

library(shiny)
library(plotly)

# ---- Model ------------------------------------------------------------------
# Agents are stored as parallel vectors inside an environment (mutable state).
# type: 1 = sheep, 2 = wolf.  Grid is periodic, Chebyshev metric (8 neighbours).

SHEEP <- 1L
WOLF  <- 2L
DX <- c(-1L, 0L, 1L, -1L, 1L, -1L, 0L, 1L)
DY <- c(-1L, -1L, -1L, 0L, 0L, 1L, 1L, 1L)

# Julia seeds do not carry over to R's RNG; seed 24 gives long-run coexistence here.
presets <- list(
  "Default (wolves die out)" = list(
    n_sheep = 100, n_wolves = 50, dims = 20, regrowth_time = 30,
    de_sheep = 4, de_wolf = 20, sheep_rep = 0.04, wolf_rep = 0.05, seed = 23182
  ),
  "Stable coexistence" = list(
    n_sheep = 140, n_wolves = 20, dims = 30, regrowth_time = 30,
    de_sheep = 5, de_wolf = 30, sheep_rep = 0.31, wolf_rep = 0.06, seed = 24
  )
)

initialize_model <- function(n_sheep = 100, n_wolves = 50, dims = 20,
                             regrowth_time = 30, de_sheep = 4, de_wolf = 20,
                             sheep_rep = 0.04, wolf_rep = 0.05, seed = 23182) {
  set.seed(seed)
  m <- new.env()
  W <- H <- as.integer(dims)
  n <- n_sheep + n_wolves
  m$W <- W
  m$H <- H

  # Agents: random positions; initial energy uniform on 0, ..., 2*delta - 1
  m$type   <- c(rep(SHEEP, n_sheep), rep(WOLF, n_wolves))
  m$x      <- sample.int(W, n, replace = TRUE)
  m$y      <- sample.int(H, n, replace = TRUE)
  m$energy <- as.numeric(c(sample.int(2 * de_sheep, n_sheep, replace = TRUE) - 1,
                           sample.int(2 * de_wolf,  n_wolves, replace = TRUE) - 1))
  m$alive  <- rep(TRUE, n)
  m$rp     <- c(sheep_rep, wolf_rep)   # reproduction probability by type
  m$de     <- c(de_sheep, de_wolf)     # energy gained from eating, by type

  # Occupancy counts per cell (needed for the sheep's "move only to empty cell")
  m$occ <- matrix(tabulate((m$y - 1L) * W + m$x, W * H), W, H)

  # Grass: half fully grown (countdown = regrowth_time), rest random countdown
  m$regrowth  <- regrowth_time
  m$fully     <- matrix(sample(c(TRUE, FALSE), W * H, replace = TRUE), W, H)
  m$countdown <- matrix(ifelse(m$fully, regrowth_time,
                               sample.int(regrowth_time, W * H, replace = TRUE) - 1L),
                        W, H)

  m$step <- 0L
  m$h_step <- 0L
  m$h_sheep <- n_sheep
  m$h_wolf <- n_wolves
  m$h_grass <- sum(m$fully)
  m
}

kill_agent <- function(m, i) {
  m$alive[i] <- FALSE
  m$occ[m$x[i], m$y[i]] <- m$occ[m$x[i], m$y[i]] - 1L
}

reproduce <- function(m, i) {
  # Offspring costs half the parent's energy; it's placed on the parent's cell
  m$energy[i] <- m$energy[i] / 2
  m$type   <- c(m$type, m$type[i])
  m$x      <- c(m$x, m$x[i])
  m$y      <- c(m$y, m$y[i])
  m$energy <- c(m$energy, m$energy[i])
  m$alive  <- c(m$alive, TRUE)
  m$occ[m$x[i], m$y[i]] <- m$occ[m$x[i], m$y[i]] + 1L
}

agent_step <- function(m, i) {
  t <- m$type[i]

  # walk!: move to a random neighbouring cell. Sheep only move to empty cells
  # (Agents.jl default ifempty = true); wolves move anywhere (ifempty = false).
  nx <- ((m$x[i] - 1L + DX) %% m$W) + 1L
  ny <- ((m$y[i] - 1L + DY) %% m$H) + 1L
  if (t == SHEEP) {
    ok <- m$occ[cbind(nx, ny)] == 0L
    nx <- nx[ok]
    ny <- ny[ok]
  }
  if (length(nx) > 0) {
    k <- sample.int(length(nx), 1)
    m$occ[m$x[i], m$y[i]] <- m$occ[m$x[i], m$y[i]] - 1L
    m$x[i] <- nx[k]
    m$y[i] <- ny[k]
    m$occ[nx[k], ny[k]] <- m$occ[nx[k], ny[k]] + 1L
  }

  # Moving costs 1 energy; die if energy < 0
  m$energy[i] <- m$energy[i] - 1
  if (m$energy[i] < 0) {
    kill_agent(m, i)
    return(invisible())
  }

  # Eat
  if (t == SHEEP) {
    if (m$fully[m$x[i], m$y[i]]) {
      m$energy[i] <- m$energy[i] + m$de[SHEEP]
      m$fully[m$x[i], m$y[i]] <- FALSE
    }
  } else {
    dinner <- which(m$alive & m$type == SHEEP & m$x == m$x[i] & m$y == m$y[i])
    if (length(dinner) > 0) {
      kill_agent(m, dinner[1])
      m$energy[i] <- m$energy[i] + m$de[WOLF]
    }
  }

  # Reproduce
  if (runif(1) <= m$rp[t]) reproduce(m, i)
}

grass_step <- function(m) {
  growing <- !m$fully
  ready   <- growing & m$countdown <= 0
  m$fully[ready]     <- TRUE
  m$countdown[ready] <- m$regrowth
  dec <- growing & !ready
  m$countdown[dec] <- m$countdown[dec] - 1L
}

model_step <- function(m) {
  # Schedulers.randomly: agents alive at the start of the step, in random order.
  # Newborns wait until the next step; agents eaten mid-step are skipped.
  ids <- which(m$alive)
  ids <- ids[sample.int(length(ids))]
  for (i in ids) if (m$alive[i]) agent_step(m, i)
  grass_step(m)

  m$step <- m$step + 1L
  m$h_step  <- c(m$h_step, m$step)
  m$h_sheep <- c(m$h_sheep, sum(m$alive & m$type == SHEEP))
  m$h_wolf  <- c(m$h_wolf,  sum(m$alive & m$type == WOLF))
  m$h_grass <- c(m$h_grass, sum(m$fully))

  # Occasionally drop dead agents so the vectors don't grow without bound
  if (length(m$alive) > 4 * sum(m$alive) + 500) {
    keep <- m$alive
    m$type <- m$type[keep]
    m$x <- m$x[keep]
    m$y <- m$y[keep]
    m$energy <- m$energy[keep]
    m$alive <- m$alive[keep]
  }
}

# ---- Plots ------------------------------------------------------------------
col_sheep_ts <- "#8B8878"   # cornsilk4, as in the Julia example
col_wolf     <- "#33334D"
col_grass    <- "#2E8B2E"

plot_grid <- function(m) {
  W <- m$W
  H <- m$H
  # Grass colour: 1 = fully grown (green); otherwise fraction regrown (brown -> green)
  growth <- ifelse(m$fully, 1, 1 - m$countdown / m$regrowth)

  s <- which(m$alive & m$type == SHEEP)
  w <- which(m$alive & m$type == WOLF)

  p <- plot_ly() |>
    add_heatmap(
      x = seq_len(W), y = seq_len(H), z = t(growth),
      zmin = 0, zmax = 1, showscale = FALSE,
      colorscale = list(c(0, "#A0522D"), c(1, "#3C9A3C")),
      hovertemplate = "cell (%{x}, %{y})<br>grass regrown: %{z:.0%}<extra></extra>"
    )

  if (length(s) > 0) {
    p <- p |> add_markers(
      x = m$x[s] - 0.15, y = m$y[s] - 0.15 * runif(length(s)),
      name = "Sheep",
      marker = list(symbol = "circle", size = 9, color = "rgba(255,255,255,0.85)",
                    line = list(color = "black", width = 1)),
      text = sprintf("Sheep<br>energy: %.1f", m$energy[s]),
      hoverinfo = "text"
    )
  }
  if (length(w) > 0) {
    p <- p |> add_markers(
      x = m$x[w] + 0.15, y = m$y[w] + 0.15 * runif(length(w)),
      name = "Wolves",
      marker = list(symbol = "triangle-up", size = 11, color = "rgba(51,51,77,0.9)",
                    line = list(color = "black", width = 1)),
      text = sprintf("Wolf<br>energy: %.1f", m$energy[w]),
      hoverinfo = "text"
    )
  }

  p |>
    layout(
      title = list(text = sprintf("Step %d  |  sheep: %d  wolves: %d  grass: %d",
                                  m$step, length(s), length(w), sum(m$fully)),
                   font = list(size = 14)),
      xaxis = list(range = c(0.5, W + 0.5), showgrid = FALSE, zeroline = FALSE,
                   title = "", constrain = "domain"),
      yaxis = list(range = c(0.5, H + 0.5), showgrid = FALSE, zeroline = FALSE,
                   title = "", scaleanchor = "x", constrain = "domain"),
      legend = list(orientation = "h", y = -0.05),
      margin = list(t = 40, l = 10, r = 10, b = 10),
      uirevision = "grid"
    ) |>
    config(displayModeBar = FALSE)
}

plot_timeseries <- function(m) {
  plot_ly(x = m$h_step) |>
    add_lines(y = m$h_sheep, name = "Sheep",  line = list(color = col_sheep_ts, width = 2)) |>
    add_lines(y = m$h_wolf,  name = "Wolves", line = list(color = col_wolf, width = 2)) |>
    add_lines(y = m$h_grass, name = "Grass (fully grown patches)",
              line = list(color = col_grass, width = 2)) |>
    layout(
      xaxis = list(title = "Step"),
      yaxis = list(title = "Population", rangemode = "tozero"),
      hovermode = "x unified",
      legend = list(orientation = "h", y = -0.2),
      margin = list(t = 20)
    )
}

# ---- UI ---------------------------------------------------------------------
ui <- fluidPage(
  titlePanel("Predator-prey dynamics: sheep, wolves and grass"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("preset", "Parameter preset", choices = names(presets)),
      sliderInput("n_sheep", "Initial sheep", 0, 400, 100, step = 5),
      sliderInput("n_wolves", "Initial wolves", 0, 200, 50, step = 5),
      sliderInput("dims", "Grid size (n x n)", 5, 50, 20),
      sliderInput("regrowth_time", "Grass regrowth time", 1, 100, 30),
      sliderInput("de_sheep", "Energy from grass (sheep)", 1, 20, 4),
      sliderInput("de_wolf", "Energy from sheep (wolf)", 1, 60, 20),
      sliderInput("sheep_rep", "Sheep reproduction prob.", 0, 0.5, 0.04, step = 0.01),
      sliderInput("wolf_rep", "Wolf reproduction prob.", 0, 0.5, 0.05, step = 0.01),
      numericInput("seed", "Random seed", 23182),
      actionButton("reset", "Reset model", class = "btn-primary", width = "100%"),
      helpText("Changing a parameter takes effect when you reset the model.")
    ),
    mainPanel(
      width = 9,
      fluidRow(
        column(3, actionButton("play", "Play", width = "100%")),
        column(3, actionButton("step1", "Step once", width = "100%")),
        column(3, numericInput("n_run", NULL, 500, min = 1, step = 100, width = "100%")),
        column(3, actionButton("run_n", "Run N steps", width = "100%"))
      ),
      fluidRow(
        column(6, sliderInput("speed", "Frame delay (ms)", 50, 1000, 200, step = 50,
                              width = "100%")),
        column(6, sliderInput("spf", "Steps per frame", 1, 20, 1, width = "100%"))
      ),
      fluidRow(
        column(6, plotlyOutput("grid", height = "520px")),
        column(6, plotlyOutput("ts", height = "520px"))
      )
    )
  )
)

# ---- Server -----------------------------------------------------------------
server <- function(input, output, session) {
  rv <- reactiveValues(m = NULL, tick = 0)
  running <- reactiveVal(FALSE)

  set_running <- function(on) {
    running(on)
    updateActionButton(session, "play", label = if (on) "Pause" else "Play")
  }

  observeEvent(input$preset, {
    p <- presets[[input$preset]]
    for (nm in c("n_sheep", "n_wolves", "dims", "regrowth_time", "de_sheep",
                 "de_wolf", "sheep_rep", "wolf_rep")) {
      updateSliderInput(session, nm, value = p[[nm]])
    }
    updateNumericInput(session, "seed", value = p$seed)
  }, ignoreInit = TRUE)

  reset_model <- function() {
    set_running(FALSE)
    rv$m <- initialize_model(
      n_sheep = input$n_sheep, n_wolves = input$n_wolves, dims = input$dims,
      regrowth_time = input$regrowth_time, de_sheep = input$de_sheep,
      de_wolf = input$de_wolf, sheep_rep = input$sheep_rep,
      wolf_rep = input$wolf_rep, seed = input$seed
    )
    rv$tick <- rv$tick + 1
  }
  observeEvent(input$reset, reset_model(), ignoreNULL = FALSE)

  observeEvent(input$play, set_running(!running()))

  observeEvent(input$step1, {
    model_step(rv$m)
    rv$tick <- rv$tick + 1
  })

  observeEvent(input$run_n, {
    set_running(FALSE)
    n <- max(1, as.integer(input$n_run))
    withProgress(message = "Running model", value = 0, {
      for (k in seq_len(n)) {
        model_step(rv$m)
        if (k %% 50 == 0) incProgress(50 / n)
      }
    })
    rv$tick <- rv$tick + 1
  })

  observe({
    if (!running()) return()
    invalidateLater(input$speed)
    isolate({
      for (k in seq_len(input$spf)) model_step(rv$m)
      rv$tick <- rv$tick + 1
    })
  })

  output$grid <- renderPlotly({
    rv$tick
    req(rv$m)
    plot_grid(rv$m)
  })

  output$ts <- renderPlotly({
    rv$tick
    req(rv$m)
    plot_timeseries(rv$m)
  })
}

shinyApp(ui, server)
