# =============================================================================
# T Cells vs. Tumour: An Agent-Based Immune-Oncology Model in R / Shiny / Plotly
# =============================================================================
#
# Biological Background:
# ----------------------
# This model translates the classic Lotka-Volterra predator-prey dynamics into
# a spatial, cell-level simulation of the tumor microenvironment:
#   - Tumour Cells (Prey):
#       Consume vascular nutrients, divide into adjacent empty space (respecting
#       contact inhibition / volume exclusion), and die from starvation (necrosis).
#   - Cytotoxic T Lymphocytes (Predators):
#       Patrol tissue, navigate toward cancer cells via chemotaxis, kill tumour cells
#       upon physical contact, and undergo antigen-driven clonal proliferation.
#       Chronic killing induces T-cell exhaustion, rendering them inactive.
#   - Vascular Nutrient Field (Resource):
#       Models glucose/oxygen delivered by microvessels. Consumed patches deplete
#       and require a countdown period to replenish.
#   - Immune Checkpoint Blockade (Therapy):
#       Simulates anti-PD-1 immunotherapy, which reinvigorates exhausted T cells
#       and restores cytotoxic killing efficiency.
#
# How to Run:
# -----------
# In R or RStudio:
#   install.packages(c("shiny", "plotly"))
#   shiny::runApp("scripts/tumour_immune_app.R")
# =============================================================================

library(shiny)
library(plotly)

# =============================================================================
# SECTION 1: MODEL CONSTANTS & DEFAULT PARAMETERS
# =============================================================================

# Cell type identifiers (using integers is fast and memory-efficient in R)
TUMOUR <- 1L  # 1 represents a tumour cell
TCELL  <- 2L  # 2 represents a cytotoxic T cell

# Moore 8-neighbourhood coordinate offsets:
# These vectors represent the 8 possible steps on a 2D grid:
# (North, South, East, West, and all 4 diagonals).
DX <- c(-1L,  0L,  1L, -1L, 1L, -1L, 0L, 1L)
DY <- c(-1L, -1L, -1L,  0L, 0L,  1L, 1L, 1L)

# Default parameter list:
# Storing parameters in a named list makes it easy to modify, pass to functions,
# and define preset scenarios.
default_params <- list(
  dims               = 40,    # Lattice dimensions (40x40 grid = 1,600 sites)
  n_tumour           = 60,    # Initial number of tumour cells
  n_tcell            = 10,    # Initial number of patrolling T cells
  resupply_time      = 12,    # Steps required for depleted nutrients to replenish
  de_tumour          = 14,    # Energy gained by a tumour cell from eating nutrient
  tumour_div         = 0.06,  # Probability of tumour cell division per step
  tumour_motility    = 0.05,  # Probability of a tumour cell migrating to an adjacent site
  de_tcell           = 12,    # Energy gained by a T cell upon killing a tumour cell
  p_kill             = 0.50,  # Base probability of a T cell killing a contacted tumour cell
  tcell_prolif       = 0.40,  # Probability of T-cell division (clonal expansion) after a kill
  exhaustion_limit   = 6,     # Max kills before a T cell becomes functionally exhausted
  recruit            = 1.0,   # Baseline T-cell recruitment rate (per 100 tumour cells)
  chemotaxis         = 0.60,  # Probability of moving toward a neighboring tumour cell
  therapy            = FALSE, # Is checkpoint inhibitor therapy active?
  therapy_start      = 200,   # Simulation step when therapy is administered
  therapy_kill_boost = 1.6,   # Multiplier for kill probability during therapy
  therapy_exh_mult   = 4,     # Multiplier for exhaustion threshold during therapy
  seed               = 1      # Random seed for reproducible runs
)

# =============================================================================
# SECTION 2: INITIALIZING THE SIMULATION ENVIRONMENT
# =============================================================================
# Why use an R "environment" (`new.env()`)?
# In R, passing large data frames or lists to functions usually makes copies
# in memory ("copy-on-modify"). For an agent-based simulation that updates
# hundreds of times per second, modifying variables inside a mutable environment
# allows us to update agent states in place, making the code run much faster!

initialize_model <- function(p = default_params) {
  set.seed(p$seed)  # Set random seed so results are reproducible
  m <- new.env()    # Create a fresh, mutable environment to hold the world state
  
  W <- H <- as.integer(p$dims)
  m$W <- W
  m$H <- H
  m$p <- p

  # --- Seed the Tumour Nodule ---
  # Tumours do not appear scattered at random; they originate as a compact nodule.
  # We compute the geometric center and rank all grid points by distance to center
  # (plus a tiny random jitter) so the initial tumour forms a realistic round cluster.
  cx <- (W + 1) / 2
  cy <- (H + 1) / 2
  cand <- expand.grid(x = seq_len(W), y = seq_len(H))
  d <- sqrt((cand$x - cx)^2 + (cand$y - cy)^2)
  cand <- cand[order(d + runif(nrow(cand), 0, 0.5)), ]
  
  nt <- min(p$n_tumour, W * H)
  tx <- cand$x[seq_len(nt)]
  ty <- cand$y[seq_len(nt)]

  # --- Seed Patrolling T Cells ---
  # T cells begin scattered randomly across the tissue domain.
  ncel <- p$n_tcell
  n_total <- nt + ncel

  # --- Parallel Vectors (Agent Attributes) ---
  # Instead of creating separate objects for each cell, we keep synchronized vectors:
  # agent `i` has type `m$type[i]`, coordinates `(m$x[i], m$y[i])`, and `m$energy[i]`.
  m$type   <- c(rep(TUMOUR, nt), rep(TCELL, ncel))
  m$x      <- c(tx, sample.int(W, ncel, replace = TRUE))
  m$y      <- c(ty, sample.int(H, ncel, replace = TRUE))
  m$energy <- as.numeric(c(
    sample.int(2 * p$de_tumour, nt, replace = TRUE) - 1,
    sample.int(2 * p$de_tcell,  ncel, replace = TRUE) - 1
  ))
  m$kills  <- integer(n_total)   # Kills performed by each agent (tracked for T cells)
  m$alive  <- rep(TRUE, n_total) # Status flag: TRUE = active, FALSE = dead/removed

  # --- Spatial Grid Matrix (Volume Exclusion) ---
  # `m$tum` is a 2D matrix where `m$tum[x, y]` stores the agent ID of the tumour cell
  # occupying site (x, y), or 0 if the site has no tumour cell.
  # This guarantees that at most ONE tumour cell occupies any grid patch.
  m$tum <- matrix(0L, W, H)
  m$tum[cbind(tx, ty)] <- seq_len(nt)

  # --- Vascular Nutrient Field ---
  # `m$fully` is TRUE if nutrient is currently present at site (x, y).
  # `m$countdown` tracks how many steps remain before depleted nutrient replenishes.
  m$fully     <- matrix(sample(c(TRUE, FALSE), W * H, replace = TRUE), W, H)
  m$countdown <- matrix(
    ifelse(m$fully, p$resupply_time, sample.int(p$resupply_time, W * H, replace = TRUE) - 1L),
    W, H
  )

  # History tracking vectors for plotting population time series
  m$step        <- 0L
  m$kill_events <- 0L
  m$h_step      <- integer(0)
  m$h_tumour    <- integer(0)
  m$h_active    <- integer(0)
  m$h_exhausted <- integer(0)
  m$h_nutrient  <- integer(0)
  
  record(m)  # Record initial census at step 0
  m
}

# =============================================================================
# SECTION 3: THERAPY & STATE HELPER FUNCTIONS
# =============================================================================

# Check if checkpoint inhibitor therapy is active at the current step
therapy_on <- function(m) {
  m$p$therapy && m$step >= m$p$therapy_start
}

# Dynamic exhaustion threshold (increased when therapy is active)
exh_limit <- function(m) {
  m$p$exhaustion_limit * (if (therapy_on(m)) m$p$therapy_exh_mult else 1)
}

# Dynamic kill probability (boosted when therapy is active)
kill_prob <- function(m) {
  min(1.0, m$p$p_kill * (if (therapy_on(m)) m$p$therapy_kill_boost else 1))
}

# Add a newly born or recruited agent to the simulation
add_agent <- function(m, type, x, y, energy, kills = 0L) {
  m$type   <- c(m$type, type)
  m$x      <- c(m$x, x)
  m$y      <- c(m$y, y)
  m$energy <- c(m$energy, energy)
  m$kills  <- c(m$kills, kills)
  m$alive  <- c(m$alive, TRUE)
  new_id   <- length(m$alive)
  if (type == TUMOUR) m$tum[x, y] <- new_id
  new_id
}

# Mark an agent as dead and clear its position from the tumour spatial matrix
remove_agent <- function(m, i) {
  m$alive[i] <- FALSE
  if (m$type[i] == TUMOUR) m$tum[m$x[i], m$y[i]] <- 0L
}

# Find all 8 neighboring coordinates with periodic (toroidal) boundary conditions:
# The `%%` modulo operator wraps cells exiting one edge around to the opposite edge.
neighbours <- function(m, i) {
  list(
    x = ((m$x[i] - 1L + DX) %% m$W) + 1L,
    y = ((m$y[i] - 1L + DY) %% m$H) + 1L
  )
}

# =============================================================================
# SECTION 4: INDIVIDUAL AGENT BEHAVIOR RULES
# =============================================================================

# --- Tumour Cell Rules ---
# 1. Motility: occasionally crawls into an unoccupied neighboring patch.
# 2. Basal Metabolism: expends 1 energy unit. If energy < 0, dies (necrosis).
# 3. Nutrient Uptake: consumes nutrient at its site if available.
# 4. Proliferation: divides into a free neighbor site (contact inhibition: no space = no division).
tumour_step <- function(m, i) {
  p <- m$p
  nb <- neighbours(m, i)
  free <- m$tum[cbind(nb$x, nb$y)] == 0L  # Find neighboring sites with no tumour cell

  # Step 1: Migration (with probability tumour_motility)
  if (any(free) && runif(1) < p$tumour_motility) {
    k <- which(free)
    k <- k[sample.int(length(k), 1)]       # Pick one random free neighboring site
    m$tum[m$x[i], m$y[i]] <- 0L            # Vacate old site
    m$x[i] <- nb$x[k]; m$y[i] <- nb$y[k]   # Move to new site
    m$tum[m$x[i], m$y[i]] <- i             # Claim new site in spatial matrix
    
    # Refresh neighborhood after moving
    nb <- neighbours(m, i)
    free <- m$tum[cbind(nb$x, nb$y)] == 0L
  }

  # Step 2: Metabolism (costs 1 unit of energy each time step)
  m$energy[i] <- m$energy[i] - 1
  if (m$energy[i] < 0) {
    remove_agent(m, i)  # Starvation: cell undergoes necrotic death
    return(invisible())
  }

  # Step 3: Feeding (absorb local glucose/oxygen from vasculature)
  if (m$fully[m$x[i], m$y[i]]) {
    m$energy[i] <- m$energy[i] + p$de_tumour
    m$fully[m$x[i], m$y[i]] <- FALSE  # Nutrient depleted at this site
  }

  # Step 4: Division (requires both available space and passing division probability)
  if (any(free) && runif(1) <= p$tumour_div) {
    k <- which(free)
    k <- k[sample.int(length(k), 1)]
    m$energy[i] <- m$energy[i] / 2  # Energy is split evenly with daughter cell
    add_agent(m, TUMOUR, nb$x[k], nb$y[k], m$energy[i])
  }
}

# --- Cytotoxic T-Cell Rules ---
# 1. Chemotaxis: moves toward neighboring tumour cells with probability `chemotaxis`;
#    otherwise performs a random walk step.
# 2. Basal Metabolism: expends 1 energy unit. If energy < 0, dies (finite lifespan).
# 3. Exhaustion Check: if total kills reach the exhaustion limit, cell cannot attack or divide.
# 4. Predation: attacks tumour cell at current site with probability `kill_prob`.
# 5. Clonal Expansion: after a kill, divides with probability `tcell_prolif`.
tcell_step <- function(m, i) {
  p <- m$p
  nb <- neighbours(m, i)

  # Step 1: Movement with Chemotaxis (sensing tumour-derived chemokine signals)
  tgt <- which(m$tum[cbind(nb$x, nb$y)] > 0L)
  k <- if (length(tgt) > 0 && runif(1) < p$chemotaxis) {
    tgt[sample.int(length(tgt), 1)]  # Step directly toward a nearby tumour cell
  } else {
    sample.int(8L, 1)                # Random walk step
  }
  m$x[i] <- nb$x[k]
  m$y[i] <- nb$y[k]

  # Step 2: Metabolism & Natural Lifespan
  m$energy[i] <- m$energy[i] - 1
  if (m$energy[i] < 0) {
    remove_agent(m, i)  # Finite lifespan: T cell dies of senescence/exhaustion
    return(invisible())
  }

  # Step 3: Exhaustion Check
  # Exhausted T cells remain physically present but cannot kill or divide.
  if (m$kills[i] >= exh_limit(m)) {
    return(invisible())
  }

  # Step 4: Contact Cytotoxicity (Attacking Tumour Cell at Current Site)
  target <- m$tum[m$x[i], m$y[i]]
  if (target > 0L && runif(1) <= kill_prob(m)) {
    remove_agent(m, target)               # Destroy tumour cell (cell lysis)
    m$kill_events <- m$kill_events + 1L
    m$energy[i]   <- m$energy[i] + p$de_tcell
    m$kills[i]    <- m$kills[i] + 1L

    # Step 5: Antigen-Driven Clonal Expansion (Proliferation)
    if (runif(1) <= p$tcell_prolif) {
      m$energy[i] <- m$energy[i] / 2     # Parent splits energy with daughter T cell
      add_agent(m, TCELL, m$x[i], m$y[i], m$energy[i], m$kills[i])
    }
  }
}

# =============================================================================
# SECTION 5: ENVIRONMENTAL DYNAMICS & SYSTEMIC RECRUITMENT
# =============================================================================

# Nutrient Regrowth: Depleted sites count down and regenerate after `resupply_time`
nutrient_step <- function(m) {
  R <- m$p$resupply_time
  growing <- !m$fully
  ready   <- growing & m$countdown <= 0
  m$fully[ready]     <- TRUE
  m$countdown[ready] <- R
  dec <- growing & !ready
  m$countdown[dec]   <- m$countdown[dec] - 1L
}

# Systemic T-Cell Recruitment:
# Fresh T cells extravasate from the bloodstream into tissue at a rate proportional
# to tumour burden (antigen load). The count is sampled from a Poisson distribution.
recruit_tcells <- function(m) {
  n_tum <- sum(m$alive & m$type == TUMOUR)
  n_new <- rpois(1, m$p$recruit * n_tum / 100)
  for (j in seq_len(n_new)) {
    add_agent(
      m, TCELL,
      sample.int(m$W, 1), sample.int(m$H, 1),
      sample.int(2 * m$p$de_tcell, 1) - 1
    )
  }
}

# Count populations across all living agents
counts <- function(m) {
  a   <- m$alive
  tc  <- a & m$type == TCELL
  exh <- tc & (m$kills >= exh_limit(m))
  c(
    tumour    = sum(a & m$type == TUMOUR),
    active    = sum(tc & !exh),
    exhausted = sum(exh),
    nutrient  = sum(m$fully)
  )
}

# Record macroscopic population census for history plots
record <- function(m) {
  cc <- counts(m)
  m$h_step      <- c(m$h_step, m$step)
  m$h_tumour    <- c(m$h_tumour, cc[["tumour"]])
  m$h_active    <- c(m$h_active, cc[["active"]])
  m$h_exhausted <- c(m$h_exhausted, cc[["exhausted"]])
  m$h_nutrient  <- c(m$h_nutrient, cc[["nutrient"]])
}

# Master simulation step: advances the clock by 1 unit of time (Delta t)
model_step <- function(m) {
  # Asynchronous scheduling: shuffle agent indices so activation order is random
  ids <- which(m$alive)
  ids <- ids[sample.int(length(ids))]
  
  for (i in ids) {
    if (!m$alive[i]) next
    if (m$type[i] == TUMOUR) tumour_step(m, i) else tcell_step(m, i)
  }
  
  nutrient_step(m)   # Regenerate depleted nutrients
  recruit_tcells(m)  # Influx of fresh T cells from circulation
  m$step <- m$step + 1L
  record(m)          # Update population census

  # Periodic garbage collection: prune inactive agents to keep vectors compact
  if (length(m$alive) > 4 * sum(m$alive) + 500) {
    keep <- m$alive
    old_ids <- which(keep)
    for (v in c("type", "x", "y", "energy", "kills", "alive")) m[[v]] <- m[[v]][keep]
    remap <- integer(length(keep))
    remap[old_ids] <- seq_along(old_ids)
    m$tum[m$tum > 0L] <- remap[m$tum[m$tum > 0L]]
  }
}

# =============================================================================
# SECTION 6: PRESETS (THE THREE "E"S OF CANCER IMMUNOEDITING + THERAPY)
# =============================================================================

presets <- list(
  "Equilibrium (immune control)" = list(),
  "Elimination (strong immune response)" = list(recruit = 8),
  "Escape (T-cell exhaustion)" = list(recruit = 1, exhaustion_limit = 2),
  "Escape + checkpoint inhibitor at step 200" = list(
    recruit = 1, exhaustion_limit = 2, therapy = TRUE, therapy_start = 200,
    therapy_kill_boost = 1.6, therapy_exh_mult = 4
  ),
  "Escape + weaker therapy (partial response)" = list(
    recruit = 1, exhaustion_limit = 2, therapy = TRUE, therapy_start = 200,
    therapy_kill_boost = 1.2, therapy_exh_mult = 2
  )
)
param_names <- setdiff(names(default_params), c("therapy", "seed"))

# =============================================================================
# SECTION 7: INTERACTIVE PLOTLY VISUALIZATIONS
# =============================================================================
# Plotly creates interactive HTML widgets that run directly inside any modern
# web browser without requiring an active R server.

col_tumour    <- "#8E2C83"  # Purple for cancer cells
col_active    <- "#1F77B4"  # Vibrant blue for active cytotoxic T cells
col_exhausted <- "#8C8C8C"  # Neutral grey for exhausted T cells
col_nutrient  <- "#C9A227"  # Warm golden amber for vascular nutrients

plot_grid <- function(m) {
  W <- m$W; H <- m$H; R <- m$p$resupply_time
  level <- ifelse(m$fully, 1, 1 - m$countdown / R)
  lim <- exh_limit(m)
  tu  <- which(m$alive & m$type == TUMOUR)
  tc  <- which(m$alive & m$type == TCELL)
  ex  <- tc[m$kills[tc] >= lim]
  ac  <- tc[m$kills[tc] < lim]

  p <- plot_ly() |>
    add_heatmap(
      x = seq_len(W), y = seq_len(H), z = t(level),
      zmin = 0, zmax = 1, showscale = FALSE,
      colorscale = list(c(0, "#5A3A2E"), c(1, "#F3DCD4")),
      hovertemplate = "site (%{x}, %{y})<br>nutrient resupplied: %{z:.0%}<extra></extra>"
    )

  if (length(tu) > 0) {
    p <- p |> add_markers(
      x = m$x[tu], y = m$y[tu], name = "Tumour cells",
      marker = list(symbol = "circle", size = 7, color = col_tumour,
                    line = list(color = "white", width = 0.5)),
      text = sprintf("Tumour cell<br>energy: %.1f", m$energy[tu]),
      hoverinfo = "text"
    )
  }
  if (length(ac) > 0) {
    p <- p |> add_markers(
      x = m$x[ac] + 0.2, y = m$y[ac] + 0.2 * runif(length(ac)), name = "Active T cells",
      marker = list(symbol = "triangle-up", size = 9, color = col_active,
                    line = list(color = "white", width = 0.8)),
      text = sprintf("Active T cell<br>energy: %.1f<br>kills: %d", m$energy[ac], m$kills[ac]),
      hoverinfo = "text"
    )
  }
  if (length(ex) > 0) {
    p <- p |> add_markers(
      x = m$x[ex] + 0.2, y = m$y[ex] + 0.2 * runif(length(ex)), name = "Exhausted T cells",
      marker = list(symbol = "triangle-up", size = 9, color = col_exhausted,
                    line = list(color = "white", width = 0.8)),
      text = sprintf("Exhausted T cell<br>energy: %.1f<br>kills: %d", m$energy[ex], m$kills[ex]),
      hoverinfo = "text"
    )
  }

  status <- if (!m$p$therapy) "" else if (therapy_on(m)) "  |  therapy ON" else
    sprintf("  |  therapy at step %d", m$p$therapy_start)

  p |>
    layout(
      title = list(
        text = sprintf("Step %d  |  tumour: %d  T cells: %d (%d exhausted)%s",
                       m$step, length(tu), length(tc), length(ex), status),
        font = list(size = 13)
      ),
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
  shapes <- list()
  if (m$p$therapy && m$step >= m$p$therapy_start) {
    shapes <- list(list(
      type = "line", x0 = m$p$therapy_start, x1 = m$p$therapy_start,
      y0 = 0, y1 = 1, yref = "paper",
      line = list(color = "#2CA02C", dash = "dot", width = 2)
    ))
  }
  
  plot_ly(x = m$h_step) |>
    add_lines(y = m$h_tumour, name = "Tumour cells",
              line = list(color = col_tumour, width = 2.5)) |>
    add_lines(y = m$h_active, name = "Active T cells",
              line = list(color = col_active, width = 2)) |>
    add_lines(y = m$h_exhausted, name = "Exhausted T cells",
              line = list(color = col_exhausted, width = 2, dash = "dash")) |>
    add_lines(y = m$h_nutrient, name = "Nutrient-available sites",
              line = list(color = col_nutrient, width = 1.5), visible = "legendonly") |>
    layout(
      xaxis = list(title = "Simulation Step"),
      yaxis = list(title = "Population Count", rangemode = "tozero"),
      hovermode = "x unified",
      shapes = shapes,
      annotations = if (length(shapes)) list(list(
        x = m$p$therapy_start, y = 1, yref = "paper", text = "checkpoint inhibitor",
        showarrow = FALSE, xanchor = "left", font = list(color = "#2CA02C", size = 11)
      )) else list(),
      legend = list(orientation = "h", y = -0.2),
      margin = list(t = 20)
    )
}

# =============================================================================
# SECTION 8: SHINY USER INTERFACE (UI) & SERVER LOGIC
# =============================================================================

rules_md <- "
### Tumour-Immune Agent-Based Model Rules
This model translates the classic predator-prey dynamics into an oncology microenvironment:
- **Tumour cells** represent prey: they feed on background nutrients and divide.
- **Cytotoxic T cells** represent predators: they hunt and kill tumour cells.
- **Nutrients** represent oxygen and glucose supplied by the vascular bed.

#### Tumour Cells (at most one per grid site)
1. **Motility:** With probability *tumour motility*, migrate to an empty neighboring site.
2. **Metabolism:** Consume 1 unit of energy each step; if energy drops below 0, the cell dies (necrosis).
3. **Feeding:** If nutrient is available at the site, take it up and gain *tumour energy from nutrient*.
4. **Division:** With probability *tumour division*, divide into an adjacent empty site; parent energy is halved between the two cells. If no adjacent space is free, division cannot occur (*contact inhibition*).

#### Cytotoxic T Cells
1. **Chemotaxis & Movement:** With probability *chemotaxis*, step toward a neighboring tumour cell; otherwise take a random step.
2. **Metabolism:** Expend 1 unit of energy per step; if energy drops below 0, the cell dies of natural senescence.
3. **Predation (Cell Lysis):** If the T cell is active (not exhausted) and shares a site with a tumour cell, it attacks with probability *kill probability*, gains *energy per kill*, and logs the kill event.
4. **Clonal Proliferation:** Following a successful kill, divide with probability *proliferation*; parent energy is split.
5. **Exhaustion:** When a T cell reaches the *exhaustion limit* of kills, it becomes exhausted: it can no longer kill or divide.

#### Systemic Recruitment
At each step, new T cells infiltrate the tissue from the circulation at a rate proportional to overall tumour burden (antigen load).

#### Checkpoint Inhibitor Therapy (Anti-PD-1)
When therapy is active, the kill probability is multiplied by the *kill boost*, and the exhaustion limit is multiplied by the *exhaustion multiplier*. This immediately reinvigorates exhausted T cells whose kill count falls below the new limit!
"

ui <- fluidPage(
  titlePanel("T Cells vs. Tumour: An Agent-Based Immune-Oncology Model"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("preset", "Clinical Scenario Preset", choices = names(presets)),
      tags$h5(tags$b("Tissue & Tumour Parameters")),
      sliderInput("dims", "Grid size (n x n)", 10, 60, 40),
      sliderInput("n_tumour", "Initial tumour cells", 0, 400, 60, step = 10),
      sliderInput("resupply_time", "Nutrient resupply time", 1, 50, 12),
      sliderInput("de_tumour", "Tumour energy from nutrient", 1, 40, 14),
      sliderInput("tumour_div", "Tumour division prob.", 0, 0.3, 0.06, step = 0.01),
      sliderInput("tumour_motility", "Tumour motility", 0, 0.5, 0.05, step = 0.01),
      tags$h5(tags$b("Immune (T-Cell) Parameters")),
      sliderInput("n_tcell", "Initial T cells", 0, 200, 10, step = 5),
      sliderInput("recruit", "Recruitment rate (per 100 tumour cells)", 0, 15, 1, step = 0.5),
      sliderInput("p_kill", "Base kill probability", 0, 1, 0.5, step = 0.05),
      sliderInput("de_tcell", "T-cell energy per kill", 1, 40, 12),
      sliderInput("tcell_prolif", "Proliferation prob. after kill", 0, 1, 0.4, step = 0.05),
      sliderInput("exhaustion_limit", "Exhaustion limit (kills)", 1, 20, 6),
      sliderInput("chemotaxis", "Chemotaxis toward tumour", 0, 1, 0.6, step = 0.05),
      tags$h5(tags$b("Immunotherapy (Checkpoint Blockade)")),
      checkboxInput("therapy", "Enable Checkpoint Inhibitor (anti-PD-1)", FALSE),
      sliderInput("therapy_start", "Start therapy at step", 0, 1000, 200, step = 10),
      sliderInput("therapy_kill_boost", "Kill probability boost (multiplier)", 1, 3, 1.6, step = 0.1),
      sliderInput("therapy_exh_mult", "Exhaustion limit multiplier", 1, 6, 4, step = 0.5),
      numericInput("seed", "Random seed", 1),
      actionButton("reset", "Reset Model", class = "btn-primary", width = "100%"),
      helpText("Parameter changes take effect when you click Reset Model.")
    ),
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel(
          "Interactive Simulation",
          br(),
          fluidRow(
            column(3, actionButton("play", "Play", width = "100%")),
            column(3, actionButton("step1", "Step Once", width = "100%")),
            column(3, numericInput("n_run", NULL, 200, min = 1, step = 50, width = "100%")),
            column(3, actionButton("run_n", "Run N Steps", width = "100%"))
          ),
          fluidRow(
            column(6, sliderInput("speed", "Frame delay (ms)", 50, 1000, 150, step = 50, width = "100%")),
            column(6, sliderInput("spf", "Steps per frame", 1, 20, 1, width = "100%"))
          ),
          fluidRow(
            column(6, plotlyOutput("grid", height = "520px")),
            column(6, plotlyOutput("ts", height = "520px"))
          )
        ),
        tabPanel("Biological Rules & Formulation", br(), markdown(rules_md))
      )
    )
  )
)

server <- function(input, output, session) {
  rv <- reactiveValues(m = NULL, tick = 0)
  running <- reactiveVal(FALSE)

  set_running <- function(on) {
    running(on)
    updateActionButton(session, "play", label = if (on) "Pause" else "Play")
  }

  observeEvent(input$preset, {
    p <- modifyList(default_params, presets[[input$preset]])
    for (nm in param_names) updateSliderInput(session, nm, value = p[[nm]])
    updateCheckboxInput(session, "therapy", value = p$therapy)
  }, ignoreInit = TRUE)

  reset_model <- function() {
    set_running(FALSE)
    p <- default_params
    for (nm in c(param_names, "therapy", "seed")) {
      if (!is.null(input[[nm]])) p[[nm]] <- input[[nm]]
    }
    rv$m <- initialize_model(p)
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
    withProgress(message = "Simulating steps...", value = 0, {
      for (k in seq_len(n)) {
        model_step(rv$m)
        if (k %% 25 == 0) incProgress(25 / n)
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
