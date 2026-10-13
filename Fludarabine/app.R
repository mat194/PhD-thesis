# Libraries -------------------------------------------------------------------
library(shiny)
library(DT)
library(shinyjs)
library(tidyverse)
library(posologyr)
library(rxode2)
library(plotly)
library(shinydashboard)
library(shinyWidgets)
library(shinyFeedback)
library(shinyalert)

options(scipen = 999)
rxode2::setRxThreads(2L)

source("models.R")

# Only one model is used by this app
MODEL      <- Fludarabine$Langenhorst_CPK2019
MODEL_NAME <- "Langenhorst_CPK2019"


# PK helper functions ---------------------------------------------------------

# eGFR as defined in Langenhorst et al. 2019 (Clin Pharmacokinet 58:627-637),
# mL/min/1.73 m2, from serum creatinine (mg/dL):
#  - adults: Cockcroft-Gault (age >= 17 y women, >= 14 y men), indexed to 1.73 m2
#  - children (< 17 y women, < 14 y men): Schwartz (bedside, 0.413 * height / SCr)
#  - capped at 140; the cap rises from 35 at birth to 140 at 1.5 y and falls by
#    8 per decade after 30 y (the "maximum eGFR" of the paper).
# NOTE: the paper cites its equations without giving the constants; the BSA
# (Du Bois) used for indexing and the Schwartz constant are assumptions.
calculate_egfr <- function(age, serum_creatinine_mg_dl, gender, weight_kg, height_cm) {
  child <- age < if (gender == "male") 14 else 17
  if (child) {
    egfr <- 0.413 * height_cm / serum_creatinine_mg_dl
  } else {
    cg <- (140 - age) * weight_kg / (72 * serum_creatinine_mg_dl)
    if (gender == "female") cg <- cg * 0.85
    bsa  <- 0.007184 * weight_kg^0.425 * height_cm^0.725      # Du Bois
    egfr <- cg * 1.73 / bsa
  }
  cap <- if (age < 1.5) 35 + (140 - 35) * age / 1.5
         else if (age <= 30) 140
         else 140 - 8 * (age - 30) / 10
  min(egfr, cap)
}


# The model uses eGFR in L/h (renal CL = eGFR x 0.782 L/h per 70 kg): the
# size-adjusted value in mL/min/1.73 m2 times 0.06; body size is handled by the
# model's own allometric scaling, so no BSA de-indexing here.
egfr_lh <- function(egfr_indexed) egfr_indexed * 0.06


# Core helpers ----------------------------------------------------------------

# Sampling grid (h after first dose) covering every dose and the 24 h after it.
# Dose times and dose + 24 h are always included so the AUC can be read exactly.
make_time_grid <- function(dose_rel_times, step = 0.5) {
  end <- max(dose_rel_times) + 24
  sort(unique(round(c(seq(0, end, by = step), dose_rel_times, dose_rel_times + 24), 1)))
}

# Event/covariate table for posologyr.
#  - with TDM -> poso_estim_map(): grid rows are dummy doses (EVID = 1, AMT = 0)
#    and TDM rows are observations (EVID = 0, DV = conc)
#  - no TDM   -> poso_simu_pop(): grid rows are blank samples (EVID = 0),
#    the doses are added as separate EVID = 1 rows
build_patient_table <- function(doses_rel, tdm_rel, egfr, bw, with_tdm) {
  grid <- make_time_grid(doses_rel$TIME)

  base <- tibble(
    ID = 1, TIME = grid, EVID = if (with_tdm) 1 else 0, AMT = 0, DUR = 0,
    DV = NA_real_, eGFR = egfr, BW = bw, OCC = 1
  )

  if (with_tdm) {
    base <- base %>%
      rows_update(doses_rel %>% transmute(TIME, AMT, DUR), by = "TIME")
    # a TDM taken exactly at a dosing time is a pre-dose sample: nudge it back
    tdm_rel <- tdm_rel %>%
      mutate(TIME = if_else(TIME %in% doses_rel$TIME, TIME - 0.1, TIME))
    base <- base %>%
      filter(!TIME %in% tdm_rel$TIME) %>%
      bind_rows(tdm_rel %>% transmute(ID = 1, TIME, EVID = 0, AMT = 0, DUR = 0,
                                      DV, eGFR = egfr, BW = bw, OCC = 1)) %>%
      arrange(TIME)
  } else {
    base <- base %>%
      bind_rows(doses_rel %>% transmute(ID = 1, TIME, EVID = 1, AMT, DUR, DV = NA_real_,
                                        eGFR = egfr, BW = bw, OCC = 1)) %>%
      arrange(TIME, EVID)   # at equal times the (blank) sample precedes the dose
  }
  data.table::as.data.table(base)
}

# AUC0-24 after each dose = AUC(t_dose + 24) - AUC(t_dose), read off the solved profile
auc24_per_dose <- function(profile, doses_rel, first_dose_time) {
  auc_at <- function(t) profile$AUC[which.min(abs(profile$time - t))]
  doses_rel %>%
    mutate(
      Dose       = row_number(),
      AUC0_24    = map_dbl(TIME, ~ auc_at(.x + 24) - auc_at(.x)),
      Time       = first_dose_time + as.difftime(TIME, units = "hours")
    )
}

N_SIM <- 1000   # simulated patients for the prediction intervals (no-TDM case)

# 90% prediction intervals (5th-95th percentile) from N_SIM simulated patients:
# concentration at every time, and AUC0-24 after every dose.
simulate_pi <- function(dat, doses_rel, model) {
  sim <- poso_simu_pop(dat, prior_model = model, n_simul = N_SIM)$model %>%
    data.table::as.data.table()
  sim[, time := round(time, 1)]

  conc <- sim[, .(lower = quantile(Cc, 0.05), upper = quantile(Cc, 0.95)), by = time]

  auc_dose <- lapply(doses_rel$TIME, function(t) {
    a0 <- sim[time == t,        .(sim.id, a0 = AUC)]
    a1 <- sim[time == t + 24,   .(sim.id, a1 = AUC)]
    v  <- (a1[a0, on = "sim.id"])[, a1 - a0]
    tibble(lower = quantile(v, 0.05, na.rm = TRUE), upper = quantile(v, 0.95, na.rm = TRUE))
  }) %>% bind_rows()

  list(conc = as_tibble(conc), auc = auc_dose)
}


run_fludarabine <- function(doses, tdm, egfr, bw, model) {
  first_dose_time <- min(doses$TIME)
  doses_rel <- doses %>%
    mutate(TIME = round(as.numeric(difftime(TIME, first_dose_time, units = "hours")), 1)) %>%
    arrange(TIME)
  tdm_rel <- tdm %>%
    mutate(TIME = round(as.numeric(difftime(TIME, first_dose_time, units = "hours")), 1))

  with_tdm <- nrow(tdm) > 0
  dat <- build_patient_table(doses_rel, tdm_rel, egfr, bw, with_tdm)

  if (with_tdm) {
    fit <- poso_estim_map(dat, model)
  } else {
    fit <- poso_simu_pop(dat, prior_model = model, n_simul = 0)
  }

  # prediction intervals: only without TDM (population simulation)
  pi <- if (with_tdm) NULL else {
    set.seed(123)
    simulate_pi(dat, doses_rel, model)
  }

  profile <- data.table::as.data.table(fit$model) %>%
    as_tibble() %>%
    distinct(time, .keep_all = TRUE) %>%
    arrange(time)

  params <- profile %>%
    slice(1) %>%
    select(any_of(c("Cl", "V1", "V2", "V3", "Q2", "Q3"))) %>%
    pivot_longer(everything(), names_to = "Parameter", values_to = "Value") %>%
    mutate(
      Parameter = recode(Parameter, Cl = "CL (L/h)", V1 = "V1 (L)", V2 = "V2 (L)",
                         V3 = "V3 (L)", Q2 = "Q2 (L/h)", Q3 = "Q3 (L/h)"),
      Value = signif(Value, 3)
    )

  list(
    profile         = profile %>%
      { if (is.null(pi)) . else left_join(., pi$conc %>% mutate(time = round(time, 1)), by = "time") } %>%
      mutate(Time = first_dose_time + as.difftime(time, units = "hours")),
    auc             = auc24_per_dose(profile, doses_rel, first_dose_time) %>%
      { if (is.null(pi)) . else bind_cols(., pi$auc) },
    has_pi          = !is.null(pi),
    params          = params,
    tdm             = tdm_rel %>% mutate(Time = first_dose_time + as.difftime(TIME, units = "hours")),
    eta             = if (with_tdm) fit$eta else NULL,   # MAP ETAs (IIV only)
    with_tdm        = with_tdm,
    first_dose_time = first_dose_time
  )
}


# Plots -----------------------------------------------------------------------

build_auc_plot <- function(res) {
  plot_ly(res$auc, x = ~Time, y = ~AUC0_24, type = "bar",
          marker = list(color = "#3c8dbc"),
          hoverinfo = "text",
          text = ~paste0("Dose ", Dose, " (", AMT, " mg)<br>",
                         format(Time, "%d/%m/%Y %H:%M"),
                         "<br>AUC0-24: ", round(AUC0_24, 1), " mg·h/L")) %>%
    layout(
      title = list(text = paste0("AUC0-24 per dose - ", MODEL_NAME,
                                 if (res$with_tdm) " (MAP-Bayesian)" else " (population typical values)")),
      xaxis = list(title = ""),
      yaxis = list(title = "AUC0-24 (mg·h/L)", rangemode = "tozero")
    )
}

build_conc_plot <- function(res) {
  p <- plot_ly()
  if (res$has_pi) {
    p <- p %>% add_ribbons(data = res$profile, x = ~Time, ymin = ~pmax(lower, 1e-4), ymax = ~upper,
                           name = "90% PI", line = list(color = "rgba(0,0,255,0.3)"),
                           fillcolor = "rgba(0,0,255,0.2)", hoverinfo = "text",
                           text = ~paste("90% PI:", round(lower, 3), "-", round(upper, 3)))
  }
  p <- p %>%
    add_lines(data = res$profile, x = ~Time, y = ~Cc, name = "Concentration",
              line = list(color = "blue"), hoverinfo = "text",
              text = ~paste("Date:", format(Time, "%b %d, %H:%M"), "<br>Estimated:", round(Cc, 3)))
  if (nrow(res$tdm) > 0) {
    p <- p %>% add_markers(data = res$tdm, x = ~Time, y = ~DV, name = "TDM",
                           marker = list(size = 12, color = "orange"), hoverinfo = "text",
                           text = ~paste("Measured:", DV))
  }
  p %>% layout(
    title = list(text = "Concentration-time profile"),
    xaxis = list(title = ""),
    yaxis = list(title = "Concentration (mg/L)", type = "log"),
    hovermode = "x unified",
    legend = list(orientation = "h", x = 0.5, y = -0.2, xanchor = "center")
  )
}


# EBE (MAP eta) of this patient against the prior N(0, omega) of each eta
build_eta_plot <- function(eta) {
  eta <- eta[!grepl("^KAPPA_", names(eta))]
  om  <- MODEL$omega[names(eta), names(eta), drop = FALSE]

  tbl <- tibble(eta_name = names(eta), eta = as.numeric(eta), omega_sd = sqrt(diag(om))) %>%
    mutate(
      eta_z   = eta / omega_sd,
      outside = abs(eta_z) > 1.96,
      lim     = pmax(4 * omega_sd, 1.15 * abs(eta)),   # keep an extreme EBE on screen
      peak    = dnorm(0, 0, omega_sd),
      lab     = sprintf("EBE %+.2f (%.1f SD)", eta, eta_z)
    )

  curves <- tbl %>%
    mutate(x = map(lim, ~ seq(-.x, .x, length.out = 400))) %>%
    tidyr::unnest(x) %>%
    mutate(density = dnorm(x, 0, omega_sd))

  ggplot(curves, aes(x, density)) +
    geom_area(fill = "#4575b4", alpha = 0.13) +
    geom_line(colour = "#4575b4", linewidth = 0.8) +
    geom_segment(data = tbl, aes(x = eta, xend = eta, y = 0, yend = peak, colour = outside),
                 linewidth = 1, inherit.aes = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35") +
    geom_label(data = tbl, aes(x = lim, y = peak, label = lab), hjust = 1, vjust = 1, size = 3.3,
               fill = "white", alpha = 0.75, linewidth = 0, inherit.aes = FALSE) +
    scale_colour_manual(
      values = c(`FALSE` = "#1b7837", `TRUE` = "#b2182b"),
      labels = c(`FALSE` = "within the prior's central 95%", `TRUE` = "outside it"),
      name = NULL, drop = FALSE) +
    facet_wrap(~ eta_name, scales = "free") +
    labs(title = "Individual MAP etas (EBE) against the prior distribution",
         subtitle = "Curve = N(0, omega) declared by the model; vertical line = this patient's EBE",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "top")
}


# UI helper -------------------------------------------------------------------

addRemoveButtons <- function(add_id, remove_id, width) {
  column(width, br(),
         actionButton(add_id, "+", class = "btn btn-success btn-xs"),
         actionButton(remove_id, "-", class = "btn btn-danger btn-xs"))
}


# UI --------------------------------------------------------------------------

ui <- dashboardPage(
  dashboardHeader(title = tagList(span(class = "logo-lg", "Fludarabine"), icon("pills", class = "fa-lg"))),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Info",      tabName = "patient_info",      icon = icon("user")),
      menuItem("Dose Observations", tabName = "dose_observations", icon = icon("pills")),
      hr()
    )
  ),

  dashboardBody(
    useShinyjs(),
    useShinyFeedback(),
    # keep-alive: a small message every 20 s so the proxy never sees an idle connection
    tags$script(HTML("setInterval(function() { if (window.Shiny && Shiny.setInputValue) Shiny.setInputValue('keepalive', Date.now()); }, 20000);")),

    tabItems(
      tabItem(tabName = "patient_info",
        fluidRow(
          column(width = 8,
            box(width = 12, solidHeader = TRUE, status = "primary",
              fluidRow(
                column(4, numericInput("age",    "Age (years)",  value = 40)),
                column(4, selectInput("gender", "Gender", choices = c("male", "female"))),
                column(4, numericInput("weight", "Weight (kg)", value = 70))
              ),
              fluidRow(column(4, numericInput("height", "Height (cm)", value = 170))),
              h4("Creatinine"),
              fluidRow(
                column(4, numericInput("creatinine", "Serum Cr (mg/dL)", value = 0.9, step = 0.1)),
                column(8, br(), uiOutput("egfr_text"))
              ),
              h4("TDM (optional)"),
              p("If no TDM is entered, the AUC is computed with the population typical values.",
                style = "color:#666;"),
              fluidRow(
                column(4, airDatepickerInput("tdm_time", "Day/time of record", value = Sys.time(),
                  timepicker = TRUE, timepickerOpts = timepickerOptions(hoursStep = 1, minutesStep = 10), width = "100%")),
                column(4, numericInput("tdm", "F-ara-A (mg/L)", value = 0, min = 0)),
                addRemoveButtons("addTDM", "removeTDM", 4)
              )
            )
          ),
          column(width = 4,
            box(title = "TDM", width = 12, solidHeader = TRUE, status = "primary", DTOutput("TableTDM"))
          )
        )
      ),

      tabItem(tabName = "dose_observations",
        fluidRow(
          column(width = 6,
            box(width = 12, solidHeader = TRUE, status = "primary",
              airDatepickerInput("dose_time", "Day/time of first administration", value = Sys.time(),
                timepicker = TRUE, timepickerOpts = timepickerOptions(hoursStep = 1, minutesStep = 5), width = "100%"),
              numericInput("dose", "Dose per administration (mg)", value = 50, min = 0),
              numericInput("dur", "Infusion duration (hours)", value = 0.5, min = 0.1, step = 0.1),
              numericInput("n_days", "Number of consecutive daily doses (q24h)", value = 1, min = 1, step = 1),
              fluidRow(addRemoveButtons("addDose", "removeDose", 6))
            )
          ),
          column(width = 6,
            box(width = 12, solidHeader = TRUE, status = "primary", DTOutput("doseTable", width = "75%"))
          )
        )
      )
    ),

    fluidRow(
      column(width = 12,
        box(width = 12, solidHeader = TRUE, status = "primary", title = "",
          actionButton("calculate_button", "Run calculation", class = "btn btn-primary btn-lg"),
          br(), br(),
          uiOutput("message"),
          plotlyOutput("aucPlot", height = "400px"),
          br(),
          plotlyOutput("concentrationPlot", height = "400px"),
          br(),
          uiOutput("box_eta"),
          uiOutput("box_results")
        )
      )
    ),
    scrollToTop = TRUE
  )
)


# Server ----------------------------------------------------------------------

server <- function(input, output, session) {

  showModal(modalDialog(
    title = "Fludarabine Precision Dosing Tool",
    p("Disclaimer: This tool is for informational purposes only and intended for use only by health care professionals. The tool is not intended to be a substitute for professional medical advice, dosing, diagnosis or treatment"),
    footer = tagList(fluidRow(
      column(3, checkboxInput("accept", "I understand", value = FALSE)),
      column(9, actionButton("confirm", "Confirm", class = "btn-primary"))
    )),
    easyClose = FALSE, size = "m"
  ))
  shinyjs::disable("confirm")
  observeEvent(input$accept, shinyjs::toggleState("confirm", isTRUE(input$accept)))
  observeEvent(input$confirm, removeModal())

  rv <- reactiveValues(
    doses = tibble(TIME = as.POSIXct(character()), AMT = numeric(), DUR = numeric()),
    tdm   = tibble(TIME = as.POSIXct(character()), DV = numeric())
  )

  # Input validation
  observe({
    for (id in c("age", "weight", "height", "creatinine")) {
      hideFeedback(id)
      v <- input[[id]]
      if (is.null(v) || is.na(v) || v <= 0) showFeedbackDanger(id, "A positive value is required")
    }
  })

  egfr <- reactive({
    req(isTruthy(input$age), isTruthy(input$creatinine), isTruthy(input$weight), isTruthy(input$height))
    calculate_egfr(input$age, input$creatinine, input$gender, input$weight, input$height)
  })
  egfr_model <- reactive(egfr_lh(egfr()))
  output$egfr_text <- renderUI(p(HTML(sprintf(
    "eGFR: <strong>%.1f</strong> mL/min/1.73 m&sup2; = <strong>%.2f</strong> L/h (as used by the model)",
    egfr(), egfr_model()))))

  # TDM
  observeEvent(input$addTDM, {
    new_entry <- tibble(TIME = as.POSIXct(input$tdm_time, format = "%d/%m/%Y %H:%M"), DV = input$tdm)
    if (!any(rv$tdm$TIME == new_entry$TIME)) rv$tdm <- bind_rows(rv$tdm, new_entry) %>% arrange(TIME)
  })
  observeEvent(input$removeTDM, if (nrow(rv$tdm) > 0) rv$tdm <- rv$tdm[-nrow(rv$tdm), ])

  # Doses: n_days consecutive q24h administrations starting at dose_time
  observeEvent(input$addDose, {
    req(isTruthy(input$dose), isTruthy(input$dur), isTruthy(input$n_days))
    t0 <- as.POSIXct(input$dose_time, format = "%d/%m/%Y %H:%M")
    new_doses <- tibble(
      TIME = t0 + as.difftime(24 * (seq_len(input$n_days) - 1), units = "hours"),
      AMT = input$dose, DUR = input$dur
    ) %>% filter(!TIME %in% rv$doses$TIME)
    rv$doses <- bind_rows(rv$doses, new_doses) %>% arrange(TIME)
  })
  observeEvent(input$removeDose, if (nrow(rv$doses) > 0) rv$doses <- rv$doses[-nrow(rv$doses), ])

  output$doseTable <- renderDT({
    rv$doses %>% mutate(TIME = format(TIME, "%d/%m/%Y %H:%M")) %>% datatable(options = list(dom = "t"))
  }, server = FALSE)

  output$TableTDM <- renderDT({
    rv$tdm %>% mutate(TIME = format(TIME, "%d/%m/%Y %H:%M")) %>% datatable(options = list(dom = "t"))
  }, server = FALSE)

  # Main calculation
  observeEvent(input$calculate_button, {
    doses <- isolate(rv$doses)
    tdm   <- isolate(rv$tdm)

    if (nrow(doses) == 0)                        { shinyalert("Error", "No doses data are provided", type = "error"); return(NULL) }
    if (!all(c(input$weight, input$age, input$creatinine, input$height) > 0, na.rm = FALSE) ||
        anyNA(c(input$weight, input$age, input$creatinine, input$height))) {
      shinyalert("Error", "Age, weight, height and creatinine must be positive numbers", type = "error"); return(NULL)
    }
    if (nrow(tdm) > 0 && min(tdm$TIME) <= min(doses$TIME)) { shinyalert("Error", "TDM data must be after the first dose", type = "error"); return(NULL) }

    shinybusy::show_modal_spinner(text = "Running estimation...")
    on.exit(shinybusy::remove_modal_spinner())

    res <- tryCatch(
      run_fludarabine(doses, tdm, egfr_model(), input$weight, MODEL),
      error = function(e) { shinyalert("Error", paste("Estimation failed:", conditionMessage(e)), type = "error"); NULL }
    )
    req(res)

    output$aucPlot           <- renderPlotly(build_auc_plot(res))
    output$concentrationPlot <- renderPlotly(build_conc_plot(res))

    output$message <- renderUI(p(
      if (res$with_tdm) "AUC0-24 after each dose was derived from the MAP-Bayesian individual profile (conditioned on the TDM)."
      else sprintf("No TDM provided: AUC0-24 after each dose was derived from the population typical values. The concentration plot band and the AUC table show the 90%% prediction interval from %d simulated patients.", N_SIM),
      " All computations were performed with the posologyr package.", style = "color:#666;"
    ))

    # EBE vs prior: only meaningful when the MAP estimation was run (TDM available)
    output$etaPlot <- renderPlot(if (res$with_tdm) build_eta_plot(res$eta))
    output$box_eta <- renderUI(
      if (res$with_tdm) {
        box(title = tagList(icon("chart-area"), " EBE vs prior distribution"), status = "primary",
            solidHeader = TRUE, width = 12, plotOutput("etaPlot", height = "320px"))
      }
    )

    output$aucTable <- renderDT({
      res$auc %>%
        transmute(Dose, `Date/time` = format(Time, "%d/%m/%Y %H:%M"), `Dose (mg)` = AMT,
                  `AUC0-24 (mg·h/L)` = round(AUC0_24, 1),
                  `90% PI` = if (res$has_pi) sprintf("%.1f - %.1f", lower, upper) else NA_character_) %>%
        { if (res$has_pi) . else select(., -`90% PI`) } %>%
        datatable(options = list(dom = "t"), rownames = FALSE)
    })
    output$paramTable <- renderDT({
      datatable(res$params, options = list(dom = "t"), rownames = FALSE)
    })

    output$box_results <- renderUI(fluidRow(
      column(6, box(title = tagList(icon("calculator"), " AUC0-24 per dose"), status = "primary",
                    solidHeader = TRUE, width = 12, DTOutput("aucTable"))),
      column(6, box(title = tagList(icon("notes-medical"), " Individual PK parameters"), status = "primary",
                    solidHeader = TRUE, width = 12,
                    p(sprintf("eGFR %.1f mL/min/1.73 m², body weight %.0f kg. %s", egfr(), input$weight,
                              if (res$with_tdm) "MAP Bayesian estimates." else "Population typical values (no TDM)."),
                      style = "color:#666;font-size:13px;"),
                    DTOutput("paramTable")))
    ))
  })
}

shinyApp(ui = ui, server = server)
