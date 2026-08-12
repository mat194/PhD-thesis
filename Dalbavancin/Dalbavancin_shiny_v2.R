# Libraries -------------------------------------------------------------------
library(shiny)
library(DT)
library(shinyjs)
library(tidyverse)
library(posologyr)
library(lubridate)
library(lotri)
library(rxode2)
library(plotly)
library(shinycssloaders)
library(shinydashboard)
library(shinyWidgets)
library(shinydashboardPlus)
library(shinyFeedback)
library(shinyBS)
library(waiter)
library(shinyalert)

options(scipen = 999)
rxode2::setRxThreads(2L)

source("models.R")


# PK helper functions ---------------------------------------------------------

calculate_bsa <- function(weight_kg, height_cm) {
  0.007184 * (weight_kg^0.425) * (height_cm^0.725)
}

calculate_bsa_mosteller <- function(weight_kg, height_cm) {
  sqrt((weight_kg * height_cm) / 3600)
}

calculate_abw <- function(height_cm, weight, gender) {
  height_in <- (height_cm / 100) * 39.37
  IBW <- if (gender == "male") 50 + 2.3 * (height_in - 60) else 45.5 + 2.3 * (height_in - 60)
  IBW + 0.4 * (weight - IBW)
}

calculate_crcl <- function(age, weight_kg, serum_creatinine_mg_dl, gender) {
  CrCl <- (140 - age) * weight_kg / (72 * serum_creatinine_mg_dl)
  if (gender == "female") CrCl <- CrCl * 0.85
  CrCl
}

calculate_egfr <- function(age, serum_creatinine_mg_dl, gender) {
  if (gender == "male") {
    kappa <- 0.9; alpha <- -0.302
    egfr <- 142 * (min(serum_creatinine_mg_dl / kappa, 1)^alpha) *
      (max(serum_creatinine_mg_dl / kappa, 1)^-1.200) * (0.9938^age)
  } else {
    kappa <- 0.7; alpha <- -0.241
    egfr <- 142 * (min(serum_creatinine_mg_dl / kappa, 1)^alpha) *
      (max(serum_creatinine_mg_dl / kappa, 1)^-1.200) * (0.9938^age) * 1.012
  }
  egfr
}

calculate_metrics <- function(data, model_name) {
  N         <- nrow(data)
  predicted <- data$Cc_tot
  observed  <- data$Cc_observed
  weight    <- data$weight
  rBias <- (1 / N) * sum((predicted - observed) / observed) * 100
  rRMSE <- sqrt((1 / N) * sum(((predicted - observed)^2) / (observed^2))) * 100
  tibble(model_name = model_name, rBias = rBias, rRMSE = rRMSE, weight = weight)
}


# UI helper components --------------------------------------------------------

dalbaAddRemoveButtons <- function(add_id, remove_id, width) {
  column(
    width,
    br(),
    actionButton(add_id, "+", class = "btn btn-success btn-xs"),
    actionButton(remove_id, "-", class = "btn btn-danger btn-xs")
  )
}

dalbaModelSpecTab <- function(model_names) {
  tabItem(
    tabName = "model_specification",
    box(
      width = 12, solidHeader = TRUE, status = "primary", collapsible = FALSE,
      fluidRow(
        column(4,
          tags$label("Models included in the automated model selection/averaging"),
          checkboxGroupButtons(
            inputId  = "model_selection",
            label    = NULL,
            selected = model_names,
            choices  = sort(model_names),
            size      = "lg",
            direction = "vertical",
            justified = TRUE,
            checkIcon = list(yes = icon("ok", lib = "glyphicon"))
          )
        ),
        column(8,
          tags$p("Select one or more models to include in the automated model selection or model averaging algorithm.")
        )
      )
    )
  )
}


# Generic add/remove observer -------------------------------------------------

dalbaRegisterAddRemove <- function(input, output, session,
                                   input_add, input_remove,
                                   input_time, input_value,
                                   var_name, var_col, rv,
                                   time_format = "%d/%m/%Y %H:%M") {
  observeEvent(input[[input_add]], {
    new_entry <- tibble(
      TIME = as.POSIXct(input[[input_time]], format = time_format),
      !!var_col := input[[input_value]]
    ) %>%
      distinct(TIME, .keep_all = TRUE)

    if (!any(rv[[var_name]]$TIME == new_entry$TIME)) {
      rv[[var_name]] <- bind_rows(rv[[var_name]], new_entry)
    }
  })

  observeEvent(input[[input_remove]], {
    if (nrow(rv[[var_name]]) > 0) {
      rv[[var_name]] <- rv[[var_name]][-nrow(rv[[var_name]]), ]
    }
  })
}


# Core processing helpers -----------------------------------------------------

run_map_estimation <- function(patient_data_for_map2, model_list, free_drug) {
  model_list_names <- as.list(names(model_list))
  m_out      <- list()
  error_log  <- list()

  for (i in seq_along(model_list)) {
    tryCatch({
      m_i      <- poso_estim_map(patient_data_for_map2, model_list[[i]], return_model = TRUE, return_ofv = TRUE)
      model_dt <- data.table::data.table(m_i$model)
      m_out[[i]] <- cbind(
        model_dt[, .(time, Cc, AUC)],
        Cc_tot     = if ("Cc_tot" %in% names(model_dt)) model_dt$Cc_tot else model_dt$Cc,
        ofv        = m_i$ofv,
        LL         = exp(-0.5 * m_i$ofv),
        model_name = model_list_names[[i]],
        model_number = i
      )
    }, error = function(e) {
      error_log[[length(error_log) + 1]] <<- list(
        model_name   = model_list_names[[i]],
        model_number = i,
        error_message = e$message
      )
    })
  }

  m_dt_out <- data.table::as.data.table(do.call(rbind, m_out))
  m_dt_out[, weight := LL / sum(unique(LL))]

  free_auc_models <- names(model_list)[sapply(model_list, function(m) isTRUE(m$auc_is_free))]
  m_dt_out[, AUC_free := data.table::fifelse(model_name %in% free_auc_models, AUC, AUC * free_drug)]

  m_dt_out[, MoA_IPRED := sum(Cc       * weight), by = time]
  m_dt_out[, MoA_AUC   := sum(AUC_free * weight), by = time]

  m_dt_out
}


run_sir_estimation <- function(patient_data_for_map, model_list, best_model, time_seq) {
  sir_map    <- poso_estim_sir(patient_data_for_map, model_list[[best_model]], n_sample = 1e3, n_resample = 1e2)
  event_table <- rxode2::as.et(patient_data_for_map)
  event_table$add_sampling(seq(0, time_seq, by = 1))

  ci_map <- poso_replace_et(
    target_model  = sir_map$model,
    prior_model   = model_list[[best_model]],
    event_table   = event_table,
    interpolation = "locf"
  )

  ci_map %>%
    group_by(time) %>%
    summarise(
      lower90 = round(quantile(Cc, probs = 0.05, na.rm = TRUE)),
      upper90 = round(quantile(Cc, probs = 0.95, na.rm = TRUE))
    ) %>%
    ungroup()
}


compute_auc_targets <- function(prova, target_mic, target_auc, from) {
  auc_data_plot <- prova %>%
    mutate(AUC = AUC / target_mic) %>%
    filter(time %% 24 == 0) %>%
    mutate(AUC24 = round(AUC - lag(AUC, default = 0))) %>%
    mutate(above_target = if_else(AUC24 > target_auc, "Above", "Below"))

  auc_data <- auc_data_plot %>% filter(time >= from)

  first_below_target_ab <- auc_data %>%
    filter(above_target == "Above") %>%
    tail(1)

  list(auc_data_plot = auc_data_plot, first_below_target_ab = first_below_target_ab)
}


build_weight_plot <- function(metrics_long) {
  plot_ly(
    data      = metrics_long %>% dplyr::filter(metric == "weight"),
    x         = ~"Weight",
    y         = ~value,
    color     = ~model_name,
    colors    = "Set2",
    type      = "bar",
    hoverinfo = "text",
    showlegend = FALSE,
    text      = ~paste(model_name, ":", round(value, 2))
  ) %>%
    layout(
      barmode = "stack",
      yaxis   = list(title = "", range = c(0, 1.01), side = "right"),
      xaxis   = list(title = ""),
      showlegend = FALSE
    )
}


build_moa_plot <- function(data_prova, data_plot_dv, p_weights) {
  plot_moa <- plot_ly() %>%
    add_lines(
      data      = data_prova,
      x         = ~Time,
      y         = ~Cc,
      color     = ~model_name,
      name      = ~model_name,
      line      = list(width = 2),
      hoverinfo = "text",
      text      = ~paste("Date:", format(Time, "%b %d, %H:%M"), "<br>Estimated:", Cc, "<br>Model:", model_name)
    ) %>%
    add_lines(
      data      = data_prova,
      x         = ~Time,
      y         = ~MoA_IPRED,
      name      = "Model Average",
      line      = list(width = 3, dash = "solid", color = "black"),
      hoverinfo = "text",
      text      = ~paste("Date:", format(Time, "%b %d, %H:%M"), "<br>Model Averaging:", MoA_IPRED)
    ) %>%
    add_markers(
      data      = data_plot_dv,
      x         = ~Time,
      y         = ~value,
      name      = "Patient Data",
      marker    = list(size = 12, color = "orange"),
      hoverinfo = "text",
      text      = ~paste("Measured:", value),
      showlegend = TRUE
    ) %>%
    layout(
      title      = list(text = "Model Averaging"),
      yaxis      = list(title = "Concentration (mg/L)", autorange = TRUE, type = "log"),
      xaxis      = list(title = "", tickformat = "%b%d"),
      hovermode  = "x unified",
      hoverlabel = list(font = list(size = 12)),
      showlegend = TRUE,
      legend     = list(orientation = "h", x = 0.5, y = -0.2, xanchor = "center", yanchor = "bottom")
    )

  subplot(
    plot_moa, p_weights,
    widths = c(0.75, 0.25),
    margin = 0.02,
    titleX = TRUE, titleY = TRUE
  ) %>%
    layout(
      showlegend = TRUE,
      legend     = list(orientation = "h", x = 0.5, y = -0.2, xanchor = "center", yanchor = "bottom", hovermode = "x unified")
    )
}


build_single_model_plot <- function(data_plot, data_plot_map, data_plot_dv, best_model) {
  plot_ly() %>%
    add_lines(
      data      = data_plot,
      x         = ~Time, y = ~Cocentration, name = "Concentration",
      line      = list(color = "blue"),
      hoverinfo = "text",
      text      = ~paste("Date:", format(Time, "%b %d, %H:%M"), "<br>Estimated:", Cocentration),
      showlegend = TRUE
    ) %>%
    add_ribbons(
      data      = data_plot_map,
      x         = ~Time, ymin = ~lower90, ymax = ~upper90, name = "90% PI",
      line      = list(color = "rgba(0,0,255,0.3)"),
      fillcolor = "rgba(0, 0, 255, 0.3)",
      hoverinfo = "text",
      text      = ~paste("Upper 90 PI:", upper90, "<br>Lower 90 PI:", lower90),
      showlegend = TRUE
    ) %>%
    add_markers(
      data      = data_plot_dv,
      x         = ~Time, y = ~value, name = "Patient Data",
      marker    = list(size = 12, color = "orange"),
      hoverinfo = "text",
      text      = ~paste("Measured:", value, "<br>Error MAP:", Error, "%"),
      showlegend = TRUE
    ) %>%
    layout(
      title      = list(text = paste0("Single model, ", best_model)),
      yaxis      = list(title = "Concentration (mg/L)", autorange = TRUE, type = "log"),
      xaxis      = list(title = "", tickformat = "%b%d"),
      hovermode  = "x unified",
      hoverlabel = list(font = list(size = 12)),
      legend     = list(orientation = "h", x = 0.5, y = -0.2, xanchor = "center", yanchor = "bottom")
    )
}


# UI --------------------------------------------------------------------------

ui <- dashboardPage(
  dashboardHeader(title = tagList(
    span(class = "logo-lg", "Dalbavancin"),
    icon("pills", class = "fa-lg")
  )),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Info",       tabName = "patient_info",        icon = icon("user")),
      menuItem("Dose Observations",  tabName = "dose_observations",   icon = icon("pills")),
      menuItem("Model Targets",      tabName = "model_targets",       icon = icon("bullseye")),
      menuItem("Model Specification",tabName = "model_specification", icon = icon("cogs")),
      menuItem("Report/Load/Save",   tabName = "report",              icon = icon("file-pdf")),
      hr()
    )
  ),

  footer = dashboardFooter(
    tags$p("To cite this application, please use the following reference:")
  ),

  dashboardBody(
    useShinyjs(),
    useShinyFeedback(),

    tabItems(

      # Patient Info Tab -------------------------------------------------------
      tabItem(tabName = "patient_info",
        fluidRow(
          column(width = 8,
            box(width = 12, solidHeader = TRUE, status = "primary",
              fluidRow(
                column(4, numericInput("age",    "Age (years)",  value = 70)),
                column(4, selectInput("gender", "Gender", choices = c("male", "female"), selected = "male")),
                column(4, numericInput("height", "Height (cm)", value = 170))
              ),
              h4("Weight"),
              fluidRow(
                column(4, airDatepickerInput("weight_time",     "Day of record",     value = Sys.Date(), width = "100%")),
                column(4, numericInput("weight",                "Weight (kg)",       value = 70)),
                dalbaAddRemoveButtons("addWeight", "removeWeight", 4)
              ),
              h4("Creatinine"),
              fluidRow(
                column(4, airDatepickerInput("creatinine_time", "Day of record",     value = Sys.Date(), width = "100%")),
                column(4, numericInput("creatinine",            "Serum Cr (mg/dL)", value = 0.9)),
                dalbaAddRemoveButtons("addCreatinine", "removeCreatinine", 4)
              ),
              h4("Albumin"),
              fluidRow(
                column(4, airDatepickerInput("albumin_time",    "Day of record",     value = Sys.Date(), width = "100%")),
                column(4, numericInput("albumin",               "Albumin (g/dL)",   value = 4)),
                dalbaAddRemoveButtons("addAlbumin", "removeAlbumin", 4)
              ),
              h4("TDM"),
              fluidRow(
                column(4, airDatepickerInput("tdm_time", "Day/time of record", value = Sys.time(),
                  timepicker = TRUE, timepickerOpts = timepickerOptions(hoursStep = 1, minutesStep = 10), width = "100%")),
                column(4, numericInput("tdm", "TDM (mg/L)", value = 0)),
                dalbaAddRemoveButtons("addTDM", "removeTDM", 4)
              )
            )
          ),
          column(width = 4,
            box(title = "Patient Data", width = 12, solidHeader = TRUE, status = "primary", collapsible = TRUE,
              DTOutput("TableWT"),    hr(),
              DTOutput("TableCreatinine"), hr(),
              DTOutput("TableAlbumin"),    hr(),
              DTOutput("TableTDM")
            )
          )
        )
      ),

      # Dose Observations Tab --------------------------------------------------
      tabItem(tabName = "dose_observations",
        fluidRow(
          column(width = 6,
            box(width = 12, solidHeader = TRUE, status = "primary",
              airDatepickerInput("dose_time", "Day of Administration", value = Sys.time(),
                timepicker = TRUE, timepickerOpts = timepickerOptions(hoursStep = 1, minutesStep = 5), width = "100%"),
              pickerInput("dose", "Dose (mg)", choices = c(500, 1000, 1500), selected = 1500),
              numericInput("dur", "Infusion rate (hours)", value = 0.5),
              fluidRow(dalbaAddRemoveButtons("addDose", "removeDose", 6))
            )
          ),
          column(width = 6,
            box(width = 12, solidHeader = TRUE, status = "primary",
              DTOutput("doseTable", width = "75%")
            )
          )
        )
      ),

      # Model Targets Tab ------------------------------------------------------
      tabItem(tabName = "model_targets",
        fluidRow(
          column(width = 12,
            box(width = 12, solidHeader = TRUE, status = "primary",
              fluidRow(
                column(3,
                  numericInput("target_auc",             "Targeted last 24h fAUC0-24h/MIC", value = 111),
                  numericInput("target_mic",             "Target MIC",                       value = 0.125),
                  numericInput("target_protein_binding", "Protein binding (%)",              value = 93),
                  actionButton("reset_target", "Reset", class = "btn-primary", style = "width: 100%; color:white;")
                ),
                column(9,
                  tags$div(
                    tags$h4("Explanation of Targets"),
                    tags$p("The 'Targeted last 24h fAUC0-24h/MIC' represents the free drug exposure over the last 24 hours relative to the MIC."),
                    tags$p("According to EUCAST, S. aureus, E. faecalis and E. faecium have MICs/ECOFF of 0.125 mg/L for dalbavancin.")
                  )
                )
              )
            )
          )
        )
      ),

      # Model Specification Tab ------------------------------------------------
      tabItem(tabName = "model_specification",
        dalbaModelSpecTab(model_names = names(Dalbavancin))
      ),

      # Report Tab -------------------------------------------------------------
      tabItem(tabName = "report",
        box(title = "Report Management", width = 12, solidHeader = TRUE, status = "primary",
          textInput("report_name", "Report Name", value = "Dalbavancin Report"),
          fluidRow(
            column(4, downloadButton("save_report", "Save Report",
              icon = icon("save"), class = "btn-primary", style = "width: 100%; color:white;")),
            column(4, downloadButton("save_state",  "Save State",
              icon = icon("save"), class = "btn-primary", style = "width: 100%; color:white;")),
            column(4, fileInput("load_state", label = NULL, buttonLabel = "Load State",
              accept = c(".zip"), width = "100%", placeholder = "Choose .rds file"))
          )
        )
      )
    ),

    # Results row (shown on all tabs after calculation) -----------------------
    fluidRow(
      column(width = 12,
        box(width = 12, solidHeader = TRUE, status = "primary", collapsible = FALSE, title = "",
          actionButton("calculate_button", "Run calculation", class = "btn btn-primary btn-lg"),
          br(), br(),
          uiOutput("message"),
          br(),
          plotlyOutput("concentrationPlot",  height = "400px"),
          br(),
          plotlyOutput("concentrationPlot2", height = "400px"),
          br(),
          uiOutput("box_auc24"),
          br(),
          uiOutput("box_parameters")
        )
      )
    ),

    scrollToTop = TRUE,
    preloader  = list(html = spin_1(), color = "gray")
  )
)


# Server ----------------------------------------------------------------------

server <- function(input, output, session) {

  # Disclaimer modal ----------------------------------------------------------

  showModal(modalDialog(
    title = "Dalbavancin Precision Dosing Tool",
    p("Disclaimer: This tool is for informational purposes only and intended for use only by health care professionals. The tool is not intended to be a substitute for professional medical advice, dosing, diagnosis or treatment"),
    footer = tagList(
      fluidRow(
        column(3, checkboxInput("accept",  "I understand",  value = FALSE)),
        column(9, actionButton("confirm", "Confirm", class = "btn-primary"))
      )
    ),
    easyClose = FALSE,
    size = "m"
  ))

  shinyjs::disable("confirm")

  observeEvent(input$accept, {
    if (input$accept) shinyjs::enable("confirm") else shinyjs::disable("confirm")
  })

  observeEvent(input$confirm, removeModal())


  # Reactive values -----------------------------------------------------------

  rv <- reactiveValues(
    weights     = tibble(TIME = as.POSIXct(character()), WT        = numeric()),
    creatinines = tibble(TIME = as.POSIXct(character()), CREATININE = numeric()),
    albumines   = tibble(TIME = as.POSIXct(character()), ALBUMIN    = numeric()),
    doses       = tibble(TIME = as.POSIXct(character()), AMT = numeric(), DUR = numeric(), ADDL = numeric(), II = numeric()),
    tdm         = tibble(TIME = as.POSIXct(character()), DV  = numeric(), DVID = character())
  )


  # Save / Load state ---------------------------------------------------------

  save_state_to_file <- function(file) {
    state <- list(
      patient_data           = patient_data()$absolute,
      doses                  = rv$doses,
      weights                = rv$weights,
      creatinines            = rv$creatinines,
      albumines              = rv$albumines,
      tdm                    = rv$tdm,
      target_auc             = input$target_auc,
      target_cum_auc         = input$target_cum_auc,
      target_mic             = input$target_mic,
      target_protein_binding = input$target_protein_binding
    )
    temp_file <- tempfile(fileext = ".rds")
    saveRDS(state, temp_file)
    zip(file, temp_file, flags = "-j")
    showNotification("State saved successfully", type = "message")
  }

  load_state_from_file <- function(file) {
    temp_dir <- file.path(tempdir(), "unzip_temp")
    dir.create(temp_dir, showWarnings = FALSE)
    unzip(file, exdir = temp_dir)
    rds_file <- list.files(temp_dir, pattern = "\\.rds$", full.names = TRUE)

    if (length(rds_file) == 1) {
      state <- readRDS(rds_file)
      rv$weights     <- state$weights
      rv$creatinines <- state$creatinines
      rv$doses       <- state$doses
      rv$tdm         <- state$tdm
      rv$albumines   <- state$albumines
      updateNumericInput(session, "target_auc",             value = state$target_auc)
      updateNumericInput(session, "target_cum_auc",         value = state$target_cum_auc)
      updateNumericInput(session, "target_mic",             value = state$target_mic)
      updateNumericInput(session, "target_protein_binding", value = state$target_protein_binding)
      showNotification("State loaded successfully", type = "message")
    } else {
      showNotification("Error: Could not find the .rds file in the zip archive", type = "error")
    }
    unlink(temp_dir, recursive = TRUE)
  }

  output$save_state <- downloadHandler(
    filename = function() paste0("state_", Sys.Date(), ".zip"),
    content  = function(file) save_state_to_file(file)
  )

  observeEvent(input$load_state, {
    req(input$load_state)
    load_state_from_file(input$load_state$datapath)
  })


  # Model selection -----------------------------------------------------------

  updated_model <- reactive({
    Dalbavancin[names(Dalbavancin) %in% input$model_selection]
  })


  # Target reset --------------------------------------------------------------

  observeEvent(input$reset_target, {
    updateNumericInput(session, "target_auc",             value = 111)
    updateNumericInput(session, "target_mic",             value = 0.125)
    updateNumericInput(session, "target_protein_binding", value = 93)
  })


  # Input validation ----------------------------------------------------------

  validate_input <- function(input_id, message) {
    observe({
      hideFeedback(input_id)
      if (is.null(input[[input_id]]) || is.na(input[[input_id]]) || input[[input_id]] == "") {
        showFeedbackDanger(input_id, message)
      }
    })
  }

  validate_input("target_auc",             "Target free-AUC is required")
  validate_input("target_mic",             "Target MIC is required")
  validate_input("target_protein_binding", "Protein binding is required")
  validate_input("age",                    "Age is required")
  validate_input("height",                 "Height is required")


  # Add/Remove observers for covariates ---------------------------------------

  covariate_handlers <- list(
    list(add = "addWeight",    remove = "removeWeight",    time = "weight_time",    value = "weight",    var = "weights",     col = "WT"),
    list(add = "addCreatinine",remove = "removeCreatinine",time = "creatinine_time",value = "creatinine",var = "creatinines", col = "CREATININE"),
    list(add = "addAlbumin",   remove = "removeAlbumin",   time = "albumin_time",   value = "albumin",   var = "albumines",   col = "ALBUMIN")
  )

  purrr::walk(covariate_handlers, function(cfg) {
    dalbaRegisterAddRemove(
      input, output, session,
      input_add    = cfg$add,
      input_remove = cfg$remove,
      input_time   = cfg$time,
      input_value  = cfg$value,
      var_name     = cfg$var,
      var_col      = cfg$col,
      rv           = rv
    )
  })

  # TDM (separate: needs DVID column)
  observeEvent(input$addTDM, {
    new_entry <- tibble(
      TIME = as.POSIXct(input$tdm_time, format = "%d/%m/%Y %H:%M"),
      DV   = input$tdm,
      DVID = "Cc"
    ) %>% distinct(TIME, .keep_all = TRUE)

    if (!any(rv$tdm$TIME == new_entry$TIME)) {
      rv$tdm <- bind_rows(rv$tdm, new_entry)
    }
  })

  observeEvent(input$removeTDM, {
    if (nrow(rv$tdm) > 0) rv$tdm <- rv$tdm[-nrow(rv$tdm), ]
  })


  # Dose add/remove -----------------------------------------------------------

  observeEvent(input$addDose, {
    new_dose <- tibble(
      TIME = as.POSIXct(input$dose_time, format = "%d/%m/%Y %H:%M"),
      AMT  = as.numeric(input$dose),
      DUR  = input$dur,
      ADDL = 0,
      II   = 0
    ) %>%
      distinct(TIME, .keep_all = TRUE) %>%
      arrange(TIME)

    if (!any(rv$doses$TIME == new_dose$TIME)) {
      rv$doses <- bind_rows(rv$doses, new_dose)
    }
  })

  observeEvent(input$removeDose, {
    if (nrow(rv$doses) > 0) rv$doses <- rv$doses[-nrow(rv$doses), ]
  })


  # Dose info reactive --------------------------------------------------------

  dose_info <- reactive({
    if (nrow(rv$doses) == 0) {
      return(list(
        dose_expanded      = tibble(),
        first_dose_time    = NA,
        last_dose_time     = NA,
        last_dur           = NA,
        dose_objects       = tibble(),
        max_time_seq       = NA,
        last_time_relative = NA,
        max_duration       = NA
      ))
    }

    first_time    <- min(rv$doses$TIME, na.rm = TRUE)
    last_time     <- max(rv$doses$TIME, na.rm = TRUE)
    last_duration <- rv$doses %>% filter(TIME == last_time)  %>% pull(DUR)
    max_duration  <- max(rv$doses$DUR, na.rm = TRUE)

    dose_objects       <- rv$doses %>%
      mutate(TIME = round(as.numeric(difftime(TIME, first_time, units = "hours")), 1)) %>%
      arrange(TIME)
    last_time_relative <- max(dose_objects$TIME, na.rm = TRUE)

    list(
      first_dose_time    = first_time,
      last_dose_time     = last_time,
      last_dur           = last_duration,
      first_dur          = rv$doses %>% filter(TIME == first_time) %>% pull(DUR),
      max_time_dose      = last_time + as.difftime(last_duration, units = "hours"),
      dose_objects       = dose_objects,
      max_time_seq       = last_time_relative + 120 * 24,
      max_time_seq2      = last_time_relative + 30 * 24,
      last_time_relative = last_time_relative,
      max_duration       = max_duration
    )
  })


  # Patient data reactive -----------------------------------------------------

  patient_data <- reactive({
    age       <- input$age
    gender    <- input$gender
    height_cm <- input$height

    combined_data <- bind_rows(rv$weights, rv$creatinines, rv$albumines) %>%
      mutate(EVID = 1, DV = as.numeric(NA)) %>%
      arrange(TIME) %>%
      fill(WT, CREATININE, ALBUMIN, .direction = "downup") %>%
      distinct() %>%
      rowwise() %>%
      mutate(
        ABW   = round(calculate_abw(height_cm, WT, gender)),
        eGFR  = round(calculate_egfr(age, CREATININE, gender)),
        BSA   = round(calculate_bsa(WT, height_cm), 2),
        BSA_M = round(calculate_bsa_mosteller(WT, height_cm), 2),
        AGFR  = round((eGFR * BSA)   / 1.73),
        AGFR_M = round((eGFR * BSA_M) / 1.73),
        CLCR  = round(calculate_crcl(age, WT, CREATININE, gender)),
        BMI   = round(WT / (height_cm / 100)^2, 1)
      ) %>%
      ungroup()

    tdm_data <- rv$tdm %>% mutate(EVID = 0)
    if (!"DVID" %in% names(tdm_data)) tdm_data <- tdm_data %>% mutate(DVID = "Cc")

    combined_data <- bind_rows(combined_data, tdm_data)
    if (!"DVID" %in% names(combined_data)) combined_data <- combined_data %>% mutate(DVID = NA_character_)

    first_dose_time <- dose_info()$first_dose_time

    build_map_data <- function(time_seq) {
      bind_rows(combined_data, rv$doses) %>%
        filter(!is.na(TIME)) %>%
        mutate(TIME = round(as.numeric(difftime(TIME, first_dose_time, units = "hours")), 1)) %>%
        complete(TIME = time_seq) %>%
        arrange(TIME) %>%
        fill(AGFR, AGFR_M, ABW, eGFR, ALBUMIN, CLCR, WT, BMI, .direction = "downup") %>%
        mutate(
          AMT  = if_else(is.na(AMT),  0, AMT),
          DUR  = if_else(is.na(DUR),  0, DUR),
          EVID = if_else(is.na(EVID), 1, EVID),
          AGE  = age
        ) %>%
        filter(TIME >= 0)
    }

    list(
      absolute             = combined_data,
      patient_data_for_map  = build_map_data(dose_info()$max_time_seq),
      patient_data_for_map2 = build_map_data(dose_info()$max_time_seq2)
    )
  })


  # Table outputs -------------------------------------------------------------

  output$doseTable <- renderDT({
    rv$doses %>%
      mutate(TIME = format(TIME, "%d/%m/%Y %H:%M")) %>%
      select(-c(II, ADDL)) %>%
      datatable(options = list(dom = "t"))
  }, server = FALSE)

  output$TableWT <- renderDT({
    rv$weights %>%
      select(TIME, "Weight" = WT) %>%
      mutate(TIME = format(TIME, "%d/%m/%Y")) %>%
      datatable(options = list(dom = "t"))
  }, server = FALSE)

  output$TableCreatinine <- renderDT({
    rv$creatinines %>%
      select(TIME, "Creat" = CREATININE) %>%
      mutate(TIME = format(TIME, "%d/%m/%Y")) %>%
      datatable(options = list(dom = "t"))
  }, server = FALSE)

  output$TableAlbumin <- renderDT({
    rv$albumines %>%
      select(TIME, "Alb" = ALBUMIN) %>%
      mutate(TIME = format(TIME, "%d/%m/%Y")) %>%
      datatable(options = list(dom = "t"))
  }, server = FALSE)

  output$TableTDM <- renderDT({
    patient_data()$absolute %>%
      filter(EVID == 0) %>%
      select(TIME, "TDM" = DV) %>%
      mutate(TIME = format(TIME, "%d/%m/%Y %H:%M")) %>%
      datatable(options = list(dom = "t"))
  }, server = FALSE)

  output$combinedTable <- renderDT({
    patient_data()$absolute %>%
      filter(EVID == 1) %>%
      select(TIME, WT, ABW, CREATININE, AGFR) %>%
      mutate(TIME = format(TIME, "%d/%m/%Y")) %>%
      datatable(options = list(dom = "t"))
  }, server = FALSE)


  # Main calculation ----------------------------------------------------------

  observeEvent(input$calculate_button, {
    shinybusy::show_modal_spinner(text = "Running model averaging...")
    run_calculation()
    shinybusy::remove_modal_spinner()
  })

  run_calculation <- function() {
    set.seed(123)

    doses                  <- isolate(rv$doses)
    weights                <- isolate(rv$weights)
    creatinines            <- isolate(rv$creatinines)
    model_list             <- isolate(updated_model())
    tdm                    <- isolate(rv$tdm)
    d_info                 <- isolate(dose_info())
    patient_data_for_map   <- isolate(patient_data()$patient_data_for_map)
    patient_data_for_map2  <- isolate(patient_data()$patient_data_for_map2)
    free_drug              <- 1 - (input$target_protein_binding / 100)
    target_mic             <- input$target_mic
    target_auc             <- input$target_auc

    # Input validation --------------------------------------------------------
    if (nrow(doses) == 0)                { shinyalert("Error", "No doses data are provided",         type = "error"); return(NULL) }
    if (nrow(weights) == 0)              { shinyalert("Error", "No weight data are provided",        type = "error"); return(NULL) }
    if (nrow(creatinines) == 0)          { shinyalert("Error", "No creatinine data provided",        type = "error"); return(NULL) }
    if (length(model_list) == 0)         { shinyalert("Error", "No model selected for estimation.",  type = "error"); return(NULL) }
    if (nrow(tdm) == 0)                  { shinyalert("Error", "No TDM data",                        type = "error"); return(NULL) }
    if (min(tdm$TIME) < min(doses$TIME)) { shinyalert("Error", "TDM data is before the first dose",  type = "error"); return(NULL) }

    # MAP estimation (all models) ---------------------------------------------
    m_dt_out <- run_map_estimation(patient_data_for_map2, model_list, free_drug)

    prova <- m_dt_out %>%
      select(time, Cc = MoA_IPRED, AUC = MoA_AUC) %>%
      distinct()

    # Metrics -----------------------------------------------------------------
    data_metrics <- m_dt_out %>%
      left_join(patient_data_for_map %>% select(Cc_observed = DV, time = TIME), by = "time") %>%
      filter(!is.na(Cc_observed)) %>%
      split(f = .$model_name)

    metrics_long <- imap_dfr(data_metrics, calculate_metrics) %>%
      pivot_longer(cols = c(rBias, rRMSE, weight), names_to = "metric", values_to = "value") %>%
      dplyr::distinct()

    best_model <- metrics_long %>%
      dplyr::filter(metric == "weight") %>%
      dplyr::filter(value == max(value)) %>%
      pull(model_name) %>%
      unique()

    # AUC targets -------------------------------------------------------------
    auc_result           <- compute_auc_targets(prova, target_mic, target_auc, d_info$last_time_relative)
    first_below_target_ab <- auc_result$first_below_target_ab

    # SIR estimation (best model only) ----------------------------------------
    shinybusy::remove_modal_spinner()
    shinybusy::show_modal_spinner(text = "Estimating best model and calculating SIR...")

    quantiles_map <- run_sir_estimation(patient_data_for_map, model_list, best_model, d_info$max_time_seq)

    # Prepare plot data -------------------------------------------------------
    result     <- m_dt_out %>% dplyr::filter(model_name == best_model)

    data_plot <- result %>%
      mutate(Cocentration = round(Cc)) %>%
      filter(!(Cocentration <= 0.5 & time > d_info$max_duration)) %>%
      mutate(Time = d_info$first_dose_time + as.difftime(time, units = "hours"))

    data_plot_map <- quantiles_map %>%
      filter(time <= max(data_plot$time)) %>%
      mutate(Time = d_info$first_dose_time + as.difftime(time, units = "hours"))

    data_plot_dv <- patient_data_for_map %>%
      select(value = DV, time = TIME) %>%
      filter(time <= max(data_plot$time)) %>%
      mutate(Time = d_info$first_dose_time + as.difftime(time, units = "hours")) %>%
      left_join(data_plot %>% select(Time, Cocentration), by = "Time") %>%
      mutate(Error = round((value - Cocentration) / value, 2) * 100)

    data_prova <- m_dt_out %>%
      filter(!(Cc <= 0.5 & time > d_info$max_duration)) %>%
      mutate(
        Time       = d_info$first_dose_time + as.difftime(time, units = "hours"),
        MoA_IPRED  = round(MoA_IPRED, 2),
        Cc         = round(Cc, 2)
      )

    # Build plots -------------------------------------------------------------
    p_weights    <- build_weight_plot(metrics_long)
    combined_plot <- build_moa_plot(data_prova, data_plot_dv, p_weights)
    plot_dalba   <- build_single_model_plot(data_plot, data_plot_map, data_plot_dv, best_model)

    # Covariates table --------------------------------------------------------
    covariates_used <- unique(unlist(lapply(model_list, function(m) m$covariates)))
    patient_tb <- isolate(patient_data()$absolute) %>%
      filter(EVID == 1) %>%
      dplyr::select(any_of(c("TIME", covariates_used))) %>%
      mutate(TIME = format(TIME, "%d/%m/%Y"))

    # Render outputs ----------------------------------------------------------
    output$concentrationPlot  <- renderPlotly(combined_plot)
    output$concentrationPlot2 <- renderPlotly(plot_dalba)

    output$message <- renderUI({
      HTML(paste0(
        "Here you can find the results for the model averaging and single model prediction, ",
        "based on the most influential (highest weighted) model.<br><br>",
        "The redosing timing is calculated based on the averaged Bayesian predictions. ",
        "The single model plot shows also a 90% prediction interval derived through the Sequential Importance Resampling algorithm.<br><br>",
        "All estimations and simulations were performed using the <strong>posologyr</strong> package."
      ))
    })

    output$box_auc24 <- renderUI({
      date_to_display  <- d_info$first_dose_time + as.difftime(round(first_below_target_ab$time / 24), units = "days")
      days_from_last   <- round(first_below_target_ab$time / 24) - round(d_info$last_time_relative / 24)

      en_months <- c("January","February","March","April","May","June",
                     "July","August","September","October","November","December")
      date_formatted <- sprintf("%d %s %d",
        as.integer(format(date_to_display, "%d")),
        en_months[as.integer(format(date_to_display, "%m"))],
        as.integer(format(date_to_display, "%Y")))

      fluidRow(
        column(width = 6, offset = 3,
          box(
            title = tagList(
              icon("clock"),
              HTML("<strong style='font-size: 18px;'>Time to redose according to last 24h <span style='font-family: cursive;'>f</span>AUC<sub>0-24</sub>/MIC</strong>")
            ),
            status = "primary", solidHeader = TRUE, collapsible = FALSE, width = 12, height = "auto",
            div(
              h4(HTML(sprintf("%s (%d days from last dose)", date_formatted, days_from_last)),
                 style = "margin-bottom: 5px; font-weight: bold; color: #333;"),
              p(HTML(sprintf("When estimated to be: %.2f mg·h/L", first_below_target_ab$AUC24)),
                style = "font-size: 16px; color: #666;"),
              p(HTML(sprintf("And total C<sub>min</sub>: %.1f mg/L", first_below_target_ab$Cc)),
                style = "font-size: 16px; color: #666;"),
              style = "padding: 20px; text-align: center; background-color: #f9f9f9; border-radius: 5px;"
            )
          )
        )
      )
    })

    output$box_parameters <- renderUI({
      box(
        title = tagList(icon("notes-medical", lib = "font-awesome"), " Parameters"),
        status = "primary", solidHeader = TRUE, collapsible = FALSE, width = 12, height = "auto",
        div(style = "padding: 20px;", DTOutput("combinedTable"))
      )
    })

    output$combinedTable <- renderDT({
      datatable(patient_tb, options = list(dom = "t"))
    })

    # Report download ---------------------------------------------------------
    output$save_report <- downloadHandler(
      filename = function() {
        report_name <- ifelse(input$report_name != "", input$report_name, "Dalbavancin Report")
        paste0(report_name, " - ", Sys.Date(), ".html")
      },
      content = function(file) {
        rmarkdown::render(
          "report_template.Rmd",
          output_file = file,
          params = list(
            plot_dalba            = plot_dalba,
            combined_plot         = combined_plot,
            first_below_target_ab = first_below_target_ab,
            doses                 = doses,
            patient_tb            = patient_tb,
            date_to_display       = d_info$first_dose_time + as.difftime(round(first_below_target_ab$time / 24), units = "days"),
            days_from_last_dose   = round(first_below_target_ab$time / 24) - round(d_info$last_time_relative / 24),
            free_drug             = free_drug,
            target_mic            = target_mic,
            target_auc            = target_auc
          ),
          envir = new.env(parent = globalenv())
        )
      }
    )
  }
}

shinyApp(ui = ui, server = server)
