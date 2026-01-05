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
library(jsonlite)
library(shinyFeedback)
library(shinyBS)
library(waiter)
library(shinyalert)

options(scipen = 999) 
rxode2::setRxThreads(2L)

# Functions ---------------------------------------------------------------
# Calculate BSA
calculate_bsa <- function(weight_kg, height_cm) {
  bsa <- 0.007184 * (weight_kg^0.425) * (height_cm^0.725)
  return(bsa)
}

# Calculate ABW
calculate_abw <- function(height_cm, weight, gender) {
  height_in <- (height_cm / 100) * 39.37
  if (gender == "male") {
    IBW <- 50 + 2.3 * (height_in - 60)
  } else {
    IBW <- 45.5 + 2.3 * (height_in - 60)
  }
  ABW <- IBW + 0.4 * (weight - IBW)
  return(ABW)
}


# Calculate CrCl
calculate_crcl <- function(age, weight_kg, serum_creatinine_mg_dl, gender) {
  CrCl <- (140 - age) * weight_kg / (72 * serum_creatinine_mg_dl)
  if (gender == "female") {
    CrCl <- CrCl * 0.85
  }
  return(CrCl)
}

# Calculate eGFR
calculate_egfr <- function(age, serum_creatinine_mg_dl, gender) {
  if (gender == "male") {
    kappa <- 0.9
    alpha <- -0.302
    min_sc_k <- min(serum_creatinine_mg_dl / kappa, 1)
    max_sc_k <- max(serum_creatinine_mg_dl / kappa, 1)
    egfr <- 142 * (min_sc_k^alpha) * (max_sc_k^-1.200) * (0.9938^age)
  } else {
    kappa <- 0.7
    alpha <- -0.241
    min_sc_k <- min(serum_creatinine_mg_dl / kappa, 1)
    max_sc_k <- max(serum_creatinine_mg_dl / kappa, 1)
    egfr <- 142 * (min_sc_k^alpha) * (max_sc_k^-1.200) * (0.9938^age) * 1.012
  }
  return(egfr)
}

# Calculate metrics for model comparison

calculate_metrics <- function(data, model_name) {
  N <- nrow(data)
  predicted <- data$Cc
  observed <- data$Cc_observed
  weight <- data$weight
  
  rBias <- (1 / N) * sum((predicted - observed) / (observed)) * 100
  rRMSE <- sqrt((1 / N) * sum(((predicted - observed)^2) / ((observed)^2))) * 100
  
  tibble(model_name = model_name, rBias = rBias, rRMSE = rRMSE, weight = weight)
}

# Helper functions --------------------------------------------------------

# Reusable add/remove button group
dalbaAddRemoveButtons <- function(add_id, remove_id, weigth) {
  column(
    weigth,
    br(),
    actionButton(add_id, "+", class = "btn btn-success btn-xs"),
    actionButton(remove_id, "-", class = "btn btn-danger btn-xs")
  )
}

#Modele selection tab
dalbaModelSpecTab <- function(model_names) {
  tabItem(
    tabName = "model_specification",
    box(
      width = 12, solidHeader = TRUE, status = "primary", collapsible = FALSE,
      fluidRow(
        column(4,
               tags$label("Models included in the automated model selection/averaging"),
               checkboxGroupButtons(
                 inputId = "model_selection",
                 label = NULL,
                 selected = model_names,
                 choices = model_names,
                 size = "lg",
                 direction = "vertical",
                 justified = TRUE,
                 checkIcon = list(
                   yes = icon("ok", 
                              lib = "glyphicon"))
               )
        ),
        column(8,
               tags$p("Select one or more models to include in the automated model selection or model averaging algorithm. These models will be used downstream in the calculations or visualizations.")
        )
      )
    )
  )
}

# Generic add/remove observer function
dalbaRegisterAddRemove <- function(input, output, session, input_add, input_remove, input_time, input_value, var_name, var_col, rv, time_format = "%d/%m/%Y %H:%M") {
  
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


# Models ------------------------------------------------------------------
Dalbavancin <- list(
  Gregoire2025 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0
      TVVc <- THETA_Vc * (ABW / 66.5) ^ 0.6
      TVCl <- THETA_Cl * (AGFR / 88.5) ^ 0.21
      Vc <- TVVc * exp(ETA_Vc)
      Cl <- TVCl * exp(ETA_Cl)
      Vp <- THETA_Vp * exp(ETA_Vp)
      Q <- THETA_Q * exp(ETA_Q)
      ke <- Cl / Vc
      k12 <- Q / Vc
      k21 <- Q / Vp
      Cc <- centr / Vc
      Cp <- periph / Vp
      d / dt(centr) <- -ke * centr - k12 * centr + k21 * periph
      d / dt(periph) <- +k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_Cl = 0.04,
      THETA_Vc = 5.62,
      THETA_Vp = 8.04,
      THETA_Q = 0.026
    ),
    covariates = c("ABW", "AGFR"),
    omega = lotri::lotri({
      ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q ~
        c(4.096472, 0, 3.61458, 0, 0, 0, 0, 0, 0, 0)
    }),
    sigma = c(additive_a = 1.67, proportional_b = 0.22)
  ),
  Cojutti2022 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0
      TVCl <- THETA_Cl * (eGFR / 93) ^ 0.0043
      V1   <- THETA_V1 * exp(ETA_V1)
      V2   <- THETA_V2 * exp(ETA_V2)
      Q    <- THETA_Q  * exp(ETA_Q)
      Cl   <- TVCl     * exp(ETA_Cl)
      ke   <- Cl / V1
      k12  <- Q / V1
      k21  <- Q / V2
      Cc   <- centr / V1
      Cp   <- periph / V2
      d / dt(centr) <- -ke * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    
    theta = c(
      THETA_Cl = 0.029,
      THETA_V1 = 6.14,
      THETA_V2 = 9.52,
      THETA_Q  = 0.026
    ),
    covariates = c("eGFR"),
    omega = lotri::lotri({
      ETA_Cl + ETA_V1 + ETA_V2 + ETA_Q ~
        c(0.0654, 0, 0.0255, 0, 0, 0.1213, 0, 0, 0, 0.2063)
    }),
    sigma = c(additive_a = 0, proportional_b = 0.3392)
  ),
  Carrothers2020 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph1(0) <- 0
      periph2(0) <- 0
      
      Cl <- THETA_Cl * (ALBUMIN / 3.7) ^ THETA_CL_ALB *
        (CLCR / 100) ^ THETA_CL_CLCR *
        (WT / 85.5) ^ THETA_CL_WT * exp(ETA_Cl)
      
      V1 <- THETA_V1 * (ALBUMIN / 3.7) ^ THETA_V1_ALB *
        (WT / 85.5) ^ THETA_V1_WT * exp(ETA_V1)
      
      V2 <- THETA_V2 * (AGE / 47) ^ THETA_V2_AGE *
        (ALBUMIN / 3.7) ^ THETA_V2_ALB *
        (WT / 85.5) ^ THETA_V2_WT * exp(ETA_V2)
      
      V3 <- THETA_V3 * (ALBUMIN / 3.7) ^ THETA_V3_ALB *
        (WT / 85.5) ^ THETA_V3_WT * exp(ETA_V3)
      
      Q2 <- THETA_Q2
      Q3 <- THETA_Q3
      
      k10 <- Cl / V1
      k12 <- Q2 / V1
      k21 <- Q2 / V2
      k13 <- Q3 / V1
      k31 <- Q3 / V3
      
      Cc <- centr / V1
      d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph1 - k13 * centr + k31 * periph2
      d / dt(periph1) <- +k12 * centr - k21 * periph1
      d / dt(periph2) <- +k13 * centr - k31 * periph2
      d / dt(AUC)     <- Cc
    }),
    
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    # Fixed effects (Table 3)
    theta = c(
      THETA_Cl      = 0.0531,
      THETA_V1      = 3.04,
      THETA_V2      = 8.78,
      THETA_V3      = 3.28,
      THETA_Q2      = 0.288,
      THETA_Q3      = 2.11,
      THETA_CL_ALB  = -0.477,
      THETA_CL_CLCR = 0.273,
      THETA_CL_WT   = 0.391,
      THETA_V1_ALB  = -0.340,
      THETA_V1_WT   = 0.683,
      THETA_V2_AGE  = 0.486,
      THETA_V2_ALB  = -0.413,
      THETA_V2_WT   = 0.365,
      THETA_V3_ALB  = -0.551,
      THETA_V3_WT   = 0.518
    ),
    
    covariates = c("ALBUMIN", "CLCR", "WT", "AGE"),
    
    # Omega matrix (IIV) — diagonal only; correlations omitted for simplicity
    omega = lotri::lotri({
      ETA_Cl + ETA_V2 + ETA_V1 + ETA_V3 ~
        c(log((22 / 100) ^ 2 + 1), 0, log((41 / 100) ^ 2 + 1), 0, 0, log((24 / 100) ^ 2 + 1), 0, 0, 0, log((74 / 100) ^ 2 + 1))
    }),
    sigma = c(additive_a = 0, proportional_b = 0.0362)
  ),
  Cojutti2021 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0
      
      CL = THETA_CL * (eGFR / 90)^THETA_CL_CLCR * exp(ETA_CL)
      V1 = THETA_V1 * exp(ETA_V1)
      V2 = THETA_V2 * exp(ETA_V2)
      Q  = THETA_Q  * exp(ETA_Q)
      
      # Rate constants
      k10 = CL / V1
      k12 = Q / V1
      k21 = Q / V2
      
      # Concentrations
      Cc = centr / V1
      
      # Differential equations
      d/dt(centr)  = -k10 * centr - k12 * centr + k21 * periph
      d/dt(periph) =  k12 * centr - k21 * periph
      d/dt(AUC)    = Cc
    }),
    
    error_model = function(f, sigma) {
      g = sigma[1] + sigma[2] * f
      return(g)
    },
    
    theta = c(
      THETA_CL = 0.106,      # L/h
      THETA_V1 = 4.53,       # L
      THETA_V2 = 31.9,       # L
      THETA_Q  = 0.51,       # L/h
      THETA_CL_CLCR = 0.5    # Exponent for CLCR effect
    ),
    
    covariates = c("eGFR"),
    
    omega = lotri::lotri({
      ETA_CL + ETA_V1 + ETA_V2 + ETA_Q ~
        c(
          log((20 / 100)^2 + 1),  
          0, log((25 / 100)^2 + 1),  
          0, 0, log((30 / 100)^2 + 1),  
          0, 0, 0, log((35 / 100)^2 + 1)  
        )
    }),
    
    sigma = c(additive_a = 0, proportional_b = 0.15)  # 15% proportional error
  ),
  Baiardi2025 = list(
    ppk_model = rxode2::rxode2({
      centr(0) = 0
      periph(0) = 0
      AUC(0) = 0
      
      CL = THETA_CL * (WT / 70)^0.75 * exp(ETA_CL)
      V1 = THETA_V1 * (WT / 70)^1.0 * exp(ETA_V1)
      V2 = THETA_V2 * (WT / 70)^1.0 * exp(ETA_V2)
      Q  = THETA_Q  * (WT / 70)^0.75 * exp(ETA_Q)
      
      k10 = CL / V1
      k12 = Q / V1
      k21 = Q / V2
  
      Cc = centr / V1
      
      d/dt(centr)  = -k10 * centr - k12 * centr + k21 * periph
      d/dt(periph) =  k12 * centr - k21 * periph
      d/dt(AUC)    = Cc
    }),
    
    error_model = function(f, sigma) {
      g = sigma[1] + sigma[2] * f
      return(g)
    },
    
    theta = c(
      THETA_CL = 0.053,  # L/h
      THETA_V1 = 3.04,   # L
      THETA_V2 = 8.78,   # L
      THETA_Q  = 0.288   # L/h
    ),
    
    covariates = c("WT"),
    
    omega = lotri::lotri({
      ETA_CL + ETA_V1 + ETA_V2 + ETA_Q ~
        c(
          log((22 / 100)^2 + 1),  
          0, log((24 / 100)^2 + 1),  
          0, 0, log((41 / 100)^2 + 1),  
          0, 0, 0, log((30 / 100)^2 + 1)  
        )
    }),
    
    sigma = c(additive_a = 0, proportional_b = 0.0362)  # 3.62% proportional error
  )
)




# UI ----------------------------------------------------------------------
ui <- dashboardPage(
  dashboardHeader(title = tagList(
    span(class = "logo-lg", "Dalbavancin"),
    icon("pills", class = "fa-lg")
  )),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Patient Info", tabName = "patient_info", icon = icon("user")),
      menuItem(
        "Dose Observations",
        tabName = "dose_observations",
        icon = icon("pills")
      ),
      menuItem(
        "Model Targets",
        tabName = "model_targets",
        icon = icon("bullseye")
      ),
      menuItem(
        "Model Specification",
        tabName = "model_specification",
        icon = icon("cogs")
      ),
      menuItem(
        "Report/Load/Save",
        tabName = "report",
        icon = icon("file-pdf")
      ),
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
      # Patient Info Tab
      tabItem(tabName = "patient_info", fluidRow(
        column(
          width = 8,
          box(
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            fluidRow(
              column(4, numericInput("age", "Age (years)", value = 70)),
              column(
                4,
                selectInput(
                  "gender",
                  "Gender",
                  choices = c("male", "female"),
                  selected = "male"
                )
              ),
              column(4, numericInput("height", "Height (cm)", value = 170))
            ),
            h4("Weight"),
            fluidRow(
              column(
                4,
                airDatepickerInput(
                  "weight_time",
                  "Day of record",
                  value = Sys.Date(),
                  width = "100%"
                )
              ),
              column(4, numericInput("weight", "Weight (kg)", value = 70)),
              dalbaAddRemoveButtons("addWeight", "removeWeight", 4)
            ),
            h4("Creatinine"),
            fluidRow(
              column(
                4,
                airDatepickerInput(
                  "creatinine_time",
                  "Day of record",
                  value = Sys.Date(),
                  width = "100%"
                )
              ),
              column(4, numericInput("creatinine", "Serum Cr (mg/dL)", value = 0.9)),
              dalbaAddRemoveButtons("addCreatinine", "removeCreatinine", 4)
            ),
            h4("Albumin"),
            fluidRow(
              column(
                4,
                airDatepickerInput(
                  "albumin_time",
                  "Day of record",
                  value = Sys.Date(),
                  width = "100%"
                )
              ),
              column(4, numericInput("albumin", "Albumin (g/dL)", value = 4)),
              dalbaAddRemoveButtons("addAlbumin", "removeAlbumin", 4)
            ),
            h4("TDM"),
            fluidRow(
              column(
                4,
                airDatepickerInput(
                  "tdm_time",
                  "Day/time of record",
                  value = Sys.time(),
                  timepicker = TRUE,
                  timepickerOpts = timepickerOptions(hoursStep = 1, minutesStep = 10),
                  width = "100%"
                )
              ),
              column(4, numericInput("tdm", "TDM (mg/L)", value = 0)),
              dalbaAddRemoveButtons("addTDM", "removeTDM", 4)
            )
          )
        ), column(
          width = 4,
          box(
            title = "Patient Data",
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            collapsible = TRUE,
            DTOutput("TableWT"),
            hr(),
            DTOutput("TableCreatinine"),
            hr(),
            DTOutput("TableAlbumin"),
            hr(),
            DTOutput("TableTDM")
          )
        )
      )),
      
      # Dose Observations Tab
      tabItem(tabName = "dose_observations", fluidRow(
        column(
          width = 6,
          box(
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            airDatepickerInput(
              "dose_time",
              "Day of Administration",
              value = Sys.time(),
              timepicker = TRUE,
              timepickerOpts = timepickerOptions(hoursStep = 1, minutesStep = 5),
              width = "100%"
            ),
            pickerInput(
              "dose",
              "Dose (mg)",
              choices = c(500, 1000, 1500),
              selected = 1500
            ),
            numericInput("dur", "Infusion rate (hours)", value = 0.5),
            fluidRow(dalbaAddRemoveButtons("addDose", "removeDose", 6))
          )
        ), column(
          width = 6,
          box(
            width = 12,
            solidHeader = TRUE,
            status = "primary",
            DTOutput("doseTable", width = "75%")
          )
        )
      )),
      
      # Model Targets Tab
      tabItem(tabName = "model_targets", fluidRow(column(
        width = 12,
        box(
          width = 12,
          solidHeader = TRUE,
          status = "primary",
          fluidRow(
            column(
              3,
              numericInput("target_auc", "Targeted last 24h fAUC0-24h/MIC", value = 111),
              numericInput("target_mic", "Target MIC", value = 0.125),
              numericInput("target_protein_binding", "Protein binding (%)", value = 93),
              actionButton("reset_target", "Reset", class = "btn-primary", style = "width: 100%;color:white;")
            ),
            column(9, tags$div(
              tags$h4("Explanation of Targets"),
              tags$p(
                "The 'Targeted last 24h fAUC0-24h/MIC' represents the free drug exposure over the last 24 hours relative to the MIC."
              ),
              tags$p(
                "The 'Targeted cumulative fAUC0-24h/MIC' measures the cumulative free drug exposure over a dosing interval relative to the MIC."
              ),
              tags$p(
                "According to EUCAST, S. aureus, E. faecalis and E. faecium have MICs/ECOFF of 0.125 mg/L for dalbavancin."
              )
            ))
          )
        )
      ))),
      
      # Model Specification Tab
      tabItem(tabName = "model_specification", dalbaModelSpecTab(model_names = names(Dalbavancin))),
      
      # Report Tab
      tabItem(
        tabName = "report",
        box(
          title = "Report Management",
          width = 12,
          solidHeader = TRUE,
          status = "primary",
          textInput("report_name", "Report Name", value = "Dalbavancin Report"),
          fluidRow(
            column(
              4,
              downloadButton(
                "save_report",
                "Save Report",
                icon = icon("save"),
                class = "btn-primary",
                style = "width: 100%;color:white;"
              )
            ),
            column(
              4,
              downloadButton(
                "save_state",
                "Save State",
                icon = icon("save"),
                class = "btn-primary",
                style = "width: 100%; color:white;"
              )
            ),
            column(
              4,
              fileInput(
                "load_state",
                label = NULL,
                buttonLabel = "Load State",
                accept = c(".zip"),
                width = "100%",
                placeholder = "Choose .rds file"
              )
            )
          )
        )
      )
    ),
    fluidRow(column(
      width = 12,
      box(
        width = 12,
        solidHeader = TRUE,
        status = "primary",
        collapsible = FALSE,
        title = "",
        actionButton("calculate_button", "Run calculation", class = "btn btn-primary btn-lg"),
        br(),
        br(),
        uiOutput("message"),
        br(),
        plotlyOutput("concentrationPlot", height = "400px"),
        br(),
        plotlyOutput("concentrationPlot2", height = "400px"),
        br(),
        uiOutput("box_auc24"),
        br(),
        uiOutput("box_parameters")
      )
    )),
    scrollToTop = TRUE,
    preloader = list(html = spin_1(), color = "gray")
  )
)


# Server ------------------------------------------------------------------

server <- function(input, output, session) {
  # Initial modal fo user agreement -----------------------------------------
  session$onSessionEnded(function() { 
    stopApp()
  })
  
  showModal(modalDialog(
    title = "Dalbavancin Precision Dosing Tool",
    p("Disclaimer: This tool is for informational purposes only and intended for use only by health care professionals. The tool is not intended to be a substitute for professional medical advice, dosing, diagnosis or treatment"),
    footer = tagList(
      fluidRow(
        column(
          3,
          checkboxInput("accept", "I understand", value = FALSE)
        ),
        column(
          9,
          actionButton("confirm", "Confirm", class = "btn-primary")
        )
      )
    ),
    easyClose = FALSE, # Prevent the modal from closing without confirmation
    size = "m" # Small size
  ))
  
  shinyjs::disable("confirm") # usign shinyjs the button is "grey" and not clickable until the user accept the disclaimer
  
  observeEvent(input$accept, {
    if (input$accept) {
      shinyjs::enable("confirm")
    } else {
      shinyjs::disable("confirm")
    }
  })
  
  observeEvent(input$confirm, {
    removeModal()
  })
  
  # Reactive values used ----------------------------------------------------
  
  rv <- reactiveValues(
    weights = tibble(TIME = as.POSIXct(character()), WT = numeric()),
    creatinines = tibble(TIME = as.POSIXct(character()), CREATININE = numeric()),
    albumines = tibble(TIME = as.POSIXct(character()), ALBUMIN = numeric()),
    doses = tibble(
      TIME = as.POSIXct(character()),
      AMT = numeric(),
      DUR = numeric(),
      ADDL = numeric(),
      II = numeric()
    ),
    tdm = tibble(TIME = as.POSIXct(character()), DV = numeric())
  )
  
  
  # Save/Load data ----------------------------------------------------------
  
  save_state <- function(file) {
    state <- list(
      patient_data = patient_data()$absolute,
      doses = rv$doses,
      weights = rv$weights,
      creatinines = rv$creatinines,
      albumines = rv$albumines,
      tdm = rv$tdm,
      target_auc = input$target_auc,
      target_cum_auc = input$target_cum_auc,
      target_mic = input$target_mic,
      target_protein_binding = input$target_protein_binding
    )
    
    temp_file <- tempfile(fileext = ".rds")
    saveRDS(state, temp_file)
    zip(file, temp_file, flags = "-j")
    showNotification("State saved successfully", type = "message")
  }
  
  load_state <- function(file) {
    # Create a temporary directory specifically for this operation
    temp_dir <- file.path(tempdir(), "unzip_temp")
    dir.create(temp_dir, showWarnings = FALSE)
    
    # Unzip the file into the temporary directory
    unzip(file, exdir = temp_dir)
    
    # Look for .rds file in the extracted files
    rds_file <- list.files(temp_dir, pattern = "\\.rds$", full.names = TRUE)
    
    if (length(rds_file) == 1) {
      # Read the state from the .rds file
      state <- readRDS(rds_file)
      # Restore the state into your reactive values and inputs
      rv$weights <- state$weights
      rv$creatinines <- state$creatinines
      rv$doses <- state$doses
      rv$tdm <- state$tdm
      updateNumericInput(session, "target_auc", value = state$target_auc)
      updateNumericInput(session, "target_cum_auc", value = state$target_cum_auc)
      updateNumericInput(session, "target_mic", value = state$target_mic)
      updateNumericInput(session, "target_protein_binding", value = state$target_protein_binding)
      
      showNotification("State loaded successfully", type = "message")
    } else {
      showNotification("Error: Could not find the .rds file in the zip archive", type = "error")
    }
    
    # Clean up the temporary directory
    unlink(temp_dir, recursive = TRUE)
  }
  
  output$save_state <- downloadHandler(
    filename = function() {
      paste0("state_", Sys.Date(), ".zip")
    },
    content = function(file) {
      save_state(file)
    }
  )
  
  observeEvent(input$load_state, {
    req(input$load_state)
    load_state(input$load_state$datapath)
  })
  
  
  # Update Model ------------------------------------------------------------
  
  updated_model <- reactive({
    selected_models <- input$model_selection
    Dalbavancin_models <- Dalbavancin[names(Dalbavancin) %in% selected_models]
    return(Dalbavancin_models)
  })
  
  
  ## Resets  -----------------------------------------------------------------
  
  observeEvent(input$reset_target, {
    updateNumericInput(session, "target_auc", value = 111)
    updateNumericInput(session, "target_mic", value = 0.125)
    updateNumericInput(session, "target_protein_binding", value = 93)
  })
  
  ## Validation of inputs ----------------------------------------------------
  validate_input <- function(input_id, message) {
    observe({
      hideFeedback(input_id)
      if (is.null(input[[input_id]]) || is.na(input[[input_id]]) || input[[input_id]] == "") {
        showFeedbackDanger(input_id, message)
      }
    })
  }
  
  validate_input("target_auc", "Target free-AUC is required")
  validate_input("target_mic", "Target MIC is required")
  validate_input("target_protein_binding", "Protein binding is required")
  validate_input("age", "Age is required")
  validate_input("height", "Height is required")
  
  # Covariates values -------------------------------------------------------
  
  input_handlers <- list(
    list(add = "addWeight", remove = "removeWeight", time = "weight_time", value = "weight", var = "weights", col = "WT"),
    list(add = "addCreatinine", remove = "removeCreatinine", time = "creatinine_time", value = "creatinine", var = "creatinines", col = "CREATININE"),
    list(add = "addAlbumin", remove = "removeAlbumin", time = "albumin_time", value = "albumin", var = "albumines", col = "ALBUMIN"),
    list(add = "addTDM", remove = "removeTDM", time = "tdm_time", value = "tdm", var = "tdm", col = "DV")
  )
  
  # Register all handlers in one loop
  purrr::walk(input_handlers, function(cfg) {
    dalbaRegisterAddRemove(
      input, output, session,
      input_add = cfg$add,
      input_remove = cfg$remove,
      input_time = cfg$time,
      input_value = cfg$value,
      var_name = cfg$var,
      var_col = cfg$col,
      rv = rv
    )
  })
  
  # Add dose event ----------------------------------------------------------
  # It is simpler in contrast to vancomycin beacuse of the way doses are added, just one at a time
  
  output$doseTable <- renderDT(
    {
      doses_tb <- rv$doses %>%
        mutate(TIME = format(TIME, "%d/%m/%Y %H:%M")) %>%
        select(-c(II, ADDL))
      datatable(doses_tb, options = list(dom = "t"))
    },
    server = FALSE
  )
  
  observeEvent(input$addDose, {
    newDose <- tibble(
      TIME = as.POSIXct(input$dose_time, format = "%d/%m/%Y %H:%M"),
      AMT = as.numeric(input$dose),
      DUR = input$dur,
      ADDL = 0,
      II = 0
    ) %>%
      distinct(TIME, .keep_all = TRUE) %>%
      arrange(TIME)
    
    if (!any(rv$doses$TIME == newDose$TIME)) {
      rv$doses <- bind_rows(rv$doses, newDose)
    }
  })
  
  observeEvent(input$removeDose, {
    if (nrow(rv$doses) > 0) {
      rv$doses <- rv$doses[-nrow(rv$doses), ]
    }
  })
  
  dose_info <- reactive({
    if (nrow(rv$doses) == 0) {
      return(list(
        dose_expanded = tibble(),
        first_dose_time = NA,
        last_dose_time = NA,
        last_dur = NA,
        dose_objects = tibble(),
        max_time_seq = NA,
        last_time_relative = NA,
        max_duration = NA
      ))
    }
    
    first_time <- min(rv$doses$TIME, na.rm = TRUE) # absolute
    last_time <- max(rv$doses$TIME, na.rm = TRUE) # absolute
    last_duration <- rv$doses %>%
      filter(TIME == last_time) %>%
      select(DUR) %>%
      pull()
    first_duration <- rv$doses %>%
      filter(TIME == first_time) %>%
      select(DUR) %>%
      pull()
    max_duration <- max(rv$doses$DUR, na.rm = TRUE)
    
    # Part for relative values
    dose_objects <- rv$doses %>%
      mutate(TIME = round(as.numeric(difftime(TIME, first_time, units = "hours")), 1)) %>%
      arrange(TIME)
    max_time_seq <- max(dose_objects$TIME, na.rm = TRUE) + 120 * 24
    max_time_seq2 <-  max(dose_objects$TIME, na.rm = TRUE) + 30 * 24
    last_time_relative <- max(dose_objects$TIME, na.rm = TRUE)
    
    
    list(
      first_dose_time = first_time,
      last_dose_time = last_time,
      last_dur = last_duration,
      first_dur = first_duration,
      max_time_dose = last_time + as.difftime(last_duration, units = "hours"),
      dose_objects = dose_objects,
      max_time_seq = max_time_seq,
      max_time_seq2 = max_time_seq2,
      last_time_relative = last_time_relative,
      max_duration = max_duration
    )
  })
  
  
  # Table outputs -----------------------------------------------------------
  output$combinedTable <- renderDT(
    {
      patient_tb <- patient_data()$absolute %>%
        filter(EVID == 1) %>%
        select(TIME, WT, ABW, CREATININE, AGFR) %>%
        mutate(TIME = format(TIME, "%d/%m/%Y"))
      datatable(patient_tb, options = list(dom = "t"))
    },
    server = FALSE
  )
  
  output$TableTDM <- renderDT(
    {
      patient_tb <- patient_data()$absolute %>%
        filter(EVID == 0) %>%
        select(TIME, "TDM" = DV) %>%
        mutate(TIME = format(TIME, "%d/%m/%Y %H:%M"))
      datatable(patient_tb, options = list(dom = "t"))
    },
    server = FALSE
  )
  
  output$TableWT <- renderDT(
    {
      patient_tb <- rv$weights %>%
        select(TIME, "Weight" = WT) %>%
        mutate(TIME = format(TIME, "%d/%m/%Y"))
      datatable(patient_tb, options = list(dom = "t"))
    },
    server = FALSE
  )
  
  output$TableCreatinine <- renderDT(
    {
      patient_tb <- rv$creatinines %>%
        select(TIME, "Creat" = CREATININE) %>%
        mutate(TIME = format(TIME, "%d/%m/%Y"))
      datatable(patient_tb, options = list(dom = "t"))
    },
    server = FALSE
  )  
  
  output$TableAlbumin <- renderDT(
    {
      patient_tb <- rv$albumines %>%
        select(TIME, "Alb" = ALBUMIN) %>%
        mutate(TIME = format(TIME, "%d/%m/%Y"))
      datatable(patient_tb, options = list(dom = "t"))
    },
    server = FALSE
  )  
  
  # Reactive datasets for covariates ----------------------------------------
  patient_data <- reactive({
    age <- input$age
    gender <- input$gender
    height_cm <- input$height
    
    combined_data <- bind_rows(rv$weights, rv$creatinines, rv$albumines) %>% mutate(EVID = 1, DV = as.numeric(NA)) 
    
    combined_data <- combined_data %>%
      arrange(TIME) %>%
      fill(WT, CREATININE, ALBUMIN, .direction = "downup") %>%
      distinct() %>%
      rowwise() %>% # this is important to correctly perform the function in each row, otherwise it will perform on the last imputed value and fill all other rows
      mutate(
        ABW = round(calculate_abw(height_cm, WT, gender)),
        eGFR = round(calculate_egfr(age, CREATININE, gender)),
        BSA = round(calculate_bsa(WT, height_cm), 2), # important! round to 2 decimal to avoid problems with the model, all BSA are the same
        AGFR = round((eGFR * BSA) / 1.73),
        CLCR = round(calculate_crcl(age, WT, CREATININE, gender))
      ) %>%
      ungroup() 
    
    combined_data <- bind_rows(
      combined_data,
      rv$tdm %>% mutate(EVID = 0))
    
    first_dose_time <- dose_info()$first_dose_time
    time_seq <- dose_info()$max_time_seq
    time_seq2 <- dose_info()$max_time_seq2
    
    patient_data_for_map <- bind_rows(combined_data, rv$doses) %>% 
      filter(!is.na(TIME)) %>%
      mutate(
        TIME = round(as.numeric(difftime(TIME, first_dose_time, units = "hours")), 1)
      ) %>%
      complete(TIME = time_seq) %>%
      arrange(TIME) %>%
      fill(AGFR, ABW, eGFR, ALBUMIN, CLCR, WT, .direction = "downup") %>%
      mutate(
        AMT = if_else(is.na(AMT), 0, AMT),
        DUR = if_else(is.na(DUR), 0, DUR),
        EVID = if_else(is.na(EVID), 1, EVID),
        AGE = age
      ) %>%
      filter(TIME >= 0)
    
    
    patient_data_for_map2 <- bind_rows(combined_data, rv$doses) %>% 
      filter(!is.na(TIME)) %>%
      mutate(
        TIME = round(as.numeric(difftime(TIME, first_dose_time, units = "hours")), 1)
      ) %>%
      complete(TIME = time_seq2) %>%
      arrange(TIME) %>%
      fill(AGFR, ABW, eGFR, ALBUMIN, CLCR, WT, .direction = "downup") %>%
      mutate(
        AMT = if_else(is.na(AMT), 0, AMT),
        DUR = if_else(is.na(DUR), 0, DUR),
        EVID = if_else(is.na(EVID), 1, EVID),
        AGE = age
      ) %>%
      filter(TIME >= 0)
    
    return(list(
      absolute = combined_data, # used for table outputs
      patient_data_for_map = patient_data_for_map, # used for MAP estimation of the best model
      patient_data_for_map2 = patient_data_for_map2 #used for MAP estimation for all models, to spare time 
    ))
  })
  
  # Core processing ---------------------------------------------------------
  
  observeEvent(input$calculate_button, {
    shinybusy::show_modal_spinner(text = "Running model averaging...")
    core_function()
    shinybusy::remove_modal_spinner()
  })
  
  core_function <- function() {
    set.seed(123)
    doses <- isolate(rv$doses)
    weights <- isolate(rv$weights)
    creatinines <- isolate(rv$creatinines)
    albumines <- isolate(rv$albumines)
    model_list <- isolate(updated_model())
    model_list_names <- as.list(names(model_list))
    tdm <- isolate(rv$tdm)
    dose_info <- isolate(dose_info())
    patient_data_for_map <- isolate(patient_data()$patient_data_for_map)
    patient_data_for_map2 <- isolate(patient_data()$patient_data_for_map2)
    patient_data_for_table <- isolate(patient_data()$absolute)
    time_seq <- dose_info$max_time_seq
    from <- dose_info$last_time_relative
    free_drug <- 1 - (input$target_protein_binding / 100)
    target_mic <- input$target_mic
    target_auc <- input$target_auc
    
    
    ## Error messages ----------------------------------------------------------
    
    if (nrow(doses) == 0) {
      shinyalert(
        title = "Error",
        text = "No doses data are provided",
        type = "error"
      )
      return(NULL)
    }
    
    if (nrow(weights) == 0) {
      shinyalert(
        title = "Error",
        text = "No weight data are provided",
        type = "error"
      )
      return(NULL)
    }
    
    if (nrow(creatinines) == 0) {
      shinyalert(
        title = "Error",
        text = "No creatinine data provided",
        type = "error"
      )
      return(NULL)
    }
    
    if (length(model_list) == 0) {
      shinyalert(
        title = "Error",
        text = "No model selected for estimation.",
        type = "error"
      )
      return(NULL)
    }
    
    if (min(tdm$TIME) < min(doses$TIME)) {
      shinyalert(
        title = "Error",
        text = "TDM data is before the first dose",
        type = "error"
      )
      return(NULL)
    }
    
    m_out <- list()
    error_log <- list()
    
    for (i in seq_along(model_list)) {
      result <- tryCatch({
        m_i <- poso_estim_map(patient_data_for_map2, model_list[[i]], return_model = TRUE, return_ofv = TRUE)
        
        # Store successful result
        m_out[[i]] <- cbind(
          data.table::data.table(m_i$model)[, .(time, Cc, AUC)],
          ofv = m_i$ofv,
          LL = exp(-0.5 * m_i$ofv),
          model_name = model_list_names[[i]],
          model_number = i
        )
        NULL  # No error
      }, error = function(e) {
        # Store error message with model name
        error_log[[length(error_log) + 1]] <<- list(
          model_name = model_list_names[[i]],
          model_number = i,
          error_message = e$message
        )
        NULL  # Prevent crash
      })
    }
    
    m_dt_out <- data.table::as.data.table(do.call(rbind, m_out))
    m_dt_out[, weight := LL / sum(unique(LL))]
    m_dt_out[, MoA_IPRED := sum(Cc * weight), by = time]
    m_dt_out[, MoA_AUC := sum(AUC * weight), by = time]
    
    
    # Get the model avaraged values of Cc and AUC to use later to calculate the
    # best time for resampling
    prova <- m_dt_out %>% 
      select(
        time, Cc = MoA_IPRED, AUC = MoA_AUC
      ) %>% 
      distinct()
    
    
  # Capire qui cosa ha senso fare per le metriche, tenere solo quella di avarege?
    # se si allora prendi prova e calcola su qullo, non sullo split
    data_metrics <- m_dt_out %>%
      left_join(patient_data_for_map %>%
                  select(Cc_observed = DV, time = TIME), by = "time") %>%
      filter(!is.na(Cc_observed)) %>%
      split(f = .$model_name)
    
    metrics_summary <- imap_dfr(data_metrics, calculate_metrics)
    
    metrics_long <- metrics_summary %>%
      pivot_longer(cols = c(rBias, rRMSE, weight), names_to = "metric", values_to = "value") %>%
      dplyr::distinct()
    
    best_model <- metrics_long %>%  dplyr::filter(metric == "weight") %>% dplyr::filter(value == max(value)) %>%
      pull(model_name) %>% unique()
    
    # Weight plot, for MoA
    p_weights <- plot_ly(
      data = metrics_long %>%
        dplyr::filter(metric == "weight"),
      x = ~"Weight",
      y = ~value,
      color = ~model_name,
      colors = "Set2",
      type = "bar",
      hoverinfo = "text",
      showlegend = FALSE,
      text = ~ paste(model_name, ":", round(value, 2))
    ) %>%
      layout(
        barmode = "stack",
        yaxis = list(title = "", range = c(0, 1.01), side = "right"),
        xaxis = list(title = ""),
        showlegend = FALSE
      )
    
    # MAP Estimations
    shinybusy::remove_modal_spinner()
    shinybusy::show_modal_spinner(text = "Estimating best model and calculating SIR...")
    
    result <- m_dt_out %>% dplyr::filter(model_name == best_model)
    sir_map <- poso_estim_sir(patient_data_for_map, model_list[[best_model]], n_sample = 1e3, n_resample = 1e2)
    event_table <- rxode2::as.et(patient_data_for_map)
    event_table$add_sampling(seq(0, time_seq, by = 1))
    
    ci_map <-  poso_replace_et(
      target_model = sir_map$model,
      prior_model = model_list[[best_model]],
      event_table = event_table,
      interpolation = "locf"
    )
    
    quantiles_map <- ci_map %>%
      group_by(time) %>%
      summarise(
        lower90 = round(quantile(Cc, probs = 0.05, na.rm = TRUE)),
        upper90 = round(quantile(Cc, probs = 0.95, na.rm = TRUE))
      ) %>%
      ungroup()
    
    data_plot <- result %>%
      mutate(Cocentration = round(Cc)) %>%
      filter(!(Cocentration <= 0.5 & time > dose_info$max_duration)) %>%
      mutate(Time = dose_info$first_dose_time + as.difftime(time, units = "hours")) # in this way I exclude low values but after infusion (considering worst infusion time)
    
    data_plot_map <- quantiles_map %>%
      filter(time <= max(data_plot$time)) %>%
      mutate(Time = dose_info$first_dose_time + as.difftime(time, units = "hours"))
    
    data_plot_dv <- patient_data_for_map %>%
      select(value = DV, time = TIME) %>%
      filter(time <= max(data_plot$time)) %>%
      mutate(Time = dose_info$first_dose_time + as.difftime(time, units = "hours")) %>%
      left_join(data_plot %>% select(Time, Cocentration), by = "Time") %>%
      mutate(
        Error = round((value - Cocentration) / value, 2) * 100, #this is the relative (to the population esitmations) prediction error 
      )
    
    plot_dalba <- plot_ly() %>%
      add_lines(
        data = data_plot,
        x = ~Time, y = ~ Cocentration, name = "Concentration",
        line = list(color = "blue"),
        hoverinfo = "text",
        text = ~ paste("Date: ", format(Time, "%b %d, %H:%M"), "<br>Estimated: ", Cocentration),
        showlegend = TRUE
      ) %>%
      add_ribbons(
        data = data_plot_map,
        x = ~Time, ymin = ~lower90, ymax = ~upper90,
        name = "90% PI",
        line = list(color = "rgba(0,0,255,0.3)"),
        fillcolor = "rgba(0, 0, 255, 0.3)",
        hoverinfo = "text",
        text = ~ paste("Upper 90 PI: ", upper90, "<br>Lower 90 PI: ", lower90),
        showlegend = TRUE
      ) %>%
      add_markers(
        data = data_plot_dv,
        x = ~Time, y = ~value, name = "Patient Data",
        marker = list(size = 12, color = "orange"),
        hoverinfo = "text",
        text = ~ paste("Measured: ", value, "<br>Error MAP: ", Error, "%"),
        showlegend = TRUE
      ) %>%
      layout(
        title = list(text = paste0("Single model, ",best_model)),
        yaxis = list(title = "Concentration (mg/L)", autorange = TRUE, type = "log"),
        xaxis = list(title = "", tickformat = "%b%d"),
        hovermode = "x unified",
        hoverlabel = list(font = list(size = 12)),
        legend = list(
          orientation = "h",
          x = 0.5,
          y = -0.2,
          xanchor = "center",
          yanchor = "bottom"
        )
      )
    
    
    # Plot MoA
    data_prova <- m_dt_out %>% 
      filter(!(Cc <= 0.5 & time > dose_info$max_duration)) %>% 
      mutate(Time = dose_info$first_dose_time + as.difftime(time, units = "hours"),
             MoA_IPRED = round(MoA_IPRED,2),
             Cc = round(Cc,2)) 
    
    plot_moa <- plot_ly() %>%
      add_lines(
        data = data_prova,
        x = ~ Time,
        y = ~ Cc,
        color = ~ model_name,
        name = ~ model_name,
        line = list(width = 2),
        hoverinfo = 'text',
        text = ~ paste("Date: ", format(Time, "%b %d, %H:%M"),"<br>Estimated:", Cc, "<br>Model:", model_name)
      ) %>%
      add_lines(
        data = data_prova,
        x = ~ Time,
        y = ~ MoA_IPRED,
        name = "Model Average",
        line = list(
          width = 3,
          dash = "solid",
          color = 'black'
        ),
        hoverinfo = 'text',
        text = ~ paste("Date: ", format(Time, "%b %d, %H:%M"),"<br>Model Averaging:", MoA_IPRED)
      ) %>%
      add_markers(
        data = data_plot_dv,
        x = ~ Time,
        y = ~ value,
        name = "Patient Data",
        marker = list(size = 12, color = "orange"),
        hoverinfo = "text",
        text = ~ paste("Measured: ", value),
        showlegend = TRUE
      ) %>%
      layout(
        title = list(text = "Model Averaging"),
        yaxis = list(title = "Concentration (mg/L)", autorange = TRUE, type = "log"),
        xaxis = list(title = "", tickformat = "%b%d"),
        hovermode = "x unified",
        hoverlabel = list(font = list(size = 12)),
        showlegend = TRUE,
        legend = list(
          orientation = "h",
          x = 0.5,
          y = -0.2,
          xanchor = "center",
          yanchor = "bottom"
        )
      )
    
    combined_plot <- subplot(
      plot_moa, 
      p_weights,
      widths = c(0.75, 0.25),
      margin = 0.02,
      titleX = TRUE,
      titleY = TRUE
    ) %>%
      layout(
        showlegend = TRUE, 
        legend = list(
          orientation = "h",
          x = 0.5,
          y = -0.2,
          xanchor = "center",
          yanchor = "bottom",
          hovermode = "x unified"
        )
      )
    
    
    # Common part, Pk parameters
    
    covariates_used <- unique(unlist(lapply(model_list, function(m) m$covariates)))
    
    patient_tb <- isolate(patient_data()$absolute) %>%
      filter(EVID == 1) %>%
      dplyr::select(any_of(c("TIME",covariates_used))) %>%
      mutate(
        TIME = format(TIME, "%d/%m/%Y")
      )
    
    # Common part to get times  and targets
    auc_data_plot <- prova %>%
      mutate(
        AUC = (AUC / target_mic * free_drug)
      ) %>%
      filter(time %% 24 == 0) %>% # important to filter time based on 24 hour time/days
      mutate(AUC24 = round(AUC - lag(AUC, default = 0))) %>%
      mutate(
        above_target = if_else(AUC24 > target_auc, "Above", "Below")
      )
    
    auc_data <- auc_data_plot %>% filter(time >= from) # use the last relative time of the last given dose
    
    # Find the last instance where above_target_ab is "Above"
    first_below_target_ab <- auc_data %>%
      filter(above_target == "Above") %>%
      tail(1)
    
    output$concentrationPlot <- renderPlotly({
      combined_plot
    })
    
    
    output$concentrationPlot2 <- renderPlotly({
      plot_dalba
    })

    
    ## Textputput --------------------------------------------------------------
    output$message <- renderUI({
      HTML(paste0(
        "Here you can find the results for the model averaging and single model prediction, ",
        "based on the most influential (highest weighted) model.<br><br>",
        
        "The redosing timing is calculated based on the averaged Bayesian predictions. ",
        "The single model plot shows also a 90% prediction interval derived through the Sequential Importance Resampling algorithm.<br><br>",
        
        "All estimations and simulations were performed using the <strong>posologyr</strong> package."
      ))
    })
    
    
    ## Box ouptut --------------------------------------------------------------
    
    output$box_auc24 <- renderUI({
      date_to_display <- dose_info$first_dose_time + as.difftime(round(first_below_target_ab$time / 24), units = "days")
      days_from_last_dose <- round(first_below_target_ab$time/24) - round(from/24)
      
      box(
        title = tagList(
          icon("clock"),
          HTML("<strong style='font-size: 18px;'>Time to redose according to last 24h <span style='font-family: cursive;'>f</span>AUC<sub>0-24</sub>/MIC</strong>")
        ),
        status = "primary",
        solidHeader = TRUE,
        collapsible = FALSE,
        width = 6,
        height = "auto",
        div(
          # Display the date along with the number of days in parentheses
          h4(
            HTML(sprintf(
              "%s (%d days from last dose)",
              format(date_to_display, "%d %B %Y"), days_from_last_dose
            )),
            style = "margin-bottom: 5px; font-weight: bold; color: #333;"
          ),
          p(
            # Using HTML to format "f" in cursive and subscript for "24"
            HTML(sprintf("When estimated to be: %.2f mg·h/L", first_below_target_ab$AUC24)),
            style = "font-size: 16px; color: #666;" # Font size and color for AUC
          ),
          p(
            # Using HTML to make "min" a subscript
            HTML(sprintf("And total C<sub>min</sub>: %.1f mg/L", first_below_target_ab$Cc)),
            style = "font-size: 16px; color: #666;" # Font size and color for Cmin
          ),
          style = "padding: 20px; text-align: center; background-color: #f9f9f9; border-radius: 5px;" # Background color and rounded corners
        )
      )
    })
    
    # Render UI for the box which contains the data table
    output$box_parameters <- renderUI({
      renderUI({
        box(
          title = tagList(icon("notes-medical", lib = "font-awesome"), " Parameters"),
          status = "primary",
          solidHeader = TRUE,
          collapsible = FALSE,
          width = 12,
          height = "auto",
          div(
            style = "padding: 20px;",
            dataTableOutput("combinedTable")
          )
        )
      })
    })
    
    output$combinedTable <- renderDataTable({
      datatable(patient_tb, options = list(dom = "t")) # 't' for table only without any control elements
    })
    
    # Save report ------------------------------------------------------------
    
    output$save_report <- downloadHandler(
      filename = function() {
        report_name <- ifelse(input$report_name != "", input$report_name, "Dalbavancin Report")
        paste(report_name, "-", Sys.Date(), ".html", sep = "")
      },
      content = function(file) {
        # Use rmarkdown to render the report
        rmarkdown::render(
          "report_template.Rmd",
          output_file = file,
          params = list(
            plot_dalba = plot_dalba,
            combined_plot = combined_plot,
            first_below_target_ab = first_below_target_ab,
            doses = doses,
            patient_tb = patient_tb,
            date_to_display = dose_info$first_dose_time + as.difftime(round(first_below_target_ab$time / 24), units = "days"),
            days_from_last_dose =  round(first_below_target_ab$time/24) - round(from/24),
            free_drug = free_drug,
            target_mic = target_mic,
            target_auc = target_auc
          ),
          envir = new.env(parent = globalenv()) #to clean environment
        )
      }
    )    
  }
}
shinyApp(ui = ui, server = server)