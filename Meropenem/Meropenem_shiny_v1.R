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
library(shinybusy)

options(scipen = 999)
rxode2::setRxThreads(2L)

`%||%` <- function(a, b) if (!is.null(a)) a else b

message(format(Sys.time(), "[%H:%M:%S]"), " Sourcing meropenem_models.R...")
source("meropenem_models.R")
message(format(Sys.time(), "[%H:%M:%S]"), " meropenem_models.R loaded — ",
        length(meropenem_models), " models")

# ── Helper functions ──────────────────────────────────────────────────────────

calculate_ibw <- function(height_cm, gender) {
  h_in <- height_cm / 2.54
  if (gender == "male") 50 + 2.3 * (h_in - 60)
  else 45.5 + 2.3 * (h_in - 60)
}

calculate_lbm <- function(weight_kg, height_cm, gender) {
  bmi <- weight_kg / (height_cm / 100)^2
  if (gender == "male") (9270 * weight_kg) / (6680 + 216 * bmi)
  else (9270 * weight_kg) / (8780 + 244 * bmi)
}

calculate_clcr <- function(age, weight_kg, scr_mg_dl, gender) {
  cl <- (140 - age) * weight_kg / (72 * scr_mg_dl)
  if (gender == "female") cl * 0.85 else cl
}

calculate_egfr_ckd_epi <- function(age, scr_mg_dl, gender) {
  if (gender == "male") { kappa <- 0.9; alpha <- -0.302 }
  else                  { kappa <- 0.7; alpha <- -0.241 }
  ratio <- scr_mg_dl / kappa
  egfr  <- 142 * min(ratio, 1)^alpha * max(ratio, 1)^(-1.200) * 0.9938^age
  if (gender == "female") egfr * 1.012 else egfr
}

calculate_egfr_mdrd <- function(age, scr_mg_dl, gender) {
  186.3 * scr_mg_dl^(-1.154) * age^(-0.203) * ifelse(gender == "female", 0.742, 1)
}

# posologyr generic add/remove observer
meroRegisterAddRemove <- function(input, output, session, input_add, input_remove,
                                  input_time, input_value, var_name, var_col, rv,
                                  time_format = "%d/%m/%Y %H:%M") {
  observeEvent(input[[input_add]], {
    new_entry <- tibble(
      TIME    = as.POSIXct(input[[input_time]], format = time_format),
      !!var_col := input[[input_value]]
    ) %>% distinct(TIME, .keep_all = TRUE)
    if (!any(rv[[var_name]]$TIME == new_entry$TIME))
      rv[[var_name]] <- bind_rows(rv[[var_name]], new_entry)
  })
  observeEvent(input[[input_remove]], {
    if (nrow(rv[[var_name]]) > 0)
      rv[[var_name]] <- rv[[var_name]][-nrow(rv[[var_name]]), ]
  })
}

meroAddRemoveButtons <- function(add_id, remove_id, width = 2) {
  column(width, br(),
    actionButton(add_id,    "+", class = "btn btn-success btn-xs"),
    actionButton(remove_id, "-", class = "btn btn-danger  btn-xs")
  )
}

# ── Model / patient-type mapping ──────────────────────────────────────────────

model_by_type <- list(
  "General ICU" = c(
    "Ehmann_2019", "Eisert_2021", "Gijsen_2022", "Hanberg_2018",
    "Jaruratanasirikul_2015", "Lee_2021", "Minichmayr_2018",
    "Padulles_Zamora_2019", "Roberts_2009", "Yoko_2020", "Troisi_2024",
    "Boonpeng_2022", "Isla_2008"
  ),
  "CRRT / CVVHDF" = c(
    "An_2023", "Burger_2018", "Kang_2022", "Shekar_2014", "Onichimowski_2020"
  ),
  "CVVH + Burns" = c("Selig_2022"),
  "ECMO"         = c("Ha_Lee_2021"),
  "ACLF"         = c("Grensemann_2020"),
  "CNS Infection"= c("Wicha_2024"),
  "Advanced Renal Replacement" = c(
    "OJeanson_2021", "Westermann_2021", "Ulldemolins_2015", "Rancic_2024"
  )
)

calculate_metrics <- function(data, model_name) {
  N         <- nrow(data)
  predicted <- data$Cc
  observed  <- data$Cc_observed
  rBias <- (1/N) * sum((predicted - observed) / observed) * 100
  rRMSE <- sqrt((1/N) * sum(((predicted - observed)^2) / observed^2)) * 100
  tibble(model_name = model_name, rBias = rBias, rRMSE = rRMSE,
         weight = data$weight[1])
}

en_months <- c("January","February","March","April","May","June",
               "July","August","September","October","November","December")

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  dashboardHeader(title = tagList(
    span(class = "logo-lg", "Meropenem"),
    icon("pills", class = "fa-lg")
  )),

  dashboardSidebar(sidebarMenu(
    menuItem("Patient Info",       tabName = "patient_info",       icon = icon("user")),
    menuItem("Dose Regimen",       tabName = "dose_regimen",       icon = icon("pills")),
    menuItem("Model Targets",      tabName = "model_targets",      icon = icon("bullseye")),
    menuItem("Model Specification",tabName = "model_specification",icon = icon("cogs")),
    menuItem("Report / Load / Save",tabName = "report",            icon = icon("file-pdf")),
    hr()
  )),

  dashboardBody(
    useShinyjs(),
    useShinyFeedback(),

    tabItems(

      # ── Patient Info ─────────────────────────────────────────────────────────
      tabItem(tabName = "patient_info", fluidRow(
        column(width = 8,
          box(width = 12, solidHeader = TRUE, status = "primary",
            fluidRow(
              column(4, numericInput("age",    "Age (years)",   value = 65)),
              column(4, selectInput("gender",  "Gender",
                                   choices = c("male","female"), selected = "male")),
              column(4, numericInput("height", "Height (cm)",   value = 170))
            ),

            h4("Patient type"),
            fluidRow(column(6,
              selectInput("patient_type", NULL,
                          choices  = names(model_by_type),
                          selected = "General ICU",
                          width    = "100%")
            )),

            # ── General ICU extra flags ──────────────────────────────────────
            conditionalPanel("input.patient_type == 'General ICU'",
              fluidRow(
                column(3, checkboxInput("flag_sepsis", "Sepsis", value = FALSE)),
                column(3, checkboxInput("flag_shock",  "Shock",  value = FALSE)),
                column(3, numericInput("diur_gen", "Urine output (mL/24h)", value = 1000))
              )
            ),

            # ── CRRT / CVVHDF ────────────────────────────────────────────────
            conditionalPanel("input.patient_type == 'CRRT / CVVHDF'",
              fluidRow(
                column(4, numericInput("qfd",  "Filter flow Qfd (L/h)", value = 2)),
                column(4, checkboxInput("flag_crrt_only", "CRRT mode (not CVVH)", value = TRUE))
              )
            ),

            # ── CVVH + Burns ─────────────────────────────────────────────────
            conditionalPanel("input.patient_type == 'CVVH + Burns'",
              fluidRow(
                column(4, numericInput("cvvh_qf",    "Ultrafiltration Qf (L/h)",  value = 2)),
                column(4, numericInput("cvvh_sc",    "Sieving coefficient Sc",     value = 0.9)),
                column(4, numericInput("cvvh_cf",    "Correction factor CF",       value = 1)),
                column(4, checkboxInput("flag_burn",  "Burns patient",             value = TRUE)),
                column(4, numericInput("tbsa",       "TBSA (%)",                   value = 30)),
                column(4, numericInput("lbm_manual", "LBM (kg, optional override)",value = NA))
              )
            ),

            # ── ECMO ─────────────────────────────────────────────────────────
            conditionalPanel("input.patient_type == 'ECMO'",
              fluidRow(
                column(4, numericInput("lpm", "ECMO flow rate LPM (L/min)", value = 4))
              )
            ),

            # ── ACLF ─────────────────────────────────────────────────────────
            conditionalPanel("input.patient_type == 'ACLF'",
              p("Acute-on-chronic liver failure — no additional covariates required.")
            ),

            # ── CNS Infection ─────────────────────────────────────────────────
            conditionalPanel("input.patient_type == 'CNS Infection'",
              fluidRow(
                column(4, numericInput("il6_csf", "CSF IL-6 (ng/L)", value = 500))
              )
            ),

            # ── Advanced Renal Replacement ────────────────────────────────────
            conditionalPanel("input.patient_type == 'Advanced Renal Replacement'",
              fluidRow(
                column(3, numericInput("diur_rrt",  "Urine output (mL/24h)",        value = 500)),
                column(3, checkboxInput("flag_dia_c",  "Continuous dialysis",        value = FALSE)),
                column(3, checkboxInput("flag_dia_sc", "Semi-continuous dialysis",   value = FALSE)),
                column(3, numericInput("t_dial_start","t_Dial at start (h)",         value = 0)),
                column(3, checkboxInput("flag_van",  "Concomitant vancomycin",       value = FALSE)),
                column(3, checkboxInput("flag_col",  "Concomitant colistin",         value = FALSE)),
                column(3, checkboxInput("flag_hta",  "Hypertension",                 value = FALSE)),
                column(3, numericInput("wbc",        "WBC (×10⁹/L)",                value = 10))
              )
            ),

            hr(),

            h4("Weight"),
            fluidRow(
              column(4, airDatepickerInput("weight_time", "Day of record",
                                          value = Sys.Date(), width = "100%")),
              column(4, numericInput("weight", "Weight (kg)", value = 70)),
              meroAddRemoveButtons("addWeight", "removeWeight", 4)
            ),

            h4("Serum creatinine"),
            fluidRow(
              column(4, airDatepickerInput("creatinine_time", "Day of record",
                                          value = Sys.Date(), width = "100%")),
              column(4, numericInput("creatinine", "Serum Cr (mg/dL)", value = 0.9)),
              meroAddRemoveButtons("addCreatinine", "removeCreatinine", 4)
            ),

            h4("Albumin"),
            fluidRow(
              column(4, airDatepickerInput("albumin_time", "Day of record",
                                          value = Sys.Date(), width = "100%")),
              column(4, numericInput("albumin", "Albumin (g/dL)", value = 4.0)),
              meroAddRemoveButtons("addAlbumin", "removeAlbumin", 4)
            ),

            h4("TDM"),
            fluidRow(
              column(4, airDatepickerInput("tdm_time", "Date/time",
                                          value = Sys.time(), timepicker = TRUE,
                                          timepickerOpts = timepickerOptions(
                                            hoursStep = 1, minutesStep = 10),
                                          width = "100%")),
              column(4, numericInput("tdm", "TDM (mg/L)", value = 0)),
              meroAddRemoveButtons("addTDM", "removeTDM", 4)
            )
          )
        ),

        column(width = 4,
          box(title = "Patient data", width = 12, solidHeader = TRUE,
              status = "primary", collapsible = TRUE,
            DTOutput("TableWT"),        hr(),
            DTOutput("TableCreatinine"),hr(),
            DTOutput("TableAlbumin"),   hr(),
            DTOutput("TableTDM")
          )
        )
      )),

      # ── Dose Regimen ─────────────────────────────────────────────────────────
      tabItem(tabName = "dose_regimen", fluidRow(
        column(width = 6,
          box(width = 12, solidHeader = TRUE, status = "primary",
            h4("Add regimen block"),
            fluidRow(
              column(6, airDatepickerInput("regimen_start", "Start date/time",
                                          value = Sys.time(), timepicker = TRUE,
                                          timepickerOpts = timepickerOptions(
                                            hoursStep = 1, minutesStep = 5),
                                          width = "100%")),
              column(6, numericInput("regimen_dose", "Dose (mg)", value = 1000))
            ),
            fluidRow(
              column(4, selectInput("regimen_interval", "Interval (h)",
                                   choices = c(6, 8, 12, 24), selected = 8)),
              column(4, numericInput("regimen_dur", "Infusion duration (h)",
                                    value = 3, min = 0.5, max = 8, step = 0.5)),
              column(4, numericInput("regimen_n", "Number of doses", value = 6))
            ),
            fluidRow(
              column(6, actionButton("addRegimen",    "Add regimen",
                                    class = "btn btn-success", style = "width:100%")),
              column(6, actionButton("removeRegimen", "Remove last",
                                    class = "btn btn-danger",  style = "width:100%"))
            )
          )
        ),
        column(width = 6,
          box(title = "Regimen blocks", width = 12, solidHeader = TRUE,
              status = "primary",
            DTOutput("regimenTable"),
            hr(),
            h5("Expanded doses"),
            DTOutput("expandedDoseTable")
          )
        )
      )),

      # ── Model Targets ─────────────────────────────────────────────────────────
      tabItem(tabName = "model_targets", fluidRow(column(width = 12,
        box(width = 12, solidHeader = TRUE, status = "primary",
          fluidRow(
            column(4,
              numericInput("target_mic",     "MIC (mg/L)",              value = 2),
              numericInput("mic_multiplier", "MIC multiplier (× MIC)",  value = 1),
              numericInput("target_pct",     "Target %fT>xMIC",         value = 40),
              numericInput("protein_binding","Protein binding (%)",      value = 2),
              actionButton("reset_target", "Reset",
                           class = "btn-primary", style = "width:100%;color:white;")
            ),
            column(8,
              tags$h4("Explanation of targets"),
              tags$p("The target is the percentage of the dosing interval during which the free ",
                     "meropenem concentration exceeds x × MIC."),
              tags$p("For static kill, a typical target is 40% fT>MIC; for bactericidal effect ",
                     "100% fT>MIC (or fT>4×MIC) is often recommended in critically ill patients."),
              tags$p("Meropenem protein binding is approximately 2%, so the free fraction is ~0.98. ",
                     "Extended infusions (3–4 h) substantially improve target attainment.")
            )
          )
        )
      ))),

      # ── Model Specification ───────────────────────────────────────────────────
      tabItem(tabName = "model_specification",
        box(width = 12, solidHeader = TRUE, status = "primary", collapsible = FALSE,
          fluidRow(
            column(4,
              tags$label("Models included in model averaging"),
              uiOutput("model_checkboxes")
            ),
            column(8,
              tags$p("Models are pre-selected based on the patient type and available covariates. ",
                     "You may override the selection here.")
            )
          )
        )
      ),

      # ── Report ────────────────────────────────────────────────────────────────
      tabItem(tabName = "report",
        box(title = "Report management", width = 12, solidHeader = TRUE,
            status = "primary",
          textInput("report_name", "Report name", value = "Meropenem Report"),
          fluidRow(
            column(4, downloadButton("save_state", "Save state",
                                    icon = icon("save"), class = "btn-primary",
                                    style = "width:100%;color:white;")),
            column(4, fileInput("load_state", label = NULL,
                                buttonLabel = "Load state", accept = ".zip",
                                width = "100%"))
          )
        )
      )
    ),

    # ── Always-visible output section ─────────────────────────────────────────
    fluidRow(column(width = 12,
      box(width = 12, solidHeader = TRUE, status = "primary", collapsible = FALSE,
        actionButton("calculate_button", "Run calculation",
                     class = "btn btn-primary btn-lg"),
        br(), br(),
        uiOutput("message"),
        br(),
        plotlyOutput("concentrationPlot",  height = "400px"),
        br(),
        plotlyOutput("concentrationPlot2", height = "400px"),
        br(),
        plotlyOutput("tmic_plot",          height = "300px"),
        br(),
        uiOutput("box_tmic_summary"),
        br(),
        uiOutput("box_dose_rec"),
        br(),
        uiOutput("box_parameters")
      )
    )),

    scrollToTop = TRUE,
    preloader  = list(html = spin_1(), color = "gray")
  )
)


# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # ── Disclaimer modal ────────────────────────────────────────────────────────
  showModal(modalDialog(
    title = "Meropenem Precision Dosing Tool",
    p("Disclaimer: This tool is for informational purposes only and intended for ",
      "use by health care professionals. It is not a substitute for professional ",
      "medical advice, diagnosis, or treatment."),
    footer = tagList(fluidRow(
      column(3, checkboxInput("accept",  "I understand", value = FALSE)),
      column(9, actionButton("confirm", "Confirm", class = "btn-primary"))
    )),
    easyClose = FALSE, size = "m"
  ))

  shinyjs::disable("confirm")
  observeEvent(input$accept,  {
    if (input$accept) shinyjs::enable("confirm") else shinyjs::disable("confirm")
  })
  observeEvent(input$confirm, removeModal())

  # ── Reactive values ─────────────────────────────────────────────────────────
  rv <- reactiveValues(
    weights     = tibble(TIME = as.POSIXct(character()), WT        = numeric()),
    creatinines = tibble(TIME = as.POSIXct(character()), CREATININE= numeric()),
    albumins    = tibble(TIME = as.POSIXct(character()), ALBUMIN   = numeric()),
    tdm         = tibble(TIME = as.POSIXct(character()), DV        = numeric(),
                         DVID = character()),
    regimens    = tibble(
      START_TIME = as.POSIXct(character()),
      AMT        = numeric(),
      INTERVAL_H = numeric(),
      DUR_H      = numeric(),
      N_DOSES    = integer()
    )
  )

  # ── Time-varying covariate handlers ─────────────────────────────────────────
  purrr::walk(list(
    list(add="addWeight",     remove="removeWeight",     time="weight_time",
         value="weight",      var="weights",     col="WT"),
    list(add="addCreatinine", remove="removeCreatinine", time="creatinine_time",
         value="creatinine",  var="creatinines", col="CREATININE"),
    list(add="addAlbumin",    remove="removeAlbumin",    time="albumin_time",
         value="albumin",     var="albumins",    col="ALBUMIN")
  ), function(cfg) {
    meroRegisterAddRemove(input, output, session,
      input_add    = cfg$add,   input_remove = cfg$remove,
      input_time   = cfg$time,  input_value  = cfg$value,
      var_name     = cfg$var,   var_col      = cfg$col,
      rv           = rv,        time_format  = "%d/%m/%Y")
  })

  observeEvent(input$addTDM, {
    new_entry <- tibble(
      TIME = as.POSIXct(input$tdm_time, format = "%d/%m/%Y %H:%M"),
      DV   = input$tdm,
      DVID = "Cc"
    ) %>% distinct(TIME, .keep_all = TRUE)
    if (!any(rv$tdm$TIME == new_entry$TIME))
      rv$tdm <- bind_rows(rv$tdm, new_entry)
  })
  observeEvent(input$removeTDM, {
    if (nrow(rv$tdm) > 0) rv$tdm <- rv$tdm[-nrow(rv$tdm), ]
  })

  # ── Regimen handlers ─────────────────────────────────────────────────────────
  observeEvent(input$addRegimen, {
    new_reg <- tibble(
      START_TIME = as.POSIXct(input$regimen_start, format = "%d/%m/%Y %H:%M"),
      AMT        = as.numeric(input$regimen_dose),
      INTERVAL_H = as.numeric(input$regimen_interval),
      DUR_H      = input$regimen_dur,
      N_DOSES    = as.integer(input$regimen_n)
    )
    rv$regimens <- bind_rows(rv$regimens, new_reg) %>% arrange(START_TIME)
  })
  observeEvent(input$removeRegimen, {
    if (nrow(rv$regimens) > 0)
      rv$regimens <- rv$regimens[-nrow(rv$regimens), ]
  })

  output$regimenTable <- renderDT({
    tb <- rv$regimens %>%
      mutate(START_TIME = format(START_TIME, "%d/%m/%Y %H:%M"))
    datatable(tb, options = list(dom = "t"), rownames = FALSE)
  }, server = FALSE)

  # Expand regimens to individual dose rows
  expanded_doses <- reactive({
    if (nrow(rv$regimens) == 0) return(tibble())
    purrr::map_dfr(seq_len(nrow(rv$regimens)), function(i) {
      r <- rv$regimens[i, ]
      times <- r$START_TIME + as.difftime(
        (seq_len(r$N_DOSES) - 1) * r$INTERVAL_H, units = "hours")
      tibble(TIME = times, AMT = r$AMT, DUR = r$DUR_H, ADDL = 0L, II = 0)
    }) %>% distinct(TIME, .keep_all = TRUE) %>% arrange(TIME)
  })

  output$expandedDoseTable <- renderDT({
    tb <- expanded_doses() %>%
      mutate(TIME = format(TIME, "%d/%m/%Y %H:%M"))
    datatable(tb, options = list(dom = "t", pageLength = 6), rownames = FALSE)
  }, server = FALSE)

  # ── Table outputs ────────────────────────────────────────────────────────────
  output$TableWT <- renderDT({
    datatable(rv$weights %>% mutate(TIME = format(TIME, "%d/%m/%Y")),
              options = list(dom = "t"), rownames = FALSE)
  }, server = FALSE)

  output$TableCreatinine <- renderDT({
    datatable(rv$creatinines %>% mutate(TIME = format(TIME, "%d/%m/%Y")),
              options = list(dom = "t"), rownames = FALSE)
  }, server = FALSE)

  output$TableAlbumin <- renderDT({
    datatable(rv$albumins %>% mutate(TIME = format(TIME, "%d/%m/%Y")),
              options = list(dom = "t"), rownames = FALSE)
  }, server = FALSE)

  output$TableTDM <- renderDT({
    datatable(rv$tdm %>% select(TIME, DV) %>%
                mutate(TIME = format(TIME, "%d/%m/%Y %H:%M")),
              options = list(dom = "t"), rownames = FALSE)
  }, server = FALSE)

  # ── Model selection based on patient type ────────────────────────────────────
  selected_model_names <- reactive({
    model_by_type[[input$patient_type]] %||%
      model_by_type[["General ICU"]]
  })

  output$model_checkboxes <- renderUI({
    nms <- selected_model_names()
    checkboxGroupButtons(
      inputId   = "model_selection",
      label     = NULL,
      selected  = nms,
      choices   = names(meropenem_models),
      size      = "sm",
      direction = "vertical",
      justified = TRUE,
      checkIcon = list(yes = icon("ok", lib = "glyphicon"))
    )
  })

  # Update checkboxes when patient type changes
  observeEvent(input$patient_type, {
    updateCheckboxGroupButtons(session, "model_selection",
                               selected = selected_model_names())
  })

  active_models <- reactive({
    message(format(Sys.time(), "[%H:%M:%S]"), " active_models: reading input$model_selection...")
    selection <- input$model_selection
    # renderUI populates model_selection asynchronously; fall back to the
    # patient-type default when the input hasn't initialised yet
    if (is.null(selection) || length(selection) == 0) {
      message(format(Sys.time(), "[%H:%M:%S]"), " active_models: input NULL, using patient-type default")
      selection <- model_by_type[[input$patient_type]] %||% model_by_type[["General ICU"]]
    }
    message(format(Sys.time(), "[%H:%M:%S]"), " active_models: selection = ",
            paste(selection, collapse = ", "))
    result <- meropenem_models[names(meropenem_models) %in% selection]
    message(format(Sys.time(), "[%H:%M:%S]"), " active_models: OK — ", length(result), " models")
    result
  })

  # ── Reset targets ────────────────────────────────────────────────────────────
  observeEvent(input$reset_target, {
    updateNumericInput(session, "target_mic",     value = 2)
    updateNumericInput(session, "mic_multiplier", value = 1)
    updateNumericInput(session, "target_pct",     value = 40)
    updateNumericInput(session, "protein_binding",value = 2)
  })

  # ── Patient data reactive ────────────────────────────────────────────────────
  patient_data <- reactive({
    message(format(Sys.time(), "[%H:%M:%S]"), " patient_data: START")
    age    <- input$age
    gender <- input$gender
    height <- input$height

    ibw <- calculate_ibw(height, gender)

    # Collect all scalar covariates once (outside rowwise)
    pt   <- input$patient_type
    diur_val <- if (pt == "General ICU") (input$diur_gen %||% 1000)
                else if (pt %in% c("Advanced Renal Replacement","CRRT / CVVHDF"))
                  (input$diur_rrt %||% 1000)
                else 1000

    cl_cvvh_val <- if (pt == "CVVH + Burns" && !is.null(input$cvvh_qf))
                     input$cvvh_qf * (input$cvvh_sc %||% 0.9) * (input$cvvh_cf %||% 1)
                   else 0

    message(format(Sys.time(), "[%H:%M:%S]"), " patient_data: building cov_rows (rowwise)...")
    cov_rows <- bind_rows(rv$weights, rv$creatinines, rv$albumins) %>%
      mutate(EVID = 1, DV = as.numeric(NA)) %>%
      arrange(TIME) %>%
      fill(WT, CREATININE, ALBUMIN, .direction = "downup") %>%
      distinct() %>%
      rowwise() %>%
      mutate(
        TBW             = WT,
        IBW             = round(ibw, 1),
        LBM             = round(calculate_lbm(WT, height, gender), 1),
        CLCR_TBW        = round(calculate_clcr(age, WT,  CREATININE, gender)),
        CLCR_CG         = CLCR_TBW,
        CLCR_CG_IBW     = round(calculate_clcr(age, ibw, CREATININE, gender)),
        eGFR_CKD_EPI    = round(calculate_egfr_ckd_epi(age, CREATININE, gender)),
        eGFR_MDRD       = round(calculate_egfr_mdrd(age, CREATININE, gender)),
        CREA_micromol_l = round(CREATININE * 88.42),
        ALB             = ALBUMIN,
        AGE             = age
      ) %>%
      ungroup() %>%
      mutate(
        CRRT    = as.integer(pt == "CRRT / CVVHDF"),
        RRT     = as.integer(pt %in% c("CRRT / CVVHDF","Advanced Renal Replacement")),
        CVVH    = as.integer(pt == "CVVH + Burns"),
        ACLF    = as.integer(pt == "ACLF"),
        SEPSIS  = as.integer(isTRUE(input$flag_sepsis)),
        shock   = as.integer(isTRUE(input$flag_shock)),
        BURN    = as.integer(isTRUE(input$flag_burn)),
        HTA     = as.integer(isTRUE(input$flag_hta)),
        VAN     = as.integer(isTRUE(input$flag_van)),
        COL     = as.integer(isTRUE(input$flag_col)),
        RFS     = as.integer(eGFR_MDRD >= 120),
        DIA_C   = as.integer(isTRUE(input$flag_dia_c)),
        DIA_SC  = as.integer(isTRUE(input$flag_dia_sc)),
        Qfd     = input$qfd     %||% 0,
        LPM     = input$lpm     %||% 0,
        TBSA    = input$tbsa    %||% 0,
        CL_CVVH = cl_cvvh_val,
        IL6_CSF = input$il6_csf %||% 500,
        WBC     = input$wbc     %||% 10,
        DIUR    = diur_val,
        t_Dial  = input$t_dial_start %||% 0
      )
    message(format(Sys.time(), "[%H:%M:%S]"), " patient_data: cov_rows OK — ", nrow(cov_rows), " rows")

    doses    <- expanded_doses()
    tdm_data <- rv$tdm %>% mutate(EVID = 0)
    message(format(Sys.time(), "[%H:%M:%S]"), " patient_data: doses = ", nrow(doses),
            ", tdm_data = ", nrow(tdm_data))

    if (nrow(doses) == 0) return(list(absolute = cov_rows))

    first_time <- min(doses$TIME)

    to_relative <- function(tbl) {
      tbl %>%
        mutate(TIME = round(as.numeric(difftime(TIME, first_time, units = "hours")), 1))
    }

    cov_rel  <- to_relative(cov_rows)
    dose_rel <- to_relative(doses)
    tdm_rel  <- to_relative(tdm_data)

    last_dose_rel   <- max(dose_rel$TIME)
    last_interval_h <- if (nrow(rv$regimens) > 0)
                         tail(rv$regimens$INTERVAL_H, 1)
                       else 8

    time_seq_fit  <- sort(unique(c(dose_rel$TIME,
                                   seq(0, last_dose_rel + 48, by = 1))))
    time_seq_pred <- sort(unique(c(dose_rel$TIME,
                                   seq(0, last_dose_rel + last_interval_h, by = 1))))
    message(format(Sys.time(), "[%H:%M:%S]"), " patient_data: last_dose_rel = ", last_dose_rel,
            " h | time_seq_fit length = ", length(time_seq_fit),
            " | time_seq_pred length = ", length(time_seq_pred))

    make_map_data <- function(time_seq) {
      bind_rows(cov_rel, dose_rel, tdm_rel) %>%
        filter(!is.na(TIME)) %>%
        complete(TIME = time_seq) %>%
        arrange(TIME) %>%
        fill(TBW, WT, IBW, LBM, CLCR_TBW, CLCR_CG, CLCR_CG_IBW,
             eGFR_CKD_EPI, eGFR_MDRD, CREA_micromol_l,
             ALB, ALBUMIN, CREATININE,
             CRRT, RRT, CVVH, ACLF, SEPSIS, shock, BURN, HTA, VAN, COL,
             RFS, DIA_C, DIA_SC, Qfd, LPM, TBSA, CL_CVVH, IL6_CSF, WBC,
             DIUR, t_Dial, AGE, .direction = "downup") %>%
        mutate(
          AMT  = if_else(is.na(AMT),  0, AMT),
          DUR  = if_else(is.na(DUR),  0, DUR),
          EVID = if_else(is.na(EVID), 1, EVID)
        ) %>%
        filter(TIME >= 0) %>%
        # Ensure dose rows sort before observations at the same time, then
        # compute OCC: increments by 1 at each real dose (EVID=1, AMT>0).
        # Required by posologyr for models with pi_matrix (IOV) terms.
        arrange(TIME, desc(AMT)) %>%
        mutate(
          OCC = cumsum(EVID == 1 & AMT > 0),
          ID  = 1L          # required by posologyr for IOV (pi_matrix) models
        )
    }

    message(format(Sys.time(), "[%H:%M:%S]"), " patient_data: building pmap_fit...")
    pmap_fit <- make_map_data(time_seq_fit)
    message(format(Sys.time(), "[%H:%M:%S]"), " patient_data: pmap_fit OK — ", nrow(pmap_fit), " rows")

    message(format(Sys.time(), "[%H:%M:%S]"), " patient_data: building pmap_pred...")
    pmap_pred <- make_map_data(time_seq_pred)
    message(format(Sys.time(), "[%H:%M:%S]"), " patient_data: pmap_pred OK — ", nrow(pmap_pred), " rows")

    message(format(Sys.time(), "[%H:%M:%S]"), " patient_data: DONE")
    list(
      absolute           = cov_rows,
      first_dose_time    = first_time,
      last_dose_rel      = last_dose_rel,
      last_interval_h    = last_interval_h,
      dose_times_rel     = dose_rel$TIME,
      dose_intervals     = dose_rel %>% select(TIME, DUR),
      patient_map_fit    = pmap_fit,
      patient_map_pred   = pmap_pred
    )
  })

  # ── Regimen interval (used for %T>MIC) ───────────────────────────────────────
  median_interval <- reactive({
    if (nrow(rv$regimens) == 0) return(8)
    median(rv$regimens$INTERVAL_H)
  })

  # ── Core calculation ─────────────────────────────────────────────────────────
  observeEvent(input$calculate_button, {
    shinybusy::show_modal_spinner(text = "Running MAP estimation…")
    core_function()
    shinybusy::remove_modal_spinner()
  })

  # Timestamped console logger — visible in RStudio console in real-time
  msg <- function(...) message(format(Sys.time(), "[%H:%M:%S]"), " ", ...)

  core_function <- function() {
    set.seed(123)

    msg("=== core_function START ===")

    msg("Isolating expanded_doses()...")
    doses <- isolate(expanded_doses())
    msg("expanded_doses() OK — ", nrow(doses), " rows")

    msg("Isolating active_models()...")
    model_list       <- isolate(active_models())
    model_list_names <- as.list(names(model_list))
    msg("active_models() OK — ", length(model_list), " models")

    msg("Isolating rv$tdm...")
    tdm <- isolate(rv$tdm)
    msg("rv$tdm OK — ", nrow(tdm), " rows")

    msg("Isolating patient_data()...")
    pd <- isolate(patient_data())
    msg("patient_data() OK")

    msg("Isolating scalar inputs...")
    free_frac    <- 1 - (isolate(input$protein_binding) / 100)
    target_mic   <- isolate(input$target_mic)
    mic_mult     <- isolate(input$mic_multiplier)
    target_pct   <- isolate(input$target_pct)
    tau          <- isolate(median_interval())
    first_time      <- pd$first_dose_time
    last_dose_r     <- pd$last_dose_rel
    last_interval_h <- pd$last_interval_h
    dose_times_r    <- pd$dose_times_rel

    msg("Inputs ready. Doses: ", nrow(doses),
        " | TDM points: ", nrow(tdm),
        " | Models selected: ", length(model_list), " (",
        paste(names(model_list), collapse = ", "), ")")

    if (nrow(doses) == 0)
      return(shinyalert("Error", "No doses entered.", type = "error"))
    if (nrow(rv$weights) == 0)
      return(shinyalert("Error", "No weight data entered.", type = "error"))
    if (nrow(rv$creatinines) == 0)
      return(shinyalert("Error", "No creatinine data entered.", type = "error"))
    if (nrow(tdm) == 0)
      return(shinyalert("Error", "No TDM data entered.", type = "error"))
    if (length(model_list) == 0)
      return(shinyalert("Error", "No models selected.", type = "error"))
    if (nrow(tdm) > 0 && min(tdm$TIME) < min(doses$TIME))
      return(shinyalert("Error", "TDM is before the first dose.", type = "error"))

    pmap_fit  <- pd$patient_map_fit
    pmap_pred <- pd$patient_map_pred
    msg("Patient data assembled. pmap_fit rows: ", nrow(pmap_fit),
        " | pmap_pred rows: ", nrow(pmap_pred))

    # Add IOV (KAPPA_*) columns required by models with pi_matrix.
    # posologyr expects these columns in the data; initialise to 0 so it
    # estimates them from the pi_matrix rather than warning about all-NA values.
    iov_cols <- unique(unlist(lapply(model_list, function(m) {
      if (!is.null(m$pi_matrix)) rownames(m$pi_matrix) else NULL
    })))
    for (kappa in iov_cols) {
      if (!kappa %in% names(pmap_fit))  pmap_fit[[kappa]]  <- 0
      if (!kappa %in% names(pmap_pred)) pmap_pred[[kappa]] <- 0
    }
    if (length(iov_cols) > 0)
      msg("IOV columns initialised to 0: ", paste(iov_cols, collapse = ", "))

    # ── MAP for all models ────────────────────────────────────────────────────
    m_out      <- list()
    eta_list   <- list()   # IIV ETAs per model — used for dose recommendation
    inits_list <- list()   # compartment amounts at last_dose_r — starting point for recommendation
    for (i in seq_along(model_list)) {
      nm <- model_list_names[[i]]
      msg("MAP [", i, "/", length(model_list), "] START: ", nm)
      withCallingHandlers(
        tryCatch({
          mi <- poso_estim_map(pmap_fit, model_list[[i]],
                               return_model = TRUE, return_ofv = TRUE)
          msg("MAP [", i, "/", length(model_list), "] DONE:  ", nm,
              " | OFV = ", round(mi$ofv, 2))
          dt <- data.table::data.table(mi$model)
          m_out[[i]] <- cbind(
            dt[, .(time, Cc, AUC)],
            ofv        = mi$ofv,
            LL         = exp(-0.5 * mi$ofv),
            model_name = nm,
            model_number = i
          )
          # Save IIV ETAs only (exclude KAPPA_* occasion-level effects)
          eta_list[[nm]] <- mi$eta[!grepl("^KAPPA_", names(mi$eta))]

          # Save compartment amounts at last_dose_r as starting conditions for
          # the dose recommendation (Option B: start from current patient state)
          tryCatch({
            state_vars <- rxode2::rxModelVars(model_list[[i]]$ppk_model)$state
            model_df   <- as.data.frame(mi$model)
            avail      <- state_vars[state_vars %in% names(model_df)]
            idx        <- which.min(abs(model_df$time - last_dose_r))
            last_state <- setNames(as.numeric(model_df[idx, avail]), avail)
            last_state[grepl("(?i)^AUC", names(last_state))] <- 0  # reset AUC to 0
            inits_list[[nm]] <- last_state
          }, error = function(e) NULL)
        }, error = function(e) {
          msg("MAP [", i, "/", length(model_list), "] ERROR: ", nm, " — ", conditionMessage(e))
          NULL
        }),
        # Suppress cosmetic posologyr warnings about KAPPA_* initialisation.
        # These arise because posologyr initialises IOV parameters from the
        # pi_matrix prior; they do not affect the MAP estimates.
        warning = function(w) {
          if (grepl("has only 'NA' values", conditionMessage(w)))
            invokeRestart("muffleWarning")
        }
      )
    }

    msg("MAP loop finished. Successful models: ",
        sum(!sapply(m_out, is.null)), "/", length(model_list))

    m_out <- Filter(Negate(is.null), m_out)
    if (length(m_out) == 0)
      return(shinyalert("Error", "All models failed. Check covariates.", type = "error"))

    msg("Computing model weights and averaging...")
    m_dt <- data.table::as.data.table(do.call(rbind, m_out))
    m_dt[, weight := LL / sum(unique(LL))]
    m_dt[, MoA_Cc  := sum(Cc * weight), by = time]
    m_dt[, MoA_AUC := sum(AUC * weight), by = time]

    prova <- m_dt %>%
      select(time, Cc = MoA_Cc, AUC = MoA_AUC) %>%
      distinct()

    # ── Metrics + best model ──────────────────────────────────────────────────
    msg("Computing metrics and selecting best model...")
    data_metrics <- m_dt %>%
      left_join(pmap_fit %>% select(Cc_observed = DV, time = TIME),
                by = "time", relationship = "many-to-many") %>%
      filter(!is.na(Cc_observed)) %>%
      split(f = .$model_name)

    metrics_summary <- imap_dfr(data_metrics, calculate_metrics)
    best_model <- metrics_summary %>%
      filter(weight == max(weight)) %>% pull(model_name) %>% unique() %>% head(1)
    msg("Best model: ", best_model)

    # ── Weights bar ───────────────────────────────────────────────────────────
    p_weights <- plot_ly(
      data = metrics_summary,
      x = ~"Weight", y = ~weight, color = ~model_name,
      colors = "Set2", type = "bar",
      hoverinfo = "text", showlegend = FALSE,
      text = ~paste(model_name, ":", round(weight, 2))
    ) %>% layout(
      barmode = "stack",
      yaxis   = list(title = "", range = c(0, 1.01), side = "right"),
      xaxis   = list(title = ""),
      showlegend = FALSE
    )

    # ── SIR for best model ────────────────────────────────────────────────────
    msg("SIR START: ", best_model,
        " (n_sample = 1000, n_resample = 100, pred horizon = ",
        round(last_dose_r + last_interval_h), " h)")
    quantiles_map <- tryCatch({
      sir_map <- poso_estim_sir(pmap_fit, model_list[[best_model]],
                                n_sample = 1e3, n_resample = 1e2)
      msg("SIR sampling done. Building event table...")
      et <- rxode2::as.et(pmap_pred)
      et$add_sampling(seq(0, last_dose_r + last_interval_h, by = 1))
      msg("poso_replace_et START...")
      ci_map <- poso_replace_et(sir_map$model, model_list[[best_model]], et,
                                 interpolation = "locf")
      msg("poso_replace_et DONE. Rows: ", nrow(ci_map))
      ci_map %>%
        group_by(time) %>%
        summarise(
          lower90 = quantile(Cc, 0.05, na.rm = TRUE),
          upper90 = quantile(Cc, 0.95, na.rm = TRUE),
          .groups = "drop"
        )
    }, error = function(e) {
      msg("SIR ERROR (", best_model, "): ", conditionMessage(e),
          " — plot will show MAP only, no CI ribbon")
      NULL
    })
    msg("Quantiles: ", if (is.null(quantiles_map)) "unavailable (SIR failed)" else "computed")

    # ── Build plot data ───────────────────────────────────────────────────────
    msg("Building plots...")
    result <- m_dt %>% filter(model_name == best_model)

    data_plot <- result %>%
      filter(time <= last_dose_r + last_interval_h) %>%
      mutate(Time = first_time + as.difftime(time, units = "hours"),
             Cc   = round(Cc, 2))

    data_plot_map <- if (!is.null(quantiles_map))
      quantiles_map %>%
        filter(time <= max(data_plot$time)) %>%
        mutate(Time = first_time + as.difftime(time, units = "hours"))
    else
      NULL

    data_plot_dv <- pmap_fit %>%
      select(value = DV, time = TIME) %>%
      filter(!is.na(value)) %>%
      mutate(Time = first_time + as.difftime(time, units = "hours"))

    # Single best-model plot
    plot_best <- plot_ly() %>%
      add_lines(data = data_plot, x = ~Time, y = ~Cc,
                name = "Concentration", line = list(color = "blue"),
                hoverinfo = "text",
                text = ~paste("Date:", format(Time, "%b %d %H:%M"),
                              "<br>Estimated:", Cc))

    if (!is.null(quantiles_map)) {
      plot_best <- plot_best %>%
        add_ribbons(data = data_plot_map %>%
                      mutate(lower90 = round(lower90, 2),
                             upper90 = round(upper90, 2)),
                    x = ~Time,
                    ymin = ~lower90, ymax = ~upper90,
                    name = "90% PI",
                    line = list(color = "rgba(0,0,255,0.3)"),
                    fillcolor = "rgba(0,0,255,0.3)",
                    hoverinfo = "text",
                    text = ~paste("90% PI:", lower90, "–", upper90, "mg/L"))
    }

    plot_best <- plot_best %>%
      add_markers(data = data_plot_dv, x = ~Time, y = ~value,
                  name = "TDM",
                  marker = list(size = 10, color = "orange"),
                  hoverinfo = "text",
                  text = ~paste("Measured:", value)) %>%
      layout(
        title  = list(text = paste0("Best model: ", best_model,
                                    if (is.null(quantiles_map)) " (no CI — SIR failed)" else "")),
        yaxis  = list(title = "Concentration (mg/L)"),
        xaxis  = list(title = "", tickformat = "%b %d"),
        hovermode = "x unified",
        legend = list(orientation = "h", x = 0.5, y = -0.2,
                      xanchor = "center", yanchor = "bottom")
      )

    # Model-averaging plot
    data_moa <- m_dt %>%
      filter(time <= last_dose_r + last_interval_h) %>%
      mutate(Time     = first_time + as.difftime(time, units = "hours"),
             MoA_Cc   = round(MoA_Cc, 2),
             Cc       = round(Cc, 2))

    plot_moa <- plot_ly() %>%
      add_lines(data = data_moa, x = ~Time, y = ~Cc,
                color = ~model_name, name = ~model_name,
                line = list(width = 1.5),
                hoverinfo = "text",
                text = ~paste("Date:", format(Time, "%b %d %H:%M"),
                              "<br>Cc:", Cc, "<br>Model:", model_name)) %>%
      add_lines(data = data_moa %>% distinct(Time, MoA_Cc),
                x = ~Time, y = ~MoA_Cc, name = "Model Average",
                line = list(width = 3, color = "black"),
                hoverinfo = "text",
                text = ~paste("Date:", format(Time, "%b %d %H:%M"),
                              "<br>MoA:", MoA_Cc)) %>%
      add_markers(data = data_plot_dv, x = ~Time, y = ~value,
                  name = "TDM", marker = list(size = 10, color = "orange")) %>%
      layout(
        title  = list(text = "Model averaging"),
        yaxis  = list(title = "Concentration (mg/L)"),
        xaxis  = list(title = "", tickformat = "%b %d"),
        hovermode = "x unified",
        legend = list(orientation = "h", x = 0.5, y = -0.2,
                      xanchor = "center", yanchor = "bottom")
      )

    combined_plot <- subplot(plot_moa, p_weights,
                             widths = c(0.75, 0.25), margin = 0.02,
                             titleX = TRUE, titleY = TRUE) %>%
      layout(showlegend = TRUE,
             legend = list(orientation = "h", x = 0.5, y = -0.2,
                           xanchor = "center", yanchor = "bottom"))

    # ── %T>MIC calculation ────────────────────────────────────────────────────
    threshold <- free_frac * mic_mult * target_mic

    tmic_data <- purrr::map_dfr(seq_along(dose_times_r), function(i) {
      t_start <- dose_times_r[i]
      t_end   <- if (i < length(dose_times_r)) dose_times_r[i + 1]
                 else t_start + tau
      pts     <- prova %>% filter(time >= t_start, time < t_end)
      if (nrow(pts) < 2) return(NULL)
      tibble(
        dose_num = i,
        Time     = first_time + as.difftime(t_start, units = "hours"),
        pct_tmic = mean(pts$Cc > threshold) * 100
      )
    })

    plot_tmic <- plot_ly(data = tmic_data, x = ~Time, y = ~pct_tmic,
                         type = "bar",
                         marker = list(color = ~ifelse(pct_tmic >= target_pct,
                                                       "steelblue", "tomato")),
                         hoverinfo = "text",
                         text = ~paste("Interval starting:", format(Time, "%b %d %H:%M"),
                                       "<br>%fT>xMIC:", round(pct_tmic, 1), "%")) %>%
      add_lines(x = ~Time, y = target_pct, name = "Target",
                line = list(color = "black", dash = "dash")) %>%
      layout(
        title = list(text = paste0("fT>", mic_mult, "×MIC (",
                                   round(free_frac * mic_mult * target_mic, 2),
                                   " mg/L) per dosing interval")),
        yaxis = list(title = "%fT>xMIC", range = c(0, 105)),
        xaxis = list(title = "", tickformat = "%b %d"),
        showlegend = FALSE
      )

    last_tmic <- if (nrow(tmic_data) > 0) tail(tmic_data, 1) else NULL

    # ── Covariate table ───────────────────────────────────────────────────────
    cov_used <- unique(unlist(lapply(model_list, function(m) m$covariates)))
    patient_tb <- isolate(patient_data()$absolute) %>%
      select(any_of(c("TIME", cov_used, "WT","CREATININE","ALBUMIN"))) %>%
      mutate(TIME = format(TIME, "%d/%m/%Y"))

    # ── Dose recommendation ───────────────────────────────────────────────────
    msg("Computing dose recommendations...")

    # Most-recent covariate values for forward simulation (exclude bookkeeping cols)
    kappa_drop <- grep("^KAPPA_", names(pmap_fit), value = TRUE)
    cov_vals <- pmap_fit %>%
      tail(1) %>%
      select(-any_of(c("TIME","AMT","DUR","EVID","DV","DVID",
                        "OCC","ID","ADDL","II", kappa_drop))) %>%
      select(where(is.numeric)) %>%
      unlist(use.names = TRUE)

    # Per-model parameter vectors: theta + IIV ETAs + patient covariates
    params_per_model <- setNames(
      lapply(names(model_list), function(nm) {
        eta <- eta_list[[nm]]
        if (is.null(eta)) return(NULL)
        c(model_list[[nm]]$theta, eta, cov_vals)
      }),
      names(model_list)
    )

    # Simulate %fT>MIC starting from the patient's current compartment levels
    # (Option B: inits = compartment amounts at last_dose_r from MAP simulation)
    sim_tmic_ss <- function(ppk_model, params_all, dose_mg, interval_h, dur_h,
                             inits = NULL) {
      n <- as.integer(ceiling(48 / interval_h))  # 2 days regardless of interval
      et <- rxode2::et(amt = dose_mg, ii = interval_h, addl = n - 1L, dur = dur_h) %>%
        rxode2::add.sampling(seq(0, n * interval_h, by = 1))
      sim <- tryCatch(
        suppressWarnings(as.data.frame(
          rxode2::rxSolve(ppk_model, params = params_all, events = et,
                          inits = inits)
        )),
        error = function(e) NULL
      )
      if (is.null(sim) || !"Cc" %in% names(sim)) return(NA_real_)
      t0   <- (n - 1L) * interval_h
      last <- sim[sim$time >= t0 & sim$time < t0 + interval_h, , drop = FALSE]
      if (nrow(last) == 0L) return(NA_real_)
      mean(last$Cc > threshold, na.rm = TRUE) * 100
    }

    # Find minimum dose (250 mg steps) within daily limit of 6 g/day.
    # Uses per-model initial conditions from last_dose_r (Option B).
    dose_search <- function(nm, interval_h, dur_h, daily_limit = 6000, step = 250) {
      p <- params_per_model[[nm]]
      if (is.null(p)) return(list(dose = NA_real_, max_pct = NA_real_))
      ppk        <- model_list[[nm]]$ppk_model
      inits      <- inits_list[[nm]]   # NULL is fine — rxSolve uses model defaults
      max_single <- floor(daily_limit / (24 / interval_h) / step) * step
      doses      <- seq(step, max(step, max_single), by = step)
      best_pct   <- 0
      for (d in doses) {
        pct <- sim_tmic_ss(ppk, p, d, interval_h, dur_h, inits = inits)
        if (!is.na(pct)) best_pct <- max(best_pct, pct)
        if (!is.na(pct) && pct >= target_pct)
          return(list(dose = d, max_pct = round(pct, 1)))
      }
      list(dose = NA_real_, max_pct = round(best_pct, 1))
    }

    # Model weights from metrics_summary (NA for failed models)
    w_named <- setNames(
      metrics_summary$weight[match(names(model_list), metrics_summary$model_name)],
      names(model_list)
    )

    # Infusion strategies
    strategies <- list(
      list(label = "Intermittent (0.5 h)", dur = 0.5, intervals = c(8, 12, 24)),
      list(label = "Extended 3 h",         dur = 3,   intervals = c(8, 12, 24)),
      list(label = "Extended 4 h",         dur = 4,   intervals = c(8, 12, 24)),
      list(label = "Continuous",           dur = NA,  intervals = 24)
    )

    rec_rows <- list()
    for (strat in strategies) {
      for (tau in strat$intervals) {
        dur_h_rec <- if (is.na(strat$dur)) tau else strat$dur
        if (!is.na(strat$dur) && strat$dur >= tau) next   # e.g. 4 h infusion q6h is fine, but skip impossible

        # Best model
        best_rec <- dose_search(best_model, tau, dur_h_rec)

        # Model average: weighted mean of per-model minimum doses
        per_m_doses <- sapply(names(model_list), function(nm)
          dose_search(nm, tau, dur_h_rec)$dose)
        valid_m <- !is.na(per_m_doses)
        if (any(valid_m)) {
          w_v       <- w_named[valid_m]; w_v <- w_v / sum(w_v, na.rm = TRUE)
          moa_dose  <- round(sum(per_m_doses[valid_m] * w_v, na.rm = TRUE) / 250) * 250
          moa_pct   <- sim_tmic_ss(model_list[[best_model]]$ppk_model,
                                    params_per_model[[best_model]],
                                    moa_dose, tau, dur_h_rec,
                                    inits = inits_list[[best_model]])
          moa_rec   <- list(dose = moa_dose, max_pct = round(moa_pct, 1))
        } else {
          moa_rec <- list(dose = NA_real_, max_pct = NA_real_)
        }

        fmt_rec <- function(rec, tau) {
          if (!is.na(rec$dose)) {
            daily <- rec$dose * (24 / tau)
            if (tau == 24)
              sprintf("%g mg/day (%.1f%%)", daily, rec$max_pct)
            else
              sprintf("%g mg (%.1f%%)", rec$dose, rec$max_pct)
          } else {
            if (!is.na(rec$max_pct) && rec$max_pct > 0)
              sprintf("Not achievable (max %.1f%%)", rec$max_pct)
            else
              "Not achievable"
          }
        }

        rec_rows[[length(rec_rows) + 1]] <- tibble(
          Infusion        = strat$label,
          Interval        = paste0("q", tau, "h"),
          `Best model`    = fmt_rec(best_rec, tau),
          `Model average` = fmt_rec(moa_rec,  tau),
          best_ok         = as.integer(!is.na(best_rec$dose)),
          moa_ok          = as.integer(!is.na(moa_rec$dose))
        )

        msg("Dose rec [", strat$label, " q", tau, "h] best=",
            if (is.na(best_rec$dose)) "NA" else paste0(best_rec$dose, " mg"),
            " | moa=",
            if (is.na(moa_rec$dose))  "NA" else paste0(moa_rec$dose,  " mg"))
      }
    }
    rec_df <- bind_rows(rec_rows)
    msg("Dose recommendations complete.")

    # ── Render outputs ────────────────────────────────────────────────────────
    msg("Rendering outputs...")
    output$concentrationPlot  <- renderPlotly(combined_plot)
    output$concentrationPlot2 <- renderPlotly(plot_best)
    output$tmic_plot          <- renderPlotly(plot_tmic)

    output$message <- renderUI(HTML(paste0(
      "MAP estimation completed using <strong>posologyr</strong>. ",
      "Best model: <strong>", best_model, "</strong> (highest likelihood weight). ",
      "The %fT>xMIC bar chart shows per-interval target attainment ",
      "(blue = target met, red = target missed)."
    )))

    output$box_tmic_summary <- renderUI({
      if (is.null(last_tmic)) return(NULL)
      colour <- if (last_tmic$pct_tmic >= target_pct) "#2e7d32" else "#c62828"
      fluidRow(
        column(width = 6, offset = 3,
          box(
            title  = tagList(icon("clock"),
                             HTML("<strong style='font-size:18px;'>Last interval %fT>xMIC</strong>")),
            status = "primary", solidHeader = TRUE, width = 12,
            div(style = paste0("text-align:center; padding:20px; color:", colour, ";"),
              h3(paste0(round(last_tmic$pct_tmic, 1), "%")),
              p(paste0("Target: ", target_pct, "% | Threshold: ",
                       round(threshold, 2), " mg/L (free)"))
            )
          )
        )
      )
    })

    output$box_parameters <- renderUI({
      box(
        title = tagList(icon("notes-medical"), " Covariates used"),
        status = "primary", solidHeader = TRUE, width = 12,
        div(style = "padding:20px;", DTOutput("combinedTable"))
      )
    })

    output$combinedTable <- renderDT({
      datatable(patient_tb, options = list(dom = "t"), rownames = FALSE)
    })

    output$box_dose_rec <- renderUI({
      box(
        title = tagList(icon("calculator"),
                        HTML("<strong style='font-size:16px;'> Minimum dose to achieve target</strong>")),
        status = "primary", solidHeader = TRUE, width = 12,
        p(paste0("Minimum dose per administration achieving ", target_pct,
                 "% fT>", mic_mult, "×MIC (",
                 round(threshold, 2), " mg/L free) at steady state. ",
                 "Daily dose limit: 6 g/day. ",
                 "Model average = weighted mean of per-model minimum doses.")),
        DTOutput("dose_rec_table")
      )
    })

    output$dose_rec_table <- renderDT({
      disp <- rec_df %>% select(Infusion, Interval, `Best model`, `Model average`,
                                 best_ok, moa_ok)
      datatable(
        disp,
        options = list(
          dom = "t", pageLength = 20,
          columnDefs = list(
            list(visible = FALSE, targets = c(4, 5))   # hide best_ok, moa_ok
          )
        ),
        rownames = FALSE
      ) %>%
        formatStyle("Best model",    "best_ok",
                    backgroundColor = styleEqual(c(1, 0), c("#d4edda", "#f8d7da"))) %>%
        formatStyle("Model average", "moa_ok",
                    backgroundColor = styleEqual(c(1, 0), c("#d4edda", "#f8d7da")))
    })

    msg("=== core_function COMPLETE ===")
  }

  # ── Save / Load ──────────────────────────────────────────────────────────────
  output$save_state <- downloadHandler(
    filename = function() paste0("mero_state_", Sys.Date(), ".zip"),
    content  = function(file) {
      state <- list(
        weights     = rv$weights,
        creatinines = rv$creatinines,
        albumins    = rv$albumins,
        tdm         = rv$tdm,
        regimens    = rv$regimens,
        target_mic  = input$target_mic,
        mic_mult    = input$mic_multiplier,
        target_pct  = input$target_pct,
        prot_bind   = input$protein_binding
      )
      tmp <- tempfile(fileext = ".rds")
      saveRDS(state, tmp)
      zip(file, tmp, flags = "-j")
      showNotification("State saved.", type = "message")
    }
  )

  observeEvent(input$load_state, {
    req(input$load_state)
    tmp_dir <- file.path(tempdir(), "mero_unzip")
    dir.create(tmp_dir, showWarnings = FALSE)
    unzip(input$load_state$datapath, exdir = tmp_dir)
    rds <- list.files(tmp_dir, "\\.rds$", full.names = TRUE)
    if (length(rds) == 1) {
      s <- readRDS(rds)
      rv$weights     <- s$weights
      rv$creatinines <- s$creatinines
      rv$albumins    <- s$albumins
      rv$tdm         <- s$tdm
      rv$regimens    <- s$regimens
      updateNumericInput(session, "target_mic",     value = s$target_mic)
      updateNumericInput(session, "mic_multiplier", value = s$mic_mult)
      updateNumericInput(session, "target_pct",     value = s$target_pct)
      updateNumericInput(session, "protein_binding",value = s$prot_bind)
      showNotification("State loaded.", type = "message")
    }
    unlink(tmp_dir, recursive = TRUE)
  })
}

shinyApp(ui = ui, server = server)
