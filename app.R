
##############################################################################
#  ADC DAR Monitoring & PK App
#  app.R  —  R Shiny Dashboard  (v1.0)
#
#  Analytes : Total ADC | Total Antibody | Free Payload
#  PK Method: Non-Compartmental Analysis (PKNCA 0.12.x)
#
#  DAR formula: avg_DAR(t) = [Total_ADC(t) / Total_Ab(t)] × Nominal_DAR
#    Rationale: Total_ADC and Total_Ab are both measured in mass conc (ug/mL).
#    Their ratio gives fractional conjugation; multiplying by nominal DAR
#    converts to average drug-to-antibody ratio in drug molecules per Ab.
##############################################################################

library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinycssloaders)
library(PKNCA)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(openxlsx)
library(scales)

# ── Constants ─────────────────────────────────────────────────────────────────
REQUIRED_COLS <- c("Time_h","Total_ADC_mean","Total_Ab_mean",
                   "FreePayload_mean","Dose_mg_kg")
OPT_SD_COLS   <- c("Total_ADC_SD","Total_Ab_SD","FreePayload_SD")

COL_ADC     <- "#0279EE"
COL_AB      <- "#75A025"
COL_PAYLOAD <- "#FF9400"
COL_DAR_FIT <- "#FD9BED"

# ── Helpers ───────────────────────────────────────────────────────────────────

validate_csv <- function(df) {
  missing <- setdiff(REQUIRED_COLS, names(df))
  if (length(missing) > 0)
    return(paste("Missing required columns:", paste(missing, collapse=", ")))
  if (any(df$Time_h < 0, na.rm=TRUE))
    return("Time_h must be >= 0.")
  if (any(df$Total_ADC_mean < 0, na.rm=TRUE) ||
      any(df$Total_Ab_mean  < 0, na.rm=TRUE) ||
      any(df$FreePayload_mean < 0, na.rm=TRUE))
    return("Concentration values must be >= 0.")
  if (any(df$Total_Ab_mean == 0, na.rm=TRUE))
    return("Total_Ab_mean contains zero values — DAR cannot be computed.")
  NULL
}

# NCA via PKNCA 0.12.x
# - Uses pk.nca() (not pk_nca())
# - No subject= arg in PKNCAconc
# - Interval start = first observed timepoint (not 0) to avoid AUC exclusion
run_nca <- function(df, analyte_col, dose_val) {
  t_start <- min(df$Time_h, na.rm=TRUE)
  t_end   <- max(df$Time_h, na.rm=TRUE)

  conc_data <- PKNCAconc(df, as.formula(paste(analyte_col, "~ Time_h")))
  dose_data <- PKNCAdose(
    data.frame(Time_h = 0, Dose = dose_val),
    Dose ~ Time_h, route = "intravascular"
  )
  pk_data <- PKNCAdata(conc_data, dose_data,
    intervals = data.frame(
      start     = t_start, end = t_end,
      auclast   = TRUE,
      cmax      = TRUE, tmax      = TRUE,
      half.life = TRUE
    )
  )
  res <- pk.nca(pk_data)
  as.data.frame(res$result) %>%
    filter(is.na(exclude) | exclude == "") %>%
    select(PPTESTCD, PPORRES) %>%
    rename(Parameter = PPTESTCD, Value = PPORRES) %>%
    mutate(Value = round(as.numeric(Value), 4))
}

# DAR half-life via exponential decay fit
dar_halflife <- function(time, dar) {
  df <- data.frame(t = time, dar = dar) %>% filter(is.finite(dar) & dar > 0)
  if (nrow(df) < 3) return(list(k=NA, t_half=NA, fitted=rep(NA, length(time))))
  tryCatch({
    fit <- nls(dar ~ A * exp(-k * t), data=df,
               start=list(A=max(df$dar), k=0.005),
               control=nls.control(maxiter=500))
    k_val <- coef(fit)[["k"]]
    list(k      = round(k_val, 6),
         t_half = round(log(2)/k_val, 2),
         fitted = predict(fit))
  }, error=function(e) list(k=NA, t_half=NA, fitted=rep(NA, length(time))))
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "black",

  dashboardHeader(
    title = tags$span("ADC DAR Monitoring-PK Tracker"),
    titleWidth = 320
  ),

  dashboardSidebar(
    width = 280,
    useShinyjs(),
    tags$head(tags$style(HTML("
      .skin-black .main-sidebar { background-color: #1a1a2e; }
      .skin-black .sidebar a    { color: #c8d6e5 !important; }
      .sidebar-menu > li.active > a { border-left: 4px solid #0279EE !important; }
      .box { border-top: 3px solid #0279EE; }
      body { font-family: Arial, Helvetica, sans-serif; }
      .info-box { min-height: 70px; }
      .info-box-icon { height: 70px; line-height: 70px; width: 80px; }
      .info-box-content { padding-top: 8px; }
    "))),

    sidebarMenu(id="sidebar_menu",
      menuItem("Upload & Configure", tabName="upload",  icon=icon("upload")),
      menuItem("PK Curves",          tabName="pk",      icon=icon("chart-line")),
      menuItem("DAR Trajectory",     tabName="dar",     icon=icon("arrow-trend-down")),
      menuItem("Payload Kinetics",   tabName="payload", icon=icon("pills")),
      menuItem("PK/PD Overlay",      tabName="pkpd",    icon=icon("layer-group")),
      menuItem("NCA Parameters",     tabName="nca",     icon=icon("table")),
      menuItem("Export",             tabName="export",  icon=icon("file-export"))
    ),

    hr(),
    div(style="padding:8px 15px; color:#aaa; font-size:11px;",
      strong(style="color:#ddd;","Study Metadata"),
      br(),
      textInput("study_id",    "Study ID",       value="ADC-PK-001"),
      textInput("compound",    "Compound Name",  value="ADC-001"),
      textInput("species",     "Species",        value="Cynomolgus monkey"),
      numericInput("nominal_dar","Nominal DAR (T0)", value=3.5, min=1, max=8, step=0.1),
      selectInput("conc_unit","Concentration Unit",
                  choices=c("ug/mL","ng/mL","nM","pM"), selected="ug/mL"),
      selectInput("yscale","Y-axis Scale",
                  choices=c("Linear"="linear","Log10"="log"), selected="log")
    )
  ),

  dashboardBody(
    tabItems(

      # ── Upload ──────────────────────────────────────────────────────────────
      tabItem(tabName="upload",
        fluidRow(
          box(width=12, title="Data Upload", status="primary", solidHeader=TRUE,
            fluidRow(
              column(6,
                fileInput("csv_file","Upload CSV File", accept=".csv",
                          buttonLabel="Browse...", placeholder="No file selected"),
                tags$small(style="color:#888;",
                  "Required: Time_h, Total_ADC_mean, Total_Ab_mean,
                   FreePayload_mean, Dose_mg_kg"),
                br(),
                tags$small(style="color:#888;",
                  "Optional SD: Total_ADC_SD, Total_Ab_SD, FreePayload_SD"),
                br(), br(),
                downloadButton("dl_template","Download CSV Template",
                               class="btn-info btn-sm")
              ),
              column(6,
                div(style="background:#f8f9fa;padding:15px;border-radius:6px;
                           border-left:4px solid #0279EE;",
                  h5(icon("info-circle")," DAR Calculation Method"),
                  p(style="font-size:12px;",
                    "avg DAR(t) = [Total_ADC(t) / Total_Ab(t)] × Nominal_DAR"),
                  p(style="font-size:11px;color:#666;",
                    "Both analytes measured in mass concentration (ug/mL).
                     Their ratio gives fractional conjugation; multiplying by
                     nominal DAR converts to drug molecules per antibody.")
                )
              )
            )
          )
        ),
        fluidRow(
          box(width=12, title="Data Preview", status="info", solidHeader=TRUE,
            uiOutput("validation_msg"),
            DTOutput("data_preview") %>% withSpinner(color="#0279EE")
          )
        )
      ),

      # ── PK Curves ───────────────────────────────────────────────────────────
      tabItem(tabName="pk",
        fluidRow(
          box(width=9, title="Multi-Analyte PK Concentration-Time Profiles",
              status="primary", solidHeader=TRUE,
            plotlyOutput("pk_plot", height="480px") %>% withSpinner(color="#0279EE")
          ),
          box(width=3, title="Plot Options", status="info", solidHeader=TRUE,
            checkboxGroupInput("analytes_shown","Show Analytes:",
              choices=c("Total ADC"="adc","Total Antibody"="ab","Free Payload"="payload"),
              selected=c("adc","ab","payload")),
            checkboxInput("show_errorbars","Show ±SD error bars", value=TRUE),
            checkboxInput("show_points",   "Show data points",    value=TRUE),
            hr(),
            tableOutput("pk_summary_tbl")
          )
        )
      ),

      # ── DAR Trajectory ──────────────────────────────────────────────────────
      tabItem(tabName="dar",
        fluidRow(
          box(width=8, title="Average DAR Over Time",
              status="primary", solidHeader=TRUE,
            plotlyOutput("dar_plot", height="430px") %>% withSpinner(color="#0279EE")
          ),
          box(width=4, title="DAR Stability Metrics", status="info", solidHeader=TRUE,
            tableOutput("dar_metrics_tbl"),
            hr(),
            div(style="font-size:11px;color:#666;",
              strong("Formula:"), br(),
              "avg DAR(t) = [ADC/Ab] × Nominal_DAR", br(), br(),
              strong("Decay fit:"), br(),
              "DAR(t) = DAR₀ × exp(−k × t)", br(),
              "DAR t½ = ln(2) / k"
            )
          )
        ),
        fluidRow(
          box(width=12, title="DAR Trajectory Table", status="info", solidHeader=TRUE,
            DTOutput("dar_table") %>% withSpinner(color="#0279EE")
          )
        )
      ),

      # ── Payload Kinetics ────────────────────────────────────────────────────
      tabItem(tabName="payload",
        fluidRow(
          box(width=8, title="Free Payload Release Kinetics",
              status="primary", solidHeader=TRUE,
            plotlyOutput("payload_plot", height="430px") %>% withSpinner(color="#0279EE")
          ),
          box(width=4, title="Payload PK Highlights", status="info", solidHeader=TRUE,
            tableOutput("payload_metrics_tbl"),
            hr(),
            div(style="font-size:11px;color:#666;",
              "Payload release rate (ΔC/Δt) estimated as finite difference
               between consecutive timepoints.")
          )
        )
      ),

      # ── PK/PD Overlay ───────────────────────────────────────────────────────
      tabItem(tabName="pkpd",
        fluidRow(
          box(width=12,
              title="PK/PD Overlay: DAR Trajectory + Free Payload (Dual Axis)",
              status="primary", solidHeader=TRUE,
            plotlyOutput("pkpd_plot", height="480px") %>% withSpinner(color="#0279EE")
          )
        ),
        fluidRow(
          box(width=12, title="Biotransformation Summary",
              status="info", solidHeader=TRUE,
            uiOutput("biotr_summary")
          )
        )
      ),

      # ── NCA Parameters ──────────────────────────────────────────────────────
      tabItem(tabName="nca",
        fluidRow(
          box(width=12,
              title="Non-Compartmental Analysis (NCA) — All Analytes",
              status="primary", solidHeader=TRUE,
            div(style="font-size:11px;color:#888;padding-bottom:8px;",
              icon("info-circle"),
              " NCA interval: first observed timepoint to last observed timepoint.
                AUClast computed by linear-log trapezoidal rule (PKNCA default).
                Half-life estimated from terminal log-linear phase."
            ),
            tabsetPanel(
              tabPanel("Total ADC",
                DTOutput("nca_adc") %>% withSpinner(color="#0279EE")),
              tabPanel("Total Antibody",
                DTOutput("nca_ab")  %>% withSpinner(color="#0279EE")),
              tabPanel("Free Payload",
                DTOutput("nca_payload") %>% withSpinner(color="#0279EE")),
              tabPanel("Comparison",
                DTOutput("nca_compare") %>% withSpinner(color="#0279EE"))
            )
          )
        )
      ),

      # ── Export ──────────────────────────────────────────────────────────────
      tabItem(tabName="export",
        fluidRow(
          box(width=6, title="PDF Report", status="primary", solidHeader=TRUE,
            p("Regulatory-ready PDF containing:"),
            tags$ul(
              tags$li("Study metadata header"),
              tags$li("Multi-analyte PK concentration-time plots"),
              tags$li("DAR trajectory with exponential decay fit"),
              tags$li("Free payload release kinetics"),
              tags$li("PK/PD overlay (dual axis)"),
              tags$li("NCA parameter tables for all 3 analytes"),
              tags$li("DAR stability summary & methods section")
            ),
            br(),
            downloadButton("dl_pdf","Download PDF Report",
                           class="btn-danger btn-lg")
          ),
          box(width=6, title="Excel Workbook", status="success", solidHeader=TRUE,
            p("Formatted Excel workbook with sheets:"),
            tags$ul(
              tags$li("Raw Data"),
              tags$li("DAR Trajectory"),
              tags$li("NCA — Total ADC"),
              tags$li("NCA — Total Antibody"),
              tags$li("NCA — Free Payload"),
              tags$li("Biotransformation Summary")
            ),
            br(),
            downloadButton("dl_excel","Download Excel Workbook",
                           class="btn-success btn-lg")
          )
        ),
        fluidRow(
          box(width=12, title="Export Status", status="info", solidHeader=TRUE,
            uiOutput("export_status")
          )
        )
      )
    )
  )
)

# ── SERVER ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # ── Load & validate ─────────────────────────────────────────────────────────
  raw_data <- reactive({
    req(input$csv_file)
    tryCatch(read.csv(input$csv_file$datapath, stringsAsFactors=FALSE),
             error=function(e) NULL)
  })

  validated_data <- reactive({
    df <- raw_data(); req(df)
    err <- validate_csv(df)
    if (!is.null(err)) return(NULL)
    df <- df %>% arrange(Time_h)
    for (col in OPT_SD_COLS) if (!col %in% names(df)) df[[col]] <- NA_real_
    df
  })

  output$validation_msg <- renderUI({
    df <- raw_data(); if (is.null(df)) return(NULL)
    err <- validate_csv(df)
    if (!is.null(err))
      div(class="alert alert-danger", icon("times-circle"), " ", err)
    else
      div(class="alert alert-success", icon("check-circle"),
          sprintf(" Data loaded: %d timepoints | Dose = %s mg/kg | Study: %s",
                  nrow(df),
                  unique(df$Dose_mg_kg)[1],
                  if("Study_ID" %in% names(df)) unique(df$Study_ID)[1] else input$study_id))
  })

  output$data_preview <- renderDT({
    df <- validated_data(); req(df)
    datatable(df, options=list(scrollX=TRUE, pageLength=15), rownames=FALSE) %>%
      formatRound(intersect(names(df),
        c("Total_ADC_mean","Total_ADC_SD","Total_Ab_mean","Total_Ab_SD",
          "FreePayload_mean","FreePayload_SD")), digits=3)
  })

  # ── DAR & derived quantities ─────────────────────────────────────────────────
  dar_data <- reactive({
    df <- validated_data(); req(df)
    df %>% mutate(
      Frac_conj  = Total_ADC_mean / Total_Ab_mean,
      DAR        = Frac_conj * input$nominal_dar,
      DAR_pct    = Frac_conj * 100,   # % of nominal DAR retained
      Payload_release_rate = c(NA, diff(FreePayload_mean) / diff(Time_h))
    )
  })

  dar_fit <- reactive({
    dd <- dar_data(); req(dd)
    dar_halflife(dd$Time_h, dd$DAR)
  })

  # ── NCA ─────────────────────────────────────────────────────────────────────
  nca_results <- reactive({
    df   <- validated_data(); req(df)
    dose <- unique(df$Dose_mg_kg)[1]
    list(
      adc     = tryCatch(run_nca(df,"Total_ADC_mean",   dose),
                         error=function(e) data.frame(Parameter="Error",Value=e$message)),
      ab      = tryCatch(run_nca(df,"Total_Ab_mean",    dose),
                         error=function(e) data.frame(Parameter="Error",Value=e$message)),
      payload = tryCatch(run_nca(df,"FreePayload_mean", dose),
                         error=function(e) data.frame(Parameter="Error",Value=e$message))
    )
  })

  # ── Template download ────────────────────────────────────────────────────────
  output$dl_template <- downloadHandler(
    filename="adc_pk_template.csv",
    content=function(file) {
      write.csv(data.frame(
        Time_h=c(0.083,1,4,24,48,96,168,336,504,672),
        Total_ADC_mean=NA_real_, Total_ADC_SD=NA_real_,
        Total_Ab_mean=NA_real_,  Total_Ab_SD=NA_real_,
        FreePayload_mean=NA_real_, FreePayload_SD=NA_real_,
        Dose_mg_kg=NA_real_, Species="", Study_ID=""
      ), file, row.names=FALSE)
    }
  )

  # ── PK Curves ───────────────────────────────────────────────────────────────
  output$pk_plot <- renderPlotly({
    df <- dar_data(); req(df)
    p  <- plot_ly()
    yt <- if (input$yscale=="log") "log" else "linear"

    add_analyte <- function(p, x, y, sd, name, color, show) {
      if (!show) return(p)
      has_sd <- !all(is.na(sd))
      mode   <- if (input$show_points) "lines+markers" else "lines"
      if (input$show_errorbars && has_sd) {
        add_trace(p, x=x, y=y, type="scatter", mode=mode, name=name,
                  line=list(color=color,width=2),
                  marker=list(color=color,size=7),
                  error_y=list(type="data",array=sd,color=color,
                               thickness=1.5,width=4))
      } else {
        add_trace(p, x=x, y=y, type="scatter", mode=mode, name=name,
                  line=list(color=color,width=2),
                  marker=list(color=color,size=7))
      }
    }

    p <- add_analyte(p,df$Time_h,df$Total_ADC_mean,  df$Total_ADC_SD,
                     "Total ADC",      COL_ADC,     "adc"     %in% input$analytes_shown)
    p <- add_analyte(p,df$Time_h,df$Total_Ab_mean,   df$Total_Ab_SD,
                     "Total Antibody", COL_AB,      "ab"      %in% input$analytes_shown)
    p <- add_analyte(p,df$Time_h,df$FreePayload_mean,df$FreePayload_SD,
                     "Free Payload",   COL_PAYLOAD, "payload" %in% input$analytes_shown)

    p %>% layout(
      title=list(text="PK Concentration-Time Profiles",font=list(family="Arial")),
      xaxis=list(title="Time (h)",gridcolor="#eee",zeroline=FALSE),
      yaxis=list(title=paste0("Concentration (",input$conc_unit,")"),
                 type=yt,gridcolor="#eee"),
      legend=list(orientation="h",y=-0.18),
      hovermode="x unified",
      plot_bgcolor="#fafafa", paper_bgcolor="#ffffff",
      font=list(family="Arial")
    )
  })

  output$pk_summary_tbl <- renderTable({
    df <- validated_data(); req(df)
    data.frame(
      Analyte=c("Total ADC","Total Ab","Free Payload"),
      Cmax   =round(c(max(df$Total_ADC_mean,na.rm=TRUE),
                      max(df$Total_Ab_mean,na.rm=TRUE),
                      max(df$FreePayload_mean,na.rm=TRUE)),3),
      Tmax_h =c(df$Time_h[which.max(df$Total_ADC_mean)],
                df$Time_h[which.max(df$Total_Ab_mean)],
                df$Time_h[which.max(df$FreePayload_mean)])
    )
  }, striped=TRUE, hover=TRUE, bordered=TRUE, digits=3)

  # ── DAR Trajectory ──────────────────────────────────────────────────────────
  output$dar_plot <- renderPlotly({
    dd  <- dar_data();  req(dd)
    fit <- dar_fit()

    p <- plot_ly() %>%
      add_trace(x=dd$Time_h, y=dd$DAR, type="scatter", mode="lines+markers",
                name="Observed avg DAR",
                line=list(color=COL_ADC,width=2.5),
                marker=list(color=COL_ADC,size=9)) %>%
      add_trace(x=range(dd$Time_h),
                y=c(input$nominal_dar,input$nominal_dar),
                type="scatter", mode="lines", name="Nominal DAR",
                line=list(color="#aaa",width=1.5,dash="dot"))

    if (!is.na(fit$t_half)) {
      fitted_dar <- dd$DAR[1] * exp(-fit$k * dd$Time_h)
      p <- p %>%
        add_trace(x=dd$Time_h, y=fitted_dar, type="scatter", mode="lines",
                  name=paste0("Exp. fit (t½=",fit$t_half,"h)"),
                  line=list(color=COL_DAR_FIT,width=2,dash="dash"))
    }

    p %>% layout(
      title=list(text="Average DAR Over Time",font=list(family="Arial")),
      xaxis=list(title="Time (h)",gridcolor="#eee"),
      yaxis=list(title="Average DAR",range=c(0,input$nominal_dar*1.25),
                 gridcolor="#eee"),
      legend=list(orientation="h",y=-0.18),
      hovermode="x unified",
      plot_bgcolor="#fafafa", paper_bgcolor="#ffffff",
      font=list(family="Arial")
    )
  })

  output$dar_metrics_tbl <- renderTable({
    dd  <- dar_data();  req(dd)
    fit <- dar_fit()
    data.frame(
      Metric=c("Nominal DAR","DAR at T-first","DAR at T-last",
               "% Nominal DAR retained","Decay rate k (1/h)","DAR Half-life (h)"),
      Value =c(input$nominal_dar,
               round(dd$DAR[1],3),
               round(tail(dd$DAR,1),3),
               round(tail(dd$DAR_pct,1),1),
               ifelse(is.na(fit$k),"N/A",as.character(round(fit$k,6))),
               ifelse(is.na(fit$t_half),"N/A",as.character(round(fit$t_half,1))))
    )
  }, striped=TRUE, hover=TRUE, bordered=TRUE)

  output$dar_table <- renderDT({
    dd <- dar_data(); req(dd)
    tbl <- dd %>%
      select(Time_h,Total_ADC_mean,Total_Ab_mean,Frac_conj,DAR,DAR_pct) %>%
      mutate(across(where(is.numeric),~round(.x,3))) %>%
      rename(`Time (h)`=Time_h,
             `Total ADC`=Total_ADC_mean,
             `Total Ab`=Total_Ab_mean,
             `Fractional Conjugation`=Frac_conj,
             `Avg DAR`=DAR,
             `% Nominal DAR`=DAR_pct)
    datatable(tbl, rownames=FALSE,
              options=list(pageLength=15,scrollX=TRUE)) %>%
      formatStyle("Avg DAR",
        background=styleColorBar(c(0,max(dd$DAR,na.rm=TRUE)),"#0279EE"),
        backgroundSize="100% 90%",
        backgroundRepeat="no-repeat",
        backgroundPosition="center")
  })

  # ── Payload Kinetics ────────────────────────────────────────────────────────
  output$payload_plot <- renderPlotly({
    dd <- dar_data(); req(dd)
    cmax_idx <- which.max(dd$FreePayload_mean)

    p <- plot_ly() %>%
      add_trace(x=dd$Time_h, y=dd$FreePayload_mean,
                type="scatter", mode="lines+markers",
                name="Free Payload",
                line=list(color=COL_PAYLOAD,width=2.5),
                marker=list(color=COL_PAYLOAD,size=9))

    if (!all(is.na(dd$FreePayload_SD))) {
      p <- p %>%
        add_trace(x=dd$Time_h, y=dd$FreePayload_mean,
                  type="scatter", mode="none", name="±SD",
                  error_y=list(type="data",array=dd$FreePayload_SD,
                               color=COL_PAYLOAD,thickness=1.5,width=4),
                  showlegend=FALSE)
    }

    p <- p %>%
      add_annotations(
        x=dd$Time_h[cmax_idx], y=dd$FreePayload_mean[cmax_idx],
        text=paste0("<b>Cmax = ",round(dd$FreePayload_mean[cmax_idx],3),
                    " ",input$conc_unit,"</b><br>Tmax = ",
                    dd$Time_h[cmax_idx]," h"),
        showarrow=TRUE, arrowhead=2, ax=50, ay=-40,
        font=list(size=11,color=COL_PAYLOAD)
      )

    p %>% layout(
      title=list(text="Free Payload Release Kinetics",font=list(family="Arial")),
      xaxis=list(title="Time (h)",gridcolor="#eee"),
      yaxis=list(title=paste0("Free Payload (",input$conc_unit,")"),
                 gridcolor="#eee"),
      hovermode="x unified",
      plot_bgcolor="#fafafa", paper_bgcolor="#ffffff",
      font=list(family="Arial")
    )
  })

  output$payload_metrics_tbl <- renderTable({
    dd <- dar_data(); req(dd)
    cmax_idx <- which.max(dd$FreePayload_mean)
    auc_trap <- sum(diff(dd$Time_h) *
                    (head(dd$FreePayload_mean,-1)+tail(dd$FreePayload_mean,-1))/2)
    data.frame(
      Metric=c("Cmax","Tmax (h)","Clast","AUClast (trap.)"),
      Value =round(c(dd$FreePayload_mean[cmax_idx],
                     dd$Time_h[cmax_idx],
                     tail(dd$FreePayload_mean,1),
                     auc_trap),4)
    )
  }, striped=TRUE, hover=TRUE, bordered=TRUE)

  # ── PK/PD Overlay ───────────────────────────────────────────────────────────
  output$pkpd_plot <- renderPlotly({
    dd  <- dar_data();  req(dd)
    fit <- dar_fit()

    p <- plot_ly() %>%
      add_trace(x=dd$Time_h, y=dd$DAR,
                type="scatter", mode="lines+markers",
                name="Avg DAR (left)",
                yaxis="y",
                line=list(color=COL_ADC,width=2.5),
                marker=list(color=COL_ADC,size=8)) %>%
      add_trace(x=dd$Time_h, y=dd$FreePayload_mean,
                type="scatter", mode="lines+markers",
                name="Free Payload (right)",
                yaxis="y2",
                line=list(color=COL_PAYLOAD,width=2.5,dash="dash"),
                marker=list(color=COL_PAYLOAD,size=8,symbol="diamond"))

    if (!is.na(fit$t_half)) {
      fitted_dar <- dd$DAR[1] * exp(-fit$k * dd$Time_h)
      p <- p %>%
        add_trace(x=dd$Time_h, y=fitted_dar,
                  type="scatter", mode="lines",
                  name=paste0("DAR fit (t½=",fit$t_half,"h)"),
                  yaxis="y",
                  line=list(color=COL_DAR_FIT,width=1.5,dash="dot"))
    }

    p %>% layout(
      title=list(text="PK/PD Overlay: DAR Trajectory & Free Payload",
                 font=list(family="Arial")),
      xaxis=list(title="Time (h)",gridcolor="#eee"),
      yaxis=list(title="Average DAR",gridcolor="#eee",
                 range=c(0,input$nominal_dar*1.25)),
      yaxis2=list(overlaying="y",side="right",
                  title=paste0("Free Payload (",input$conc_unit,")"),
                  showgrid=FALSE,zeroline=FALSE),
      legend=list(orientation="h",y=-0.18),
      hovermode="x unified",
      plot_bgcolor="#fafafa", paper_bgcolor="#ffffff",
      font=list(family="Arial")
    )
  })

  output$biotr_summary <- renderUI({
    dd  <- dar_data();  req(dd)
    fit <- dar_fit()
    fluidRow(
      column(4,
        div(class="info-box",
          span(class="info-box-icon bg-blue", icon("dna")),
          div(class="info-box-content",
            span(class="info-box-text","DAR Half-life"),
            span(class="info-box-number",
                 ifelse(is.na(fit$t_half),"N/A",paste(fit$t_half,"h")))
          )
        )
      ),
      column(4,
        div(class="info-box",
          span(class="info-box-icon bg-orange", icon("pills")),
          div(class="info-box-content",
            span(class="info-box-text","Payload Cmax"),
            span(class="info-box-number",
                 paste(round(max(dd$FreePayload_mean,na.rm=TRUE),3),input$conc_unit))
          )
        )
      ),
      column(4,
        div(class="info-box",
          span(class="info-box-icon bg-green", icon("percent")),
          div(class="info-box-content",
            span(class="info-box-text","% DAR Retained (Tlast)"),
            span(class="info-box-number",
                 paste0(round(tail(dd$DAR_pct,1),1),"%"))
          )
        )
      )
    )
  })

  # ── NCA Tables ──────────────────────────────────────────────────────────────
  make_nca_dt <- function(df_r) {
    renderDT({
      df <- df_r(); req(df)
      datatable(df, rownames=FALSE,
                options=list(pageLength=20,dom="t",scrollX=TRUE)) %>%
        formatRound("Value", digits=4)
    })
  }
  output$nca_adc     <- make_nca_dt(reactive(nca_results()$adc))
  output$nca_ab      <- make_nca_dt(reactive(nca_results()$ab))
  output$nca_payload <- make_nca_dt(reactive(nca_results()$payload))

  output$nca_compare <- renderDT({
    nca <- nca_results(); req(nca)
    merged <- nca$adc %>% rename(`Total ADC`=Value) %>%
      full_join(nca$ab      %>% rename(`Total Antibody`=Value), by="Parameter") %>%
      full_join(nca$payload %>% rename(`Free Payload`=Value),   by="Parameter")
    datatable(merged, rownames=FALSE,
              options=list(pageLength=20,dom="t",scrollX=TRUE)) %>%
      formatRound(c("Total ADC","Total Antibody","Free Payload"), digits=4)
  })

  # ── Excel Export ─────────────────────────────────────────────────────────────
  output$dl_excel <- downloadHandler(
    filename=function() paste0(input$study_id,"_ADC_PK_",Sys.Date(),".xlsx"),
    content=function(file) {
      dd  <- dar_data();  req(dd)
      nca <- nca_results(); req(nca)
      fit <- dar_fit()
      wb  <- createWorkbook()

      hs <- createStyle(fontColour="#FFFFFF", fgFill="#1a1a2e",
                        halign="CENTER", textDecoration="Bold",
                        border="Bottom", borderColour="#0279EE",
                        fontSize=11, fontName="Arial")
      row_style <- createStyle(fontName="Arial", fontSize=10,
                               border="TopBottomLeftRight",
                               borderColour="#dddddd")

      add_sheet <- function(wb, name, df) {
        addWorksheet(wb, name)
        writeData(wb, name, df, headerStyle=hs)
        if (nrow(df)>0)
          addStyle(wb, name, row_style,
                   rows=2:(nrow(df)+1), cols=1:ncol(df), gridExpand=TRUE)
        setColWidths(wb, name, cols=1:ncol(df), widths="auto")
      }

      add_sheet(wb, "Raw Data", validated_data())

      dar_tbl <- dd %>%
        select(Time_h,Total_ADC_mean,Total_Ab_mean,Frac_conj,DAR,DAR_pct,
               FreePayload_mean,Payload_release_rate) %>%
        rename(`Time (h)`=Time_h,`Total ADC`=Total_ADC_mean,
               `Total Antibody`=Total_Ab_mean,
               `Fractional Conjugation`=Frac_conj,
               `Avg DAR`=DAR,`% Nominal DAR`=DAR_pct,
               `Free Payload`=FreePayload_mean,
               `Payload Release Rate (dC/dt)`=Payload_release_rate)
      add_sheet(wb, "DAR Trajectory", dar_tbl)
      add_sheet(wb, "NCA - Total ADC",      nca$adc)
      add_sheet(wb, "NCA - Total Antibody", nca$ab)
      add_sheet(wb, "NCA - Free Payload",   nca$payload)

      summary_df <- data.frame(
        Parameter=c("Study ID","Compound","Species","Dose (mg/kg)",
                    "Nominal DAR","DAR at T-first","DAR at T-last",
                    "% Nominal DAR Retained","DAR Decay Rate k (1/h)",
                    "DAR Half-life (h)","Payload Cmax","Payload Tmax (h)",
                    "Analysis Date","Software","NCA Method"),
        Value=c(input$study_id, input$compound, input$species,
                unique(validated_data()$Dose_mg_kg)[1],
                input$nominal_dar,
                round(dd$DAR[1],3),
                round(tail(dd$DAR,1),3),
                round(tail(dd$DAR_pct,1),1),
                ifelse(is.na(fit$k),"N/A",round(fit$k,6)),
                ifelse(is.na(fit$t_half),"N/A",round(fit$t_half,1)),
                round(max(dd$FreePayload_mean,na.rm=TRUE),4),
                dd$Time_h[which.max(dd$FreePayload_mean)],
                as.character(Sys.Date()),
                paste0("R ",R.version$major,".",R.version$minor,
                       " | PKNCA ",packageVersion("PKNCA")),
                "NCA (linear-log trapezoidal, PKNCA)")
      )
      add_sheet(wb, "Biotransformation Summary", summary_df)
      saveWorkbook(wb, file, overwrite=TRUE)
    }
  )

  # ── PDF Export ───────────────────────────────────────────────────────────────
  output$dl_pdf <- downloadHandler(
    filename=function() paste0(input$study_id,"_ADC_PK_Report_",Sys.Date(),".pdf"),
    content=function(file) {
      dd  <- dar_data();  req(dd)
      nca <- nca_results(); req(nca)
      fit <- dar_fit()

      params_list <- list(
        study_id    = input$study_id,
        compound    = input$compound,
        species     = input$species,
        dose        = unique(validated_data()$Dose_mg_kg)[1],
        nominal_dar = input$nominal_dar,
        conc_unit   = input$conc_unit,
        dar_data    = dd,
        nca_adc     = nca$adc,
        nca_ab      = nca$ab,
        nca_payload = nca$payload,
        dar_fit_k   = fit$k,
        dar_fit_t12 = fit$t_half,
        col_adc     = COL_ADC,
        col_ab      = COL_AB,
        col_payload = COL_PAYLOAD,
        col_dar_fit = COL_DAR_FIT
      )

      rmd_path <- file.path(getwd(), "report_template.Rmd")
      tmp_out  <- tempfile(fileext=".pdf")
      rmarkdown::render(input=rmd_path, output_file=tmp_out,
                        params=params_list,
                        envir=new.env(parent=globalenv()), quiet=TRUE)
      file.copy(tmp_out, file)
    }
  )

  output$export_status <- renderUI({
    df <- validated_data()
    if (is.null(df))
      div(class="alert alert-warning", icon("exclamation-triangle"),
          " Please upload a valid CSV file before exporting.")
    else
      div(class="alert alert-success", icon("check-circle"),
          " Data loaded. PDF and Excel exports are ready.")
  })
}

shinyApp(ui=ui, server=server)

