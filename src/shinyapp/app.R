library(shiny)
library(ggplot2)
library(randomForestSRC)

## ------------- load the fitted RSF once -----------------------------------
rsf_wt <- readRDS("shinyapp/rsf_wt.RDS")   # adjust path if needed

## ------------- user interface ---------------------------------------------
ui <- fluidPage(
  titlePanel("Patient-specific VTE PD Plot"),
  sidebarLayout(
    sidebarPanel(
      # core injury / physiology --------------------------------------------
      numericInput("gcs",  "TBI Highest Total GCS:",    15, min = 3,  max = 15),
      numericInput("iss",  "Injury Severity Score:",    10, min = 0),
      numericInput("age",  "Age (years):",              50, min = 0),
      
      # binary clinical history ---------------------------------------------
      checkboxInput("hxAnticoag", "Hx of Anticoagulant Therapy",            FALSE),
      checkboxInput("abxPreVTE",  "Antibiotic Therapy before VTE",         FALSE),
      checkboxInput("ctrlSurg",   "Hemorrhage-control Surgery before VTE", FALSE),
      checkboxInput("ventDay1",   "Ventilator Assisted (Day 1)",           FALSE),
      checkboxInput("icuEarly",   "ICU Early Intervention",                TRUE),
      
      # transfusion volume (mL within first 4 h) ----------------------------
      numericInput("blood4ml", "Packed RBC (mL, first 4 h):", 0, min = 0),
      numericInput("plate4ml", "Platelets (mL, first 4 h):",  0, min = 0),
      
      # demographics ---------------------------------------------------------
      selectInput("sex", "Sex:", c("Male" = 0, "Female" = 1)),
      
      actionButton("plotPD", "Generate PD Plot")
    ),
    mainPanel(
      plotOutput("pdPlot", height = 400),
      verbatimTextOutput("peakInfo")
    )
  )
)

## ------------- server logic -----------------------------------------------
server <- function(input, output, session) {
  
  ## gather covariates from the UI into a *named list*
  base_covariates <- reactive({
    list(
      TBIHIGHESTTOTALGCS           = input$gcs,
      ISS                          = input$iss,
      AGEyears                     = input$age,
      Hx_AnticoagulantTherapy      = as.numeric(input$hxAnticoag),
      ANTIBIOTICTHERAPY_Before_VTE = as.numeric(input$abxPreVTE),
      HMRRHGCTRLSURG_Before_VTE    = as.numeric(input$ctrlSurg),
      Ventilator_Assisted_Day1     = as.numeric(input$ventDay1),
      ICU_Early_Intervention       = as.numeric(input$icuEarly),
      BLOOD4ML                     = input$blood4ml,
      PLATELETS4ML                 = input$plate4ml,
      SEX                          = as.numeric(input$sex)
    )
  })
  
  ## build newdata grid + predictions each time the button is pressed
  pd_data <- eventReactive(input$plotPD, ignoreNULL = FALSE, {
    req(rsf_wt)                            # fail fast if model missing
    
    vte_seq  <- seq(0, 168, length.out = 300)   # 0 → 7 days
    bcov_df  <- as.data.frame(base_covariates(), stringsAsFactors = FALSE)
    
    # replicate baseline row for every hour, then add the varying column
    newdata  <- bcov_df[rep(1, length(vte_seq)), , drop = FALSE]
    newdata$VTEPROPHYLAXISHRS <- vte_seq
    
    # ensure column order matches training set
    newdata  <- newdata[, names(rsf_wt$xvar), drop = FALSE]
    
    # predict survival to the *last* time-point
    pred     <- predict(rsf_wt, newdata = newdata)
    surv_vec <- pred$survival[, ncol(pred$survival)]
    best_idx <- which.max(surv_vec)
    
    list(seq = vte_seq,
         surv = surv_vec,
         best_idx = best_idx,
         newdata = newdata)
  })
  
  ## --------- plot ---------------------------------------------------------
  output$pdPlot <- renderPlot({
    res <- pd_data()
    df  <- data.frame(VTEHR = res$seq, Surv = res$surv)
    
    best_hr   <- df$VTEHR[res$best_idx]
    best_surv <- df$Surv [res$best_idx]
    
    ggplot(df, aes(VTEHR, Surv)) +
      geom_line(linetype = "dashed", linewidth = 0.8) +
      geom_vline(xintercept = best_hr, linetype = "dotted") +
      geom_point(aes(best_hr, best_surv), size = 3) +
      labs(title = sprintf("Predicted Survival at %.0f hrs",
                           max(rsf_wt$time.interest)),
           x = "Time to first VTE prophylaxis (hrs)",
           y = "Survival probability") +
      theme_minimal(base_size = 13)
  })
  
  ## --------- peak text with 95 % CI --------------------------------------
  output$peakInfo <- renderPrint({
    res <- pd_data()
    try({
      allpred  <- predict(rsf_wt, newdata = res$newdata, predict.all = TRUE)
      surv_all <- if (is.list(allpred$survival))
        simplify2array(allpred$survival)
      else
        allpred$survival
      
      # grab survival at the final time-point for every tree
      tp   <- if (length(dim(surv_all)) == 3) dim(surv_all)[2] else 1
      vals <- if (length(dim(surv_all)) == 3)
        surv_all[res$best_idx, tp, ]
      else
        surv_all[res$best_idx, ]
      
      ci <- quantile(vals, c(0.025, 0.975), na.rm = TRUE)
      
      cat(sprintf("Peak at %.1f hrs  →  Survival = %.3f  (95%% CI %.3f – %.3f)",
                  res$seq[res$best_idx], res$surv[res$best_idx], ci[1], ci[2]))
    }, silent = TRUE)
  })
}

shinyApp(ui, server)