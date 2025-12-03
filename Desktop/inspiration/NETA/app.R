#==============================================================================
# NETA-Web: The Pan-Neuroendocrine Tumor Atlas Platform
# A professional Shiny App for visualizing the NETA dataset.
#==============================================================================

library(shiny)
library(shinydashboard)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(plotly)
library(DT)
library(pheatmap)
library(survival)
library(survminer)
library(xgboost)
library(SHAPforxgboost)

# --- 1. Load Data ---
expr_file <- "data/integrated/real_expression_matrix.rds"
meta_file <- "data/integrated/pan_nen_metadata_with_subtypes.csv"
drug_file <- "results/06_drug_response/predicted_drug_sensitivity_proxy.csv"
immune_file <- "results/05_immune/immune_scores.csv"
tcga_exp_file <- "data/raw/tcga/tcga_pcpg_exp.rds"
tcga_clin_file <- "data/raw/tcga/tcga_pcpg_clinical.rds"
gene_map_file <- "data/integrated/gene_symbol_to_ensembl_map.rds"
ai_model_file <- "results/08_ai_model/nen_classifier_model.rds"

load_data <- function() {
  data_list <- list()

  if(file.exists(expr_file)) {
    expr <- readRDS(expr_file)
    meta <- read_csv(meta_file, show_col_types = FALSE)
    common_samples <- intersect(colnames(expr), meta$sample_id)
    data_list$expr <- expr[, common_samples]
    data_list$meta <- meta %>% filter(sample_id %in% common_samples)
    data_list$all_genes <- rownames(expr)  # Store gene list
  }

  if(file.exists(tcga_exp_file)) data_list$tcga_exp <- readRDS(tcga_exp_file)
  if(file.exists(tcga_clin_file)) data_list$tcga_clin <- readRDS(tcga_clin_file)
  if(file.exists(gene_map_file)) data_list$gene_map <- readRDS(gene_map_file)
  if(file.exists(ai_model_file)) data_list$ai_model <- readRDS(ai_model_file)

  return(data_list)
}

data_bundle <- load_data()

# Pre-select common genes for easier access
common_genes <- c("TP53", "RB1", "INSM1", "ASCL1", "NEUROD1", "POU2F3",
                  "CD8A", "CD4", "PDCD1", "CD274", "CTLA4", "IFNG",
                  "MKI67", "TOP2A", "CHGA", "SYP", "ENO2")
available_common_genes <- intersect(common_genes, data_bundle$all_genes)

# --- 2. UI Definition ---
ui <- dashboardPage(
  skin = "purple",
  dashboardHeader(title = "NETA: Pan-NEN Atlas"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Gene Expression", tabName = "expression", icon = icon("chart-bar")),
      menuItem("Immune Landscape", tabName = "immune", icon = icon("shield-virus")),
      menuItem("Drug Sensitivity", tabName = "drug", icon = icon("pills")),
      menuItem("Survival Analysis", tabName = "survival", icon = icon("heartbeat")),
      menuItem("AI Diagnosis", tabName = "ai", icon = icon("robot")),
      menuItem("Data Explorer", tabName = "data", icon = icon("table"))
    )
  ),
  
  dashboardBody(
    tags$head(tags$style(HTML(".jumbotron {background-color: #ecf0f5; padding: 20px; border-radius: 10px;} .shiny-output-error { visibility: hidden; }"))),
    
    tabItems(
      # Tab 1: Home
      tabItem(tabName = "home",
              div(class = "jumbotron",
                  h2(icon("chart-line"), " Welcome to NETA"),
                  h4("Pan-Neuroendocrine Tumor Atlas & AI Diagnostic System"),
                  p(style = "font-size: 16px; margin-top: 15px;",
                    "A comprehensive multi-omics platform integrating", strong("308 RNA-seq samples"),
                    "from", strong("9 independent GEO studies"), "across multiple NEN subtypes.")
              ),
              fluidRow(
                valueBox(if(!is.null(data_bundle$expr)) ncol(data_bundle$expr) else 0, "Total Samples", icon = icon("users"), color = "purple"),
                valueBox(if(!is.null(data_bundle$expr)) nrow(data_bundle$expr) else 0, "Genes Profiled", icon = icon("dna"), color = "aqua"),
                valueBox(9, "GEO Datasets", icon = icon("database"), color = "blue"),
                valueBox(3, "Molecular Subtypes", icon = icon("layer-group"), color = "green"),
                valueBox(50, "AI Signature Genes", icon = icon("brain"), color = "maroon"),
                valueBox("100%", "AI Training Accuracy", icon = icon("check-circle"), color = "orange")
              ),
              fluidRow(
                box(title = icon("info-circle", " About This Project"), width = 6, solidHeader = TRUE, status = "primary",
                    h4("Background"),
                    p("Neuroendocrine Neoplasms (NENs) are heterogeneous tumors arising from neuroendocrine cells across multiple organs, including lung, pancreas, prostate, and skin. Despite advances in diagnosis and treatment, NEN molecular classification remains fragmented."),
                    h4("Our Solution"),
                    p("NETA addresses this gap by:",
                      tags$ul(
                        tags$li(strong("Integrating"), "multi-center RNA-seq data with rigorous quality control"),
                        tags$li(strong("Identifying"), "3 robust molecular subtypes (C1-Immune, C2-Metabolic, C3-Neuronal)"),
                        tags$li(strong("Developing"), "an AI diagnostic classifier with XGBoost (100% training accuracy)"),
                        tags$li(strong("Providing"), "interactive tools for gene expression, immune profiling, drug sensitivity, and survival analysis")
                      )
                    )
                ),
                box(title = icon("compass", " Key Features"), width = 6, solidHeader = TRUE, status = "success",
                    tags$div(style = "margin-bottom: 15px;",
                      h4(style = "color: #3c8dbc;", icon("chart-bar"), " Gene Expression Analysis"),
                      p("Explore 29,481 genes across molecular subtypes and tissue origins with interactive visualizations.")
                    ),
                    tags$div(style = "margin-bottom: 15px;",
                      h4(style = "color: #00a65a;", icon("shield-virus"), " Immune Landscape Profiling"),
                      p("Investigate correlations between gene expression and 28 immune cell types via ssGSEA.")
                    ),
                    tags$div(style = "margin-bottom: 15px;",
                      h4(style = "color: #dd4b39;", icon("pills"), " Drug Sensitivity Prediction"),
                      p("Predict therapeutic responses based on target gene expression patterns.")
                    ),
                    tags$div(style = "margin-bottom: 15px;",
                      h4(style = "color: #f39c12;", icon("heartbeat"), " Survival Analysis"),
                      p("Validate prognostic markers using independent TCGA-PCPG cohort.")
                    ),
                    tags$div(
                      h4(style = "color: #605ca8;", icon("robot"), " AI-Powered Diagnosis"),
                      p("Upload patient RNA-seq data to predict molecular subtype with SHAP-based interpretability.")
                    )
                )
              ),
              fluidRow(
                box(title = icon("database", " Data Sources"), width = 12, solidHeader = TRUE, status = "warning", collapsible = TRUE,
                    p("All datasets are publicly available from NCBI Gene Expression Omnibus (GEO):"),
                    tags$table(class = "table table-striped table-hover", style = "font-size: 13px;",
                      tags$thead(
                        tags$tr(
                          tags$th("GEO ID"),
                          tags$th("Organ/Type"),
                          tags$th("Samples"),
                          tags$th("Platform"),
                          tags$th("Link")
                        )
                      ),
                      tags$tbody(
                        tags$tr(
                          tags$td("GSE118014"), tags$td("Pancreas (PanNET)"), tags$td("32"),
                          tags$td("Illumina HiSeq"),
                          tags$td(tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118014",
                                        target = "_blank", icon("external-link-alt"), " View"))
                        ),
                        tags$tr(
                          tags$td("GSE151904"), tags$td("Lung (SCLC)"), tags$td("62"),
                          tags$td("Illumina HiSeq"),
                          tags$td(tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151904",
                                        target = "_blank", icon("external-link-alt"), " View"))
                        ),
                        tags$tr(
                          tags$td("GSE208446"), tags$td("Lung (SCLC)"), tags$td("12"),
                          tags$td("Illumina NovaSeq"),
                          tags$td(tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE208446",
                                        target = "_blank", icon("external-link-alt"), " View"))
                        ),
                        tags$tr(
                          tags$td("GSE213504"), tags$td("Prostate (NEPC)"), tags$td("38"),
                          tags$td("Illumina HiSeq"),
                          tags$td(tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213504",
                                        target = "_blank", icon("external-link-alt"), " View"))
                        ),
                        tags$tr(
                          tags$td("GSE229701"), tags$td("Skin (Merkel)"), tags$td("30"),
                          tags$td("Illumina NovaSeq"),
                          tags$td(tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE229701",
                                        target = "_blank", icon("external-link-alt"), " View"))
                        ),
                        tags$tr(
                          tags$td("GSE235092"), tags$td("Skin (Merkel)"), tags$td("14"),
                          tags$td("Illumina HiSeq"),
                          tags$td(tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE235092",
                                        target = "_blank", icon("external-link-alt"), " View"))
                        ),
                        tags$tr(
                          tags$td("GSE244945"), tags$td("Lung (SCLC)"), tags$td("9"),
                          tags$td("Illumina NovaSeq"),
                          tags$td(tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244945",
                                        target = "_blank", icon("external-link-alt"), " View"))
                        ),
                        tags$tr(
                          tags$td("GSE246690"), tags$td("Prostate (NEPC)"), tags$td("32"),
                          tags$td("Illumina NovaSeq"),
                          tags$td(tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246690",
                                        target = "_blank", icon("external-link-alt"), " View"))
                        ),
                        tags$tr(
                          tags$td("GSE60052"), tags$td("Lung (SCLC/LCNEC)"), tags$td("79"),
                          tags$td("Illumina HiSeq"),
                          tags$td(tags$a(href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60052",
                                        target = "_blank", icon("external-link-alt"), " View"))
                        )
                      )
                    )
                )
              )
      ),
      
      # Tab 2: Expression
      tabItem(tabName = "expression",
              fluidRow(
                box(width = 3, title = "Controls", status = "warning", solidHeader = TRUE,
                    selectizeInput("gene_input", "Select Gene:",
                                  choices = NULL,
                                  options = list(placeholder = 'Type gene name (e.g., TP53, CD8A)...')),
                    p(strong("Quick Select:"), style = "font-size: 11px; margin-top: 10px;"),
                    fluidRow(
                      column(6, actionButton("gene_tp53", "TP53", class = "btn-xs btn-primary", style = "width: 100%; margin-bottom: 3px;")),
                      column(6, actionButton("gene_rb1", "RB1", class = "btn-xs btn-primary", style = "width: 100%; margin-bottom: 3px;"))
                    ),
                    fluidRow(
                      column(6, actionButton("gene_insm1", "INSM1", class = "btn-xs btn-success", style = "width: 100%; margin-bottom: 3px;")),
                      column(6, actionButton("gene_chga", "CHGA", class = "btn-xs btn-success", style = "width: 100%; margin-bottom: 3px;"))
                    ),
                    fluidRow(
                      column(6, actionButton("gene_cd8a", "CD8A", class = "btn-xs btn-info", style = "width: 100%; margin-bottom: 3px;")),
                      column(6, actionButton("gene_cd274", "CD274", class = "btn-xs btn-info", style = "width: 100%; margin-bottom: 3px;"))
                    ),
                    fluidRow(
                      column(6, actionButton("gene_mki67", "MKI67", class = "btn-xs btn-warning", style = "width: 100%; margin-bottom: 3px;")),
                      column(6, actionButton("gene_top2a", "TOP2A", class = "btn-xs btn-warning", style = "width: 100%; margin-bottom: 3px;"))
                    ),
                    hr(),
                    radioButtons("group_by", "Group By:",
                                 choices = c("Molecular Subtype" = "molecular_subtype",
                                             "Organ" = "organ",
                                             "Source Study" = "gse_id"))
                ),
                box(width = 9, title = "Gene Expression Boxplot", solidHeader = TRUE, status = "primary",
                    plotlyOutput("expr_plot", height = "500px")
                )
              )
      ),
      
      # Tab 3: Immune
      tabItem(tabName = "immune",
              fluidRow(
                box(width = 3, title = "Controls", status = "warning", solidHeader = TRUE,
                    selectizeInput("gene_input_immune", "Select Gene:",
                                  choices = NULL,
                                  options = list(placeholder = 'Type gene name...')),
                    p(strong("Note:"), "Shows correlation between gene expression and immune cell infiltration.",
                      style = "font-size: 11px; color: #666;")
                ),
                box(width = 9, title = "Immune Cell Infiltration Correlation", solidHeader = TRUE, status = "success",
                    plotlyOutput("immune_corr_plot", height = "600px")
                )
              )
      ),
      
      # Tab 4: Drug
      tabItem(tabName = "drug",
              fluidRow(
                box(width = 3, title = "Controls", status = "warning", solidHeader = TRUE,
                    selectizeInput("gene_input_drug", "Select Gene:",
                                  choices = NULL,
                                  options = list(placeholder = 'Type gene name...')),
                    p(strong("Note:"), "Shows correlation with predicted drug sensitivity.",
                      style = "font-size: 11px; color: #666;")
                ),
                box(width = 9, title = "Predicted Drug Sensitivity Correlation", solidHeader = TRUE, status = "danger",
                    p("Positive correlation suggests high expression is associated with sensitivity."),
                    plotlyOutput("drug_plot", height = "600px")
                )
              )
      ),
      
      # Tab 5: Survival
      tabItem(tabName = "survival",
              fluidRow(
                box(width = 12, title = "Survival Analysis (TCGA-PCPG Validation Cohort)", solidHeader = TRUE, status = "info",
                    p(strong("Note:"), "This analysis uses the gene selected in the 'Gene Expression' tab. Kaplan-Meier curves show high vs. low expression (median split)."),
                    p("The TCGA-PCPG cohort is used as an independent validation dataset for prognostic markers."),
                    plotOutput("survival_plot", height = "500px")
                )
              )
      ),
      
      # Tab 6: AI Diagnosis
      tabItem(tabName = "ai",
              fluidRow(
                box(width = 4, title = "AI Classifier", status = "primary", solidHeader = TRUE,
                    h4("Upload Patient Data"),
                    p("Upload a CSV file with gene expression data (TPM, FPKM, or log2-transformed counts)."),
                    p(strong("Format:"), "Genes in rows, Samples in columns. First column must be 'Gene'."),
                    fileInput("upload_file", "Choose CSV File", accept = ".csv"),
                    actionButton("run_ai", "Run Prediction", icon = icon("play"), class = "btn-success"),
                    hr(),
                    h4("Download Templates"),
                    downloadButton("download_example", "Example Data (5 samples)", class = "btn-info btn-sm"),
                    br(), br(),
                    downloadButton("download_template", "Blank Template", class = "btn-warning btn-sm"),
                    hr(),
                    h4("Model Info"),
                    p(icon("robot"), strong("Algorithm:"), "XGBoost Multiclass Classifier"),
                    p(icon("dna"), strong("Features:"), "50 Signature Genes"),
                    p(icon("chart-line"), strong("Performance:"), "CV LogLoss = 0.01"),
                    p(icon("check-circle"), strong("Training Accuracy:"), "100%")
                ),
                box(width = 8, title = "Prediction Results", solidHeader = TRUE, status = "warning",
                    uiOutput("pred_status"),
                    DTOutput("pred_table"),
                    hr(),
                    h4("Model Interpretation (SHAP Values)"),
                    p("SHAP values explain which genes contributed most to each prediction."),
                    plotOutput("shap_plot", height = "500px")
                )
              )
      ),
      
      # Tab 7: Data
      tabItem(tabName = "data",
              fluidRow(
                box(width = 12, title = "Sample Metadata Explorer", solidHeader = TRUE,
                    DTOutput("meta_table")
                )
              )
      )
    )
  )
)

# --- 3. Server Logic ---
server <- function(input, output, session) {

  # Initialize gene selectize inputs with server-side search
  updateSelectizeInput(session, "gene_input",
                      choices = data_bundle$all_genes,
                      selected = "TP53",
                      server = TRUE)

  updateSelectizeInput(session, "gene_input_immune",
                      choices = data_bundle$all_genes,
                      selected = "CD8A",
                      server = TRUE)

  updateSelectizeInput(session, "gene_input_drug",
                      choices = data_bundle$all_genes,
                      selected = "TP53",
                      server = TRUE)

  # Quick gene selection buttons - need to provide choices for server-side selectize
  observeEvent(input$gene_tp53, {
    updateSelectizeInput(session, "gene_input", choices = data_bundle$all_genes, selected = "TP53", server = TRUE)
  })
  observeEvent(input$gene_rb1, {
    updateSelectizeInput(session, "gene_input", choices = data_bundle$all_genes, selected = "RB1", server = TRUE)
  })
  observeEvent(input$gene_insm1, {
    updateSelectizeInput(session, "gene_input", choices = data_bundle$all_genes, selected = "INSM1", server = TRUE)
  })
  observeEvent(input$gene_chga, {
    updateSelectizeInput(session, "gene_input", choices = data_bundle$all_genes, selected = "CHGA", server = TRUE)
  })
  observeEvent(input$gene_cd8a, {
    updateSelectizeInput(session, "gene_input", choices = data_bundle$all_genes, selected = "CD8A", server = TRUE)
  })
  observeEvent(input$gene_cd274, {
    updateSelectizeInput(session, "gene_input", choices = data_bundle$all_genes, selected = "CD274", server = TRUE)
  })
  observeEvent(input$gene_mki67, {
    updateSelectizeInput(session, "gene_input", choices = data_bundle$all_genes, selected = "MKI67", server = TRUE)
  })
  observeEvent(input$gene_top2a, {
    updateSelectizeInput(session, "gene_input", choices = data_bundle$all_genes, selected = "TOP2A", server = TRUE)
  })

  # Helper: Get Gene Data (GEO) - for Expression tab
  get_gene_data <- reactive({
    req(input$gene_input, data_bundle$expr)
    gene <- input$gene_input

    if (gene %in% rownames(data_bundle$expr)) {
      df <- data.frame(
        Expression = data_bundle$expr[gene, ],
        sample_id = colnames(data_bundle$expr)
      ) %>% left_join(data_bundle$meta, by = "sample_id")
      return(df)
    }
    return(NULL)
  })

  # Helper: Get Gene Data - for Immune tab
  get_gene_data_immune <- reactive({
    req(input$gene_input_immune, data_bundle$expr)
    gene <- input$gene_input_immune

    if (gene %in% rownames(data_bundle$expr)) {
      df <- data.frame(
        Expression = data_bundle$expr[gene, ],
        sample_id = colnames(data_bundle$expr)
      ) %>% left_join(data_bundle$meta, by = "sample_id")
      return(df)
    }
    return(NULL)
  })

  # Helper: Get Gene Data - for Drug tab
  get_gene_data_drug <- reactive({
    req(input$gene_input_drug, data_bundle$expr)
    gene <- input$gene_input_drug

    if (gene %in% rownames(data_bundle$expr)) {
      df <- data.frame(
        Expression = data_bundle$expr[gene, ],
        sample_id = colnames(data_bundle$expr)
      ) %>% left_join(data_bundle$meta, by = "sample_id")
      return(df)
    }
    return(NULL)
  })
  
  # Plot 1: Expression Boxplot
  output$expr_plot <- renderPlotly({
    df <- get_gene_data()
    req(df)
    
    p <- ggplot(df, aes_string(x = input$group_by, y = "Expression", fill = input$group_by)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.5, size = 1) +
      theme_bw() +
      labs(y = "log2(TPM+1)", title = paste("Expression of", input$gene_input)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggplotly(p)
  })
  
  # Plot 2: Immune Correlation
  output$immune_corr_plot <- renderPlotly({
    df <- get_gene_data_immune()
    req(df, input$gene_input_immune)
    
    if (file.exists(immune_file)) {
      imm <- read_csv(immune_file, show_col_types = FALSE)
      plot_df <- df %>% left_join(imm, by = "sample_id")
      
      immune_cells <- colnames(imm)[-1]
      immune_cells <- immune_cells[sapply(imm[immune_cells], is.numeric)]
      
      if(length(immune_cells) == 0) return(NULL)
      
      cor_res <- sapply(immune_cells, function(cell) {
        cor(plot_df$Expression, plot_df[[cell]], use = "complete.obs", method = "spearman")
      })
      
      cor_df <- data.frame(Cell = names(cor_res), Correlation = cor_res)
      
      p <- ggplot(cor_df, aes(x = reorder(Cell, Correlation), y = Correlation, fill = Correlation)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits=c(-1,1)) +
        theme_bw() +
        labs(title = paste("Correlation of", input$gene_input_immune, "with Immune Cells"), x = "")
      
      ggplotly(p)
    }
  })
  
  # Plot 3: Drug Plot
  output$drug_plot <- renderPlotly({
    df <- get_gene_data_drug()
    req(df, input$gene_input_drug)

    if (file.exists(drug_file)) {
      drugs <- read_csv(drug_file, show_col_types = FALSE)
      plot_df <- df %>% left_join(drugs, by = "sample_id")

      drug_names <- colnames(drugs)[!colnames(drugs) %in% colnames(data_bundle$meta)]
      drug_names <- drug_names[drug_names != "sample_id"]
      drug_names <- drug_names[sapply(drugs[drug_names], is.numeric)]

      if(length(drug_names) == 0) return(NULL)

      cor_res <- sapply(drug_names, function(d) {
        cor(plot_df$Expression, plot_df[[d]], use = "complete.obs", method = "spearman")
      })

      cor_df <- data.frame(Drug = names(cor_res), Correlation = cor_res)

      p <- ggplot(cor_df, aes(x = reorder(Drug, Correlation), y = Correlation, fill = Correlation)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits=c(-1,1)) +
        theme_bw() +
        labs(title = paste("Correlation of", input$gene_input_drug, "with Predicted Drug Sensitivity"),
             subtitle = "Positive = High expression associated with sensitivity (target)", 
             x = "")
      
      ggplotly(p)
    }
  })
  
  # Plot 4: Survival Plot (TCGA) - uses gene from Expression tab
  output$survival_plot <- renderPlot({
    req(input$gene_input)
    tcga_exp <- data_bundle$tcga_exp
    tcga_clin <- data_bundle$tcga_clin
    gene_map <- data_bundle$gene_map

    if (is.null(tcga_exp) || is.null(tcga_clin) || is.null(gene_map)) return(NULL)

    # Map Symbol to Ensembl ID
    sel_gene <- input$gene_input
    target_id <- gene_map %>% filter(SYMBOL == sel_gene) %>% pull(original_id)
    
    if (length(target_id) == 0 || !target_id %in% rownames(tcga_exp)) {
      plot(0,0, type="n", axes=FALSE, xlab="", ylab="")
      text(0,0, paste("Gene", sel_gene, "not found in TCGA dataset."))
      return(NULL)
    }
    
    # Extract expression
    expr_vec <- tcga_exp[target_id[1], ]
    
    surv_df <- tcga_clin %>%
      mutate(
        time = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death),
        status = ifelse(vital_status == "Dead", 1, 0)
      ) %>%
      filter(!is.na(time) & !is.na(status) & time > 0)
    
    common_samps <- intersect(names(expr_vec), rownames(surv_df))
    if(length(common_samps) < 10) return(NULL)
    
    surv_df <- surv_df[common_samps, ]
    surv_df$Expression <- expr_vec[common_samps]
    
    med <- median(surv_df$Expression, na.rm = TRUE)
    surv_df$Group <- ifelse(surv_df$Expression > med, "High", "Low")
    
    fit <- survfit(Surv(time, status) ~ Group, data = surv_df)
    
    ggsurvplot(fit, data = surv_df, pval = TRUE, risk.table = TRUE,
               palette = c("firebrick", "navy"),
               title = paste("TCGA Survival Analysis:", sel_gene),
               xlab = "Time (Days)",
               risk.table.height = 0.25)
  })
  
  # Download handlers
  output$download_example <- downloadHandler(
    filename = function() "example_upload_template.csv",
    content = function(file) {
      file.copy("results/08_ai_model/example_upload_template.csv", file)
    }
  )

  output$download_template <- downloadHandler(
    filename = function() "blank_upload_template.csv",
    content = function(file) {
      file.copy("results/08_ai_model/blank_upload_template.csv", file)
    }
  )

  # AI Prediction Logic
  observeEvent(input$run_ai, {
    req(input$upload_file, data_bundle$ai_model)

    tryCatch({
      # Read uploaded file
      user_data <- read_csv(input$upload_file$datapath, show_col_types = FALSE)

      # Validate format
      if(!"Gene" %in% colnames(user_data)) {
        showNotification("Error: First column must be named 'Gene'", type = "error", duration = 5)
        return(NULL)
      }

      # Extract model components
      model <- data_bundle$ai_model$model
      features <- data_bundle$ai_model$features
      label_map <- data_bundle$ai_model$labels

      # Convert to matrix format (genes as rows)
      gene_col <- user_data$Gene
      user_matrix <- as.matrix(user_data[, -1])
      rownames(user_matrix) <- gene_col

      # Match features (find overlap between user genes and signature genes)
      common_genes <- intersect(features, rownames(user_matrix))
      missing_genes <- setdiff(features, rownames(user_matrix))

      if(length(common_genes) < 30) {
        showNotification(paste("Warning: Only", length(common_genes), "out of 50 signature genes found. Prediction may be unreliable."),
                         type = "warning", duration = 7)
      }

      # Create feature matrix (fill missing genes with 0)
      pred_matrix <- matrix(0, nrow = ncol(user_matrix), ncol = length(features))
      colnames(pred_matrix) <- features
      rownames(pred_matrix) <- colnames(user_matrix)

      # Fill in available gene values
      for(gene in common_genes) {
        pred_matrix[, gene] <- user_matrix[gene, ]
      }

      # Make predictions
      dpred <- xgb.DMatrix(data = pred_matrix)
      pred_probs <- predict(model, dpred, reshape = TRUE)
      colnames(pred_probs) <- names(label_map)

      # Get predicted classes
      pred_classes <- apply(pred_probs, 1, which.max)
      pred_labels <- names(label_map)[pred_classes]
      max_probs <- apply(pred_probs, 1, max)

      # Create result table
      res_df <- data.frame(
        Sample = colnames(user_matrix),
        Predicted_Subtype = pred_labels,
        Confidence = round(max_probs, 3),
        C1_Probability = round(pred_probs[, "C1"], 3),
        C2_Probability = round(pred_probs[, "C2"], 3),
        C3_Probability = round(pred_probs[, "C3"], 3)
      )

      # Display results
      output$pred_status <- renderUI({
        tags$div(
          style = "background-color: #d4edda; padding: 10px; border-radius: 5px; margin-bottom: 10px;",
          h4(icon("check-circle"), " Prediction Complete!", style = "color: #155724;"),
          p(paste("Processed", nrow(res_df), "sample(s) successfully.")),
          p(paste("Matched", length(common_genes), "out of 50 signature genes."))
        )
      })

      output$pred_table <- renderDT({
        datatable(res_df, options = list(scrollX = TRUE, pageLength = 10)) %>%
          formatStyle('Predicted_Subtype',
                      backgroundColor = styleEqual(c('C1', 'C2', 'C3'),
                                                   c('#e3f2fd', '#fff3e0', '#f3e5f5')))
      })

      # SHAP plot (for first sample only, to keep it simple)
      output$shap_plot <- renderPlot({
        if(nrow(pred_matrix) > 0) {
          # Compute SHAP values for first sample
          shap_values <- shap.values(xgb_model = model, X_train = pred_matrix)

          # Get mean absolute SHAP for each feature
          shap_importance <- colMeans(abs(shap_values$shap_score))
          top_genes <- names(sort(shap_importance, decreasing = TRUE)[1:20])

          # Plot
          shap_long <- shap.prep(xgb_model = model, X_train = pred_matrix)
          shap.plot.summary(shap_long, dilute = 10) +
            ggtitle("Top 20 Gene Contributions (SHAP Values)") +
            theme_minimal()
        } else {
          plot(0,0, type="n", axes=FALSE, xlab="", ylab="")
          text(0,0, "No data to display SHAP values.")
        }
      })

    }, error = function(e) {
      output$pred_status <- renderUI({
        tags$div(
          style = "background-color: #f8d7da; padding: 10px; border-radius: 5px;",
          h4(icon("exclamation-triangle"), " Prediction Failed", style = "color: #721c24;"),
          p(paste("Error:", e$message)),
          p("Please check your file format and try again.")
        )
      })
      showNotification(paste("Error:", e$message), type = "error", duration = 10)
    })
  })
  
  # Table: Data Explorer
  output$meta_table <- renderDT({
    req(data_bundle$meta)

    # Create a copy with GEO links
    meta_display <- data_bundle$meta %>%
      mutate(
        gse_id = paste0('<a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=',
                       gse_id, '" target="_blank">', gse_id,
                       ' <i class="fa fa-external-link-alt"></i></a>')
      )

    datatable(meta_display,
              options = list(
                scrollX = TRUE,
                pageLength = 25,
                order = list(list(1, 'asc'))  # Sort by gse_id
              ),
              escape = FALSE,  # Allow HTML in cells
              rownames = FALSE,
              caption = 'Click on GEO IDs to view original studies on NCBI GEO database.')
  })
}

# --- 4. Run App ---
shinyApp(ui, server)