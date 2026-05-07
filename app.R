# Call the libraries needed for the shiny app and and all the steps.
library(shiny)
library(bslib)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(DT)
library(dplyr)
library(multcompView)
library(rstatix)
library(ggpubr)
library(factoextra)
library(tibble)
library(agricolae)
library(pheatmap)
library(grid)
library(writexl)
library(shinyjqui)
ui <- page_sidebar(
  tags$head(tags$link(rel = "stylesheet", 
                      href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.15.4/css/all.min.css")),
  title = "Gene Expression Analysis (from qPCR or Array chip)",
  sidebar = sidebar(
    fileInput(
      "data", 
      "Upload gene expression data (CSV or TSV file)", 
      accept = c(".csv", ".tsv", "text/csv", "text/tab-separated-values")
    ),
    downloadButton("download_selected_genes_pdf", "Download Selected Gene Plots (PDF)"),
    
    
    selectInput(
      "impute_method", 
      "Imputation Method:",
      choices = c("Mean" = "mean", "Median" = "median", "Zero" = "zero"),
      selected = "mean"
    ),
    
    selectInput("outlier_replace_method", "Outlier Replacement Method:",
                choices = c("Group Mean" = "mean", "Group Median" = "median"),
                selected = "mean"
    ),
    checkboxInput("outlier_replace_global", 
                  "Use global (all samples for gene) replacement instead of per-treatment group?", 
                  value = FALSE),
    checkboxInput(
      inputId = "use_fdr",
      label = "Use FDR correction for significant gene selection (Benjamini-Hochberg)",
      value = FALSE
    )
    ,
    actionButton("show_about", "About / Instructions", icon = icon("info-circle"), class = "btn-info", style = "margin-bottom: 10px; width: 100%;"),
    
    actionButton("show_outlier_replacement", "Show Outlier Replacement Summary"),
    
    ## functions of bottoms that will show up
    
    radioButtons("figure_type", "Choose Figure type:",
                 choices = c("PCA (All Genes)", 
                             "Biplot PCA (All Genes)", 
                             "PCA (Significant Genes)", 
                             "Biplot PCA (Significant Genes)", 
                             "Boxplot", 
                             "Boxplot_KWT",
                             "Boxplot with Dots", 
                             "Violin Plot", 
                             "Boxplot & Violin Plot",
                             "Outlier Boxplot",
                             "Heatmap (All Genes)",
                             "Custom qPCR-style Plot",
                             "Custom qPCR-style Plot (Parametric)",
                             "Custom qPCR-style Plot (Raw)",
                             "Custom qPCR-style Plot (Parametric, Raw)",
                             "Heatmap (Significant Genes)"),
                 selected = "PCA (All Genes)"),
    
    # For PCA plots: choose PC axes
    conditionalPanel(
      condition = "input.figure_type.includes('PCA')",
      selectInput("pc_x", "Select PC for X-axis", choices = NULL),
      selectInput("pc_y", "Select PC for Y-axis", choices = NULL)
    ),
    
    # For non-PCA plots (including both Heatmaps): choose gene(s)
    conditionalPanel(
      condition = "['Boxplot','Boxplot_KWT','Boxplot with Dots','Violin Plot','Boxplot & Violin Plot','Custom qPCR-style Plot',
        'Custom qPCR-style Plot (Parametric)','Custom qPCR-style Plot (Raw)',
        'Custom qPCR-style Plot (Parametric, Raw)','Outlier Boxplot','Heatmap (All Genes)',
        'Heatmap (Significant Genes)'].includes(input.figure_type)",
      selectizeInput("selected_genes", "Select gene(s) for plot", choices = NULL, multiple = TRUE)
    ),
    
    # Axis label inputs for non-PCA, non-heatmap
    conditionalPanel(
      condition = "!input.figure_type.includes('PCA') && !input.figure_type.includes('Heatmap')",
      textInput("x_label", "X-axis Label:", value = "Treatment"),
      textInput("y_label", "Y-axis Label:", value = "Expression Level")
    ),
    
    # Figure title and text size controls (applies to all plots)
    textInput("figure_title", "Figure Title:", value = ""),
    numericInput("title_size", "Title Text Size:", value = 16, min = 8, max = 40, step = 1),
    numericInput("axis_title_size", "Axis Title Text Size:", value = 14, min = 8, max = 40, step = 1),
    numericInput("axis_text_size", "Axis Text Size:", value = 12, min = 8, max = 40, step = 1),
    
    # Heatmap-specific: clustering + grouping
    conditionalPanel(
      condition = "input.figure_type.includes('Heatmap')",
      checkboxInput("cluster_rows", "Cluster rows?", value = TRUE),
      checkboxInput("cluster_cols", "Cluster columns?", value = TRUE),
      checkboxInput("do_scale",   "Scale data (z-score)?", value = TRUE),
      selectInput("group_col",   "Select grouping column (optional):", choices = c("None")),
      checkboxInput("aggregate_heatmap", "Aggregate by group means?", value = FALSE),
      checkboxInput("transpose_heatmap",
                    "Transpose heatmap (genes ↔ samples)?",
                    value = FALSE)
      
    ),
    
    radioButtons(
      "sig_gene_method",
      "Significant Genes: which test?",
      choices = c(
        "Automatic (normality)" = "auto",
        "Always ANOVA" = "anova",
        "Always Kruskal-Wallis" = "kwt"
      ),
      selected = "auto"
    ),
    
    
    
    downloadButton("download_figure", "Download Figure"),
    downloadButton("download_all_figures_pdf", "Download All Figures (PDF)"),
    downloadButton("download_excel_stats", "Download Excel Summary"),
    
    selectInput("bg_color", "Background Color:", 
                choices = list("White" = "white", "Gray" = "gray90", "Black" = "black"),
                selected = "white"),
    
    actionButton("show_sig_genes", "Show Significant Genes", 
                 class = "btn-primary", 
                 style = "margin-top: 10px; width: 100%;"),
    actionButton("show_outliers", "Show Outliers", 
                 class = "btn-secondary", 
                 style = "margin-top: 10px; width: 100%;"),
    actionButton("show_imputed", "Show Imputed Genes", 
                 class = "btn-warning", 
                 style = "margin-top: 10px; width: 100%;"),
    actionButton("show_normality", "Show Normality Test Results", 
                 class = "btn-info", 
                 style = "margin-top: 10px; width: 100%;"),
    actionButton("replace_outliers", "Replace Outliers", 
                 class = "btn-danger", 
                 style = "margin-top: 10px; width: 100%;"),
    actionButton("show_modified_data", "Show Modified Data", 
                 class = "btn-success", 
                 style = "margin-top: 10px; width: 100%;"),
    actionButton("show_original_data", "Show Original Data", 
                 class = "btn-secondary", 
                 style = "margin-top: 10px; width: 100%;"),
    checkboxInput("global_outlier_detection", 
                  "Use global outlier detection (across all groups)", 
                  value = TRUE)
  ),
  splitLayout(
    cellWidths = c("65%", "35%"),   # Or use "50%, 50%" if you want equal
    card(
      card_header("Plot Output"),
      plotOutput("main_plot", height = "calc(80vh - 60px)")
    ),
    card(
      card_header("Analysis Information"),
      card_body(
        accordion(
          accordion_panel("PCA Summary", verbatimTextOutput("pca_summary")),
          accordion_panel("Data Preview", DTOutput("data_preview")),
          accordion_panel("Significant Genes", DTOutput("sig_genes_table")),
          accordion_panel("Outliers Detected", DTOutput("outlier_table"))
        )
      )
    )
  ),
  
  # Footer card
  card(
    full_screen = FALSE,
    card_footer(
      HTML(
        paste(
          "Developed by Janan Gawra",
          "<a href='https://www.linkedin.com/in/janangawra/' target='_blank'>",
          "<i class='fab fa-linkedin'></i> LinkedIn Profile</a>",
          sep = " | "
        )
      ),
      style = "text-align: right; font-size: 0.9em; color: blue;"
    )
  )
)
##################################
server <- function(input, output, session) {
  
  ###### Data Reading with Missing Data Handling ######
  raw_df <- reactive({
    req(input$data)
    tryCatch({
      df <- read.csv(input$data$datapath, na.strings = c("", "NA", "na", "NaN"))
      df
    }, error = function(e) {
      showNotification("Error reading file: Make sure it is a CSV with the right format.", type = "error")
      NULL
    })
  })
  
  
  # Expression data (numeric only) with imputation by group means
  expr_data <- reactive({
    df <- raw_df()
    expression_data <- df[, -c(1:2), drop = FALSE]
    expression_data[] <- lapply(expression_data, as.numeric)
    rownames(expression_data) <- df[[1]]
    treatments <- df[[2]]
    attr(expression_data, "treatment") <- treatments
    
    impute_method <- input$impute_method
    imputed_data <- expression_data
    imputed_list <- list()
    
    for (gene in colnames(expression_data)) {
      for (tr in unique(treatments)) {
        idx <- which(treatments == tr)
        missing_idx <- which(is.na(expression_data[idx, gene]))
        group_vals <- expression_data[idx, gene]
        if (length(missing_idx) > 0) {
          # Choose method
          impute_value <- NA
          if (impute_method == "mean") {
            impute_value <- mean(group_vals, na.rm = TRUE)
          } else if (impute_method == "median") {
            impute_value <- median(group_vals, na.rm = TRUE)
          } else if (impute_method == "zero") {
            impute_value <- 0
          }
          if (!is.na(impute_value)) {
            imputed_data[idx[missing_idx], gene] <- impute_value
            imputed_list[[length(imputed_list) + 1]] <- data.frame(
              Gene = gene, Treatment = tr, MissingCount = length(missing_idx),
              ImputedWith = impute_method, Value = impute_value,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    imputed_summary <- if(length(imputed_list)) do.call(rbind, imputed_list) else
      data.frame(Gene=character(), Treatment=character(), MissingCount=integer(), ImputedWith=character(), Value=numeric())
    attr(imputed_data, "imputed_summary") <- imputed_summary
    imputed_data
  })
  
  
  # Reactive flag for outlier replacement
  replace_outliers_flag <- reactiveVal(FALSE)
  observeEvent(input$replace_outliers, {
    replace_outliers_flag(TRUE)
  })
  
  # Modified data reactive: if outlier replacement is selected, replace outlier values with group mean
  mod_data <- reactive({
    df_expr <- expr_data()
    treatments <- attr(df_expr, "treatment")
    outliers <- detect_outliers()
    modified <- df_expr
    
    replace_method <- input$outlier_replace_method
    use_global <- isTRUE(input$outlier_global)
    replacement_fun <- if (replace_method == "mean") mean else median
    
    replaced <- data.frame(Gene = character(), Sample = character(), OldValue = numeric(), NewValue = numeric(), stringsAsFactors = FALSE)
    
    if (isTruthy(input$replace_outliers) && nrow(outliers) > 0) {
      for (i in seq_len(nrow(outliers))) {
        gene <- outliers$Gene[i]
        sample_id <- outliers$Sample[i]
        treatment <- outliers$Treatment[i]
        sample_row <- which(rownames(modified) == sample_id)
        
        if (length(sample_row) == 1) {
          # --- NEW: Get ALL outlier samples for this gene (and group) ---
          if (use_global) {
            # All outlier samples for this gene (global replacement)
            outlier_samples_gene <- outliers$Sample[outliers$Gene == gene]
            values <- df_expr[, gene]
            clean_values <- values[!(rownames(df_expr) %in% outlier_samples_gene)]
            new_val <- replacement_fun(clean_values, na.rm = TRUE)
          } else {
            # All outlier samples for this gene IN THIS GROUP
            group_idx <- which(treatments == treatment)
            group_sample_names <- rownames(df_expr)[group_idx]
            outlier_samples_group <- outliers$Sample[outliers$Gene == gene & outliers$Treatment == treatment]
            group_vals <- df_expr[group_idx, gene]
            clean_group_vals <- group_vals[!(group_sample_names %in% outlier_samples_group)]
            new_val <- replacement_fun(clean_group_vals, na.rm = TRUE)
          }
          old_val <- modified[sample_row, gene]
          modified[sample_row, gene] <- new_val
          
          replaced <- rbind(replaced, data.frame(
            Gene = gene, Sample = sample_id, OldValue = old_val, NewValue = new_val, stringsAsFactors = FALSE
          ))
        }
      }
    }
    attr(modified, "treatment") <- treatments
    replaced_outliers(replaced)   # update the tracked changes
    return(modified)
  })
  
  normality_data <- reactive({
    df <- mod_data()
    treatments <- attr(df, "treatment")
    
    results <- lapply(colnames(df), function(gene) {
      values <- df[[gene]]
      tibble(
        Gene = gene,
        W_statistic = shapiro.test(values)$statistic,
        p_value = shapiro.test(values)$p.value,
        Normal = ifelse(shapiro.test(values)$p.value >= 0.05, "Yes", "No")
      )
    })
    
    bind_rows(results)
  })
  
  analysis_type_data <- reactive({
    df <- mod_data()
    treatments <- attr(df, "treatment")
    
    results <- lapply(colnames(df), function(gene) {
      values <- df[[gene]]
      is_normal <- shapiro.test(values)$p.value >= 0.05
      test_used <- if (is_normal) "ANOVA + Tukey" else "Kruskal-Wallis + Dunn"
      
      tibble(
        Gene = gene,
        Normal_Distribution = ifelse(is_normal, "Yes", "No"),
        Test_Used = test_used
      )
    })
    
    bind_rows(results)
  })
  
  pairwise_results <- reactive({
    df <- mod_data()
    treatments <- factor(attr(df, "treatment"))
    all_results <- lapply(colnames(df), function(gene) {
      values <- df[[gene]]
      plot_data <- data.frame(Expression = values, Treatment = treatments)
      # t-test
      t_res <- tryCatch({
        rstatix::pairwise_t_test(Expression ~ Treatment, data = plot_data, p.adjust.method = "BH") %>%
          mutate(Gene = gene, Test = "t_test") %>%
          select(Gene, Test, group1, group2, p, p.adj, p.adj.signif)
      }, error=function(e) NULL)
      # Wilcoxon
      w_res <- tryCatch({
        rstatix::pairwise_wilcox_test(Expression ~ Treatment, data = plot_data, p.adjust.method = "BH") %>%
          mutate(Gene = gene, Test = "wilcox") %>%
          select(Gene, Test, group1, group2, p, p.adj, p.adj.signif)
      }, error=function(e) NULL)
      # Tukey (agricolae)
      tk_res <- tryCatch({
        fit <- aov(Expression ~ Treatment, data = plot_data)
        hsd <- agricolae::HSD.test(fit, "Treatment")
        comp <- hsd$comparison
        if (!is.null(comp) && is.matrix(comp) && nrow(comp) > 1) {
          # convert to data.frame
          comp_df <- as.data.frame(as.table(comp))
          names(comp_df) <- c("group1", "group2", "p")
          comp_df$Gene <- gene
          comp_df$Test <- "Tukey"
          comp_df$p.adj <- comp_df$p
          comp_df$p.adj.signif <- ifelse(comp_df$p < 0.05, "*", "ns")
          comp_df[,c("Gene","Test","group1","group2","p","p.adj","p.adj.signif")]
        } else {
          NULL # will be skipped
        }
      }, error = function(e) NULL)
      
      # Dunn (rstatix)
      dn_res <- tryCatch({
        dunn <- rstatix::dunn_test(Expression ~ Treatment, data = plot_data, p.adjust.method = "BH")
        if (nrow(dunn) > 0) {
          dunn %>%
            mutate(Gene = gene, Test = "Dunn") %>%
            select(Gene, Test, group1, group2, p, p.adj, p.adj.signif)
        } else {
          NULL
        }
      }, error=function(e) NULL)
      # Combine
      bind_rows(t_res, w_res, tk_res, dn_res)
    })
    bind_rows(all_results)
  })
  
  pairwise_tukey_results <- reactive({
    df <- mod_data()
    treatments <- factor(attr(df, "treatment"))
    
    # Only keep genes considered significant under current selection
    sig_genes <- find_significant_genes()$Gene
    
    all_results <- lapply(sig_genes, function(gene) {
      plot_data <- data.frame(Expression = df[[gene]], Treatment = treatments)
      
      tukey_res <- tryCatch({
        model <- aov(Expression ~ Treatment, data = plot_data)
        res <- TukeyHSD(model)
        broom::tidy(res) %>%
          tidyr::separate(contrast, into = c("group1", "group2"), sep = "-") %>%
          mutate(
            Gene = gene,
            Test = "Tukey",
            p = adj.p.value,
            p.adj = adj.p.value,
            p.adj.signif = case_when(
              adj.p.value <= 0.001 ~ "***",
              adj.p.value <= 0.01  ~ "**",
              adj.p.value <= 0.05  ~ "*",
              TRUE ~ "ns"
            )
          ) %>%
          select(Gene, Test, group1, group2, p, p.adj, p.adj.signif)
      }, error = function(e) {
        tibble(Gene=gene, Test="Tukey", group1=NA, group2=NA, p=NA, p.adj=NA, p.adj.signif=NA)
      })
      
      tukey_res
    })
    
    bind_rows(all_results)
  })
  
  
  
  # --- Add here! ---
  all_genes_pvalues <- reactive({
    df <- mod_data()
    treatments <- attr(df, "treatment")
    anova_p <- kwt_p <- rep(NA_real_, ncol(df))
    
    for (i in seq_along(colnames(df))) {
      values <- df[[i]]
      # ANOVA
      anova_p[i] <- tryCatch({
        summary(aov(values ~ treatments))[[1]]$"Pr(>F)"[1]
      }, error=function(e) NA)
      # Kruskal-Wallis
      kwt_p[i] <- tryCatch({
        kruskal.test(values ~ treatments)$p.value
      }, error=function(e) NA)
    }
    
    # FDR (Benjamini-Hochberg) adjustments
    anova_p_adj <- p.adjust(anova_p, method = "BH")
    kwt_p_adj <- p.adjust(kwt_p, method = "BH")
    
    tibble(
      Gene = colnames(df),
      ANOVA_p = anova_p,
      ANOVA_p_adj = anova_p_adj,
      KWT_p = kwt_p,
      KWT_p_adj = kwt_p_adj
    )
  })
  
  
  all_genes_pvalues_fdr <- reactive({
    dfm <- mod_data()
    tr  <- attr(dfm, "treatment")
    results <- tibble(Gene=colnames(dfm), p_value=NA_real_, method=NA_character_)
    
    for(i in seq_along(colnames(dfm))){
      g <- colnames(dfm)[i]
      vals <- dfm[[g]]
      norm <- shapiro.test(vals)$p.value >= .05
      if(norm){
        pval <- summary(aov(vals ~ tr))[[1]]$"Pr(>F)"[1]
        results$p_value[i] <- pval
        results$method[i] <- "ANOVA"
      } else {
        pval <- kruskal.test(vals ~ tr)$p.value
        results$p_value[i] <- pval
        results$method[i] <- "Kruskal-Wallis"
      }
    }
    results$adj_p_BH <- p.adjust(results$p_value, method = "BH")
    results$adj_p_Bonferroni <- p.adjust(results$p_value, method = "bonferroni")
    results
  })
  
  # Identify factor/character columns for grouping
  observe({
    df <- raw_df()
    fc <- names(df)[sapply(df, function(x) is.character(x) || is.factor(x))]
    choices_vec <- c("None", fc)
    named_choices <- setNames(as.list(choices_vec), choices_vec)
    
    updateSelectInput(
      session, "group_col",
      choices  = named_choices,
      selected = "None"
    )
  })
  
  
  ###### Helper functions ######
  get_treatment_colors <- reactive({
    req(mod_data())
    treatments <- unique(attr(mod_data(), "treatment"))
    n_treatments <- length(treatments)
    colors <- c("steelblue", colorRampPalette(c("#ffcccc", "#ff0000"))(n_treatments - 1))
    names(colors) <- treatments
    return(colors)
  })
  
  find_significant_genes <- reactive({
    dfm    <- mod_data()
    tr     <- attr(dfm, "treatment")
    genes  <- colnames(dfm)
    nG     <- length(genes)
    
    # 1) Precompute *all* ANOVA & KW p-values
    anova_p <- map_dbl(genes, ~ {
      vals <- dfm[[.x]]
      tryCatch(summary(aov(vals ~ tr))[[1]]$"Pr(>F)"[1], error = function(e) NA_real_)
    })
    kwt_p <- map_dbl(genes, ~ {
      vals <- dfm[[.x]]
      tryCatch(kruskal.test(vals ~ tr)$p.value, error = function(e) NA_real_)
    })
    
    # 2) Normality only needed for "auto" mode
    shapiro_p <- map_dbl(genes, ~ {
      vals <- dfm[[.x]]
      tryCatch(shapiro.test(vals)$p.value, error = function(e) NA_real_)
    })
    norm <- shapiro_p >= 0.05
    
    use_fdr <- isTRUE(input$use_fdr)
    method  <- input$sig_gene_method  # "auto", "anova" or "kwt"
    
    if (method == "anova") {
      raw_p <- anova_p
      adj_p <- if(use_fdr) p.adjust(anova_p, "BH") else anova_p
      keep  <- raw_p < 0.05
      keep  <- if(use_fdr) adj_p < 0.05 else raw_p < 0.05
      
    } else if (method == "kwt") {
      raw_p <- kwt_p
      adj_p <- if(use_fdr) p.adjust(kwt_p, "BH") else kwt_p
      keep  <- raw_p < 0.05
      keep  <- if(use_fdr) adj_p < 0.05 else raw_p < 0.05
      
    } else {  # auto
      raw_p <- ifelse(norm, anova_p, kwt_p)
      adj_p <- ifelse(norm, p.adjust(anova_p, "BH"), p.adjust(kwt_p, "BH"))
      keep  <- raw_p < 0.05
      if (use_fdr) keep <- adj_p < 0.05
    }
    
    
    tibble(
      Gene  = genes,
      Method = toupper(method),
      RawP   = raw_p,
      AdjP   = adj_p
    ) %>%
      filter(keep) %>%
      arrange(if(use_fdr) AdjP else RawP)
  })
  
  
  
  detect_outliers <- reactive({
    req(expr_data())
    df_expr    <- expr_data()
    treatments <- attr(df_expr, "treatment")
    
    outliers <- data.frame(
      Gene      = character(),
      Sample    = character(),
      Treatment = character(),
      Value     = numeric(),
      stringsAsFactors = FALSE
    )
    
    for (gene in colnames(df_expr)) {
      values <- df_expr[, gene]
      
      if (input$global_outlier_detection) {
        ### ✅ GLOBAL IQR METHOD (across all groups)
        lower_bounds <- c()
        upper_bounds <- c()
        
        for (tr in unique(treatments)) {
          idx   <- which(treatments == tr)
          vals  <- df_expr[idx, gene]
          
          if (length(vals) >= 3 && !all(is.na(vals))) {
            Q1     <- quantile(vals, 0.25, na.rm = TRUE)
            Q3     <- quantile(vals, 0.75, na.rm = TRUE)
            IQRval <- Q3 - Q1
            lower  <- Q1 - 1.5 * IQRval
            upper  <- Q3 + 1.5 * IQRval
            lower_bounds <- c(lower_bounds, lower)
            upper_bounds <- c(upper_bounds, upper)
          }
        }
        
        global_lower <- min(lower_bounds, na.rm = TRUE)
        global_upper <- max(upper_bounds, na.rm = TRUE)
        
        out_idx <- which(values < global_lower | values > global_upper)
        if (length(out_idx) > 0) {
          outliers <- rbind(outliers,
                            data.frame(
                              Gene      = gene,
                              Sample    = rownames(df_expr)[out_idx],
                              Treatment = treatments[out_idx],
                              Value     = values[out_idx],
                              stringsAsFactors = FALSE
                            )
          )
        }
        
      } else {
        ### ✅ PER-GROUP IQR METHOD (treatment-specific)
        for (tr in unique(treatments)) {
          idx   <- which(treatments == tr)
          vals  <- df_expr[idx, gene, drop = FALSE]
          
          Q1     <- quantile(vals, 0.25, na.rm = TRUE)
          Q3     <- quantile(vals, 0.75, na.rm = TRUE)
          IQRval <- Q3 - Q1
          lower  <- Q1 - 1.5 * IQRval
          upper  <- Q3 + 1.5 * IQRval
          
          out_idx <- idx[which(df_expr[idx, gene] < lower | df_expr[idx, gene] > upper)]
          if (length(out_idx) > 0) {
            outliers <- rbind(outliers,
                              data.frame(
                                Gene      = gene,
                                Sample    = rownames(df_expr)[out_idx],
                                Treatment = treatments[out_idx],
                                Value     = df_expr[out_idx, gene],
                                stringsAsFactors = FALSE
                              )
            )
          }
        }
      }
    }
    
    return(outliers)
  })
  
  
  
  # Reactive normality test (Shapiro–Wilk) per gene and treatment group
  # Insert the new extended normality_results:
  normality_results <- reactive({
    req(mod_data())
    df_expr    <- mod_data()
    treatments <- attr(df_expr, "treatment")
    norm_summary <- data.frame(
      Gene       = character(),
      Treatment  = character(),
      SampleSize = numeric(),
      p_value    = numeric(),
      Normal     = character(),
      stringsAsFactors = FALSE
    )
    
    for (gene in colnames(df_expr)) {
      # by-group
      for (tr in unique(treatments)) {
        idx    <- which(treatments == tr)
        values <- df_expr[idx, gene]
        n      <- length(values)
        if (n >= 3 && !all(is.na(values))) {
          test      <- shapiro.test(values)
          p_val     <- test$p.value
          normality <- ifelse(p_val > 0.05, "Yes", "No")
        } else {
          p_val     <- NA
          normality <- "Insufficient Data"
        }
        norm_summary <- rbind(norm_summary,
                              data.frame(
                                Gene       = gene,
                                Treatment  = tr,
                                SampleSize = n,
                                p_value    = p_val,
                                Normal     = normality,
                                stringsAsFactors = FALSE
                              ))
      }
      
      # overall (all samples)
      all_vals <- df_expr[, gene]
      n_all    <- length(all_vals)
      if (n_all >= 3 && !all(is.na(all_vals))) {
        test_all      <- shapiro.test(all_vals)
        p_all         <- test_all$p.value
        normal_all    <- ifelse(p_all > 0.05, "Yes", "No")
      } else {
        p_all      <- NA
        normal_all <- "Insufficient Data"
      }
      norm_summary <- rbind(norm_summary,
                            data.frame(
                              Gene       = gene,
                              Treatment  = "All Samples",
                              SampleSize = n_all,
                              p_value    = p_all,
                              Normal     = normal_all,
                              stringsAsFactors = FALSE
                            ))
    }
    
    return(norm_summary)
  })
  
  ###### Observers for modals ######
  observeEvent(input$show_sig_genes, {
    sig_genes_info <- find_significant_genes()
    if(nrow(sig_genes_info) == 0) {
      showModal(modalDialog(
        title = "No Significant Genes",
        "No genes were found to be significantly different between treatments (p < 0.05)",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    } else {
      showModal(modalDialog(
        title = "Significant Genes",
        DTOutput("sig_genes_modal"),
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      output$sig_genes_modal <- renderDT({
        sig_genes_info$P_Value <- format(sig_genes_info$P_Value, scientific = TRUE, digits = 3)
        datatable(sig_genes_info,
                  options = list(pageLength = 10, scrollX = TRUE),
                  caption = "Significant Genes (p < 0.05)")
      })
    }
  })
  
  observeEvent(input$show_outliers, {
    outlier_data <- detect_outliers()
    if(nrow(outlier_data) == 0) {
      showModal(modalDialog(
        title = "No Outliers Detected",
        "No outliers were detected in the data.",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    } else {
      showModal(modalDialog(
        title = "Outliers Detected",
        DTOutput("outliers_modal"),
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      output$outliers_modal <- renderDT({
        datatable(outlier_data, 
                  options = list(pageLength = 10, scrollX = TRUE),
                  caption = "Outliers Detected")
      })
    }
  })
  
  observeEvent(input$show_imputed, {
    imputed_summary <- attr(expr_data(), "imputed_summary")
    if(nrow(imputed_summary) == 0) {
      showModal(modalDialog(
        title = "No Imputation Needed",
        "There were no missing values to impute.",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    } else {
      showModal(modalDialog(
        title = "Imputed Genes Summary",
        DTOutput("imputed_modal"),
        size = "l",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
      output$imputed_modal <- renderDT({
        datatable(imputed_summary, 
                  options = list(pageLength = 10, scrollX = TRUE),
                  caption = "Genes with missing values replaced by group means")
      })
    }
  })
  
  observeEvent(input$show_normality, {
    showModal(modalDialog(
      title = "Normality Test Results",
      DTOutput("normality_modal"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
    output$normality_modal <- renderDT({
      nr <- normality_results()
      nr$p_value <- format(nr$p_value, scientific = TRUE, digits = 3)
      datatable(nr, 
                options = list(pageLength = 10, scrollX = TRUE),
                caption = "Normality Test Results for each Gene and Treatment")
    })
  })
  
  
  # Observer for showing Modified Data (with outliers replaced)
  observeEvent(input$show_modified_data, {
    showModal(modalDialog(
      title = "Modified Data (Outliers Replaced)",
      DTOutput("modified_data_modal"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
    output$modified_data_modal <- renderDT({
      datatable(mod_data(), options = list(scrollX = TRUE, pageLength = 10))
    })
  })
  
  replaced_outliers <- reactiveVal(data.frame(
    Gene = character(), 
    Sample = character(), 
    OldValue = numeric(), 
    NewValue = numeric(), 
    stringsAsFactors = FALSE
  ))
  
  
  
  # Observer for showing Modified Data (with outliers replaced)
  observeEvent(input$show_outlier_replacement, {
    showModal(modalDialog(
      title = "Outlier Replacement Summary",
      DTOutput("outlier_replacement_modal"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
    output$outlier_replacement_modal <- renderDT({
      datatable(replaced_outliers(), options = list(scrollX = TRUE, pageLength = 10))
    })
  })
  
  
  # Observer for showing Original Data (with outliers intact)
  observeEvent(input$show_original_data, {
    showModal(modalDialog(
      title = "Original Data (With Outliers)",
      DTOutput("original_data_modal"),
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
    output$original_data_modal <- renderDT({
      datatable(raw_df(), options = list(scrollX = TRUE, pageLength = 10))
    })
  })
  
  ###### PCA Calculation ######
  pca_result <- reactive({
    req(mod_data())
    df_expr <- mod_data()
    if(grepl("Significant", input$figure_type)) {
      sig_genes <- find_significant_genes()$Gene
      if(length(sig_genes) == 0) return(NULL)
      scaled_data <- scale(df_expr[, sig_genes, drop = FALSE])
    } else {
      scaled_data <- scale(df_expr)
    }
    prcomp(scaled_data)
  })
  
  
  # Populate PC axes and gene selection
  observe({
    req(mod_data())
    updateSelectizeInput(session, "selected_genes", choices = colnames(mod_data()))
    updateSelectInput(session, "pc_x", choices = paste0("PC", 1:ncol(mod_data())), selected = "PC1")
    updateSelectInput(session, "pc_y", choices = paste0("PC", 1:ncol(mod_data())), selected = "PC2")
  })
  
  
  ###### Plot Generation ######
  generate_plot <- function() {
    req(input$figure_type)
    figure_type <- input$figure_type
    
    # 1) PCA Plots
    if (grepl("PCA", figure_type)) {
      pr <- pca_result()
      if(is.null(pr)) {
        return(
          ggplot() + 
            annotate("text", x = 0.5, y = 0.5, label = "No significant genes found") +
            theme_void()
        )
      }
      df_PCA <- data.frame(pr$x)
      df_PCA$Treatment <- attr(mod_data(), "treatment")
      
      if(grepl("Biplot", figure_type)) {
        p <- fviz_pca_biplot(
          pr, repel = TRUE, label = "all", 
          habillage = df_PCA$Treatment, 
          addEllipses = TRUE,
          palette = get_treatment_colors(), 
          title = "PCA Plot"
        )
      } else {
        p <- fviz_pca_ind(
          pr, 
          habillage = df_PCA$Treatment,
          addEllipses = TRUE,
          palette = get_treatment_colors(),
          title = "PCA Plot"
        )
      }
      if(nchar(input$figure_title) > 0) {
        p <- p + labs(title = input$figure_title)
      }
      p <- p + theme(
        plot.title = element_text(size = input$title_size),
        axis.title = element_text(size = input$axis_title_size),
        axis.text  = element_text(size = input$axis_text_size)
      )
      return(p)
      
    } else if (figure_type %in% c(
      "Heatmap (All Genes)",
      "Heatmap (Significant Genes)",
      "Custom Heatmap"
    )) {
      
      df_expr        <- mod_data()
      do_cluster_rows <- isTRUE(input$cluster_rows)
      do_cluster_cols <- isTRUE(input$cluster_cols)
      do_scale        <- isTRUE(input$do_scale)
      do_agg          <- isTRUE(input$aggregate_heatmap)
      group_col       <- input$group_col
      
      # 1) choose genes based on mode
      if (figure_type == "Heatmap (Significant Genes)") {
        genes <- find_significant_genes()$Gene
        # if user also selected some manually, intersect:
        if (!is.null(input$selected_genes) && length(input$selected_genes))
          genes <- intersect(genes, input$selected_genes)
        default_title <- "Heatmap of Significant Genes"
        
      } else if (figure_type == "Custom Heatmap") {
        genes <- input$selected_genes
        if (length(genes) == 0) {
          return(
            ggplot() +
              annotate("text", x = .5, y = .5,
                       label = "Please select ≥1 gene for Custom Heatmap") +
              theme_void()
          )
        }
        default_title <- paste("Custom Heatmap of", paste(genes, collapse = ", "))
        
      } else {
        # All Genes
        genes <- input$selected_genes
        if (length(genes) == 0) genes <- colnames(df_expr)
        default_title <- if (length(genes) == ncol(df_expr))
          "Heatmap of All Genes"
        else
          paste("Heatmap of", paste(genes, collapse = ", "))
      }
      
      # 2) optional aggregation by grouping column
      if (do_agg && group_col != "None") {
        group_vec <- factor(raw_df()[[group_col]])
        df_m      <- cbind(GroupVar = group_vec,
                           as.data.frame(df_expr)[, genes, drop = FALSE])
        df_means  <- df_m %>%
          group_by(GroupVar) %>%
          summarise(across(all_of(genes), ~ mean(.x, na.rm = TRUE)),
                    .groups = "drop")
        mat       <- as.matrix(df_means[ , -1, drop = FALSE])
        rownames(mat) <- df_means$GroupVar
        
      } else {
        mat <- as.matrix(df_expr[, genes, drop = FALSE])
        rownames(mat) <- rownames(df_expr)
      }
      
      # 3) scale if requested
      if (do_scale) mat <- scale(mat)
      
      # 4) transpose if requested
      if (isTRUE(input$transpose_heatmap)) {
        mat             <- t(mat)
        tmp             <- do_cluster_rows
        do_cluster_rows <- do_cluster_cols
        do_cluster_cols <- tmp
      }
      
      # 5) draw
      title <- if (nzchar(input$figure_title)) input$figure_title else default_title
      hm <- pheatmap(
        mat,
        cluster_rows = do_cluster_rows,
        cluster_cols = do_cluster_cols,
        main         = title,
        color        = colorRampPalette(c("navy", "white", "firebrick3"))(50),
        fontsize     = input$axis_text_size,
        fontsize_row = input$axis_text_size,
        fontsize_col = input$axis_text_size
      )
      return(hm$gtable)
    }
    
    
    
    # 3) Other Non-PCA Plots
    else {
      df_expr <- mod_data()
      req(input$selected_genes)
      selected_genes <- input$selected_genes
      multi <- length(selected_genes) > 1
      
      if(multi) {
        plot_data <- df_expr %>% 
          dplyr::select(all_of(selected_genes)) %>% 
          mutate(Treatment = factor(attr(df_expr, "treatment"))) %>%
          pivot_longer(cols = all_of(selected_genes), names_to = "Gene", values_to = "Expression")
      } else {
        plot_data <- data.frame(
          Expression = as.numeric(df_expr[[selected_genes]]),
          Treatment  = factor(attr(df_expr, "treatment"))
        )
      }
      
      default_title <- if(multi) {
        paste("Expression of", paste(selected_genes, collapse = ", "))
      } else {
        paste("Expression of", selected_genes)
      }
      final_title <- if(nchar(input$figure_title) > 0) input$figure_title else default_title
      
      add_common_theme <- function(p) {
        p + scale_fill_manual(values = get_treatment_colors()) +
          theme_minimal() +
          labs(
            title = final_title, 
            x = input$x_label, 
            y = input$y_label
          ) +
          theme(
            plot.background = element_rect(fill = input$bg_color),
            plot.title      = element_text(size = input$title_size),
            axis.title      = element_text(size = input$axis_title_size),
            axis.text       = element_text(size = input$axis_text_size)
          )
      }
      
      if (figure_type == "Boxplot") {
        if(!multi) {
          stats_data <- plot_data %>%
            group_by(Treatment) %>%
            summarise(
              mean = mean(Expression),
              sem  = sd(Expression) / sqrt(n()),
              .groups = 'drop'
            )
          model <- aov(Expression ~ Treatment, data = plot_data)
          hsd <- HSD.test(model, "Treatment")
          stats_data$letters <- hsd$groups[as.character(stats_data$Treatment), "groups"]
          has_sig_diff <- length(unique(stats_data$letters)) > 1
          
          p <- ggplot(stats_data, aes(x = Treatment, y = mean, fill = Treatment)) +
            geom_bar(stat = "identity", alpha = 0.7) +
            geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2)
          if(has_sig_diff) {
            p <- p + geom_text(
              aes(y = mean + sem + max(sem)/10, label = letters),
              size = 4, hjust = -0.5
            )
          }
        } else {
          stats_data <- plot_data %>%
            group_by(Treatment, Gene) %>%
            summarise(
              mean = mean(Expression),
              sem  = sd(Expression) / sqrt(n()),
              .groups = 'drop'
            )
          sig_list <- lapply(unique(plot_data$Gene), function(g) {
            df_g <- plot_data %>% filter(Gene == g)
            model <- aov(Expression ~ Treatment, data = df_g)
            hsd <- HSD.test(model, "Treatment")
            tibble(
              Treatment = rownames(hsd$groups), 
              letters   = hsd$groups[,"groups"], 
              Gene      = g
            )
          })
          sig_letters <- bind_rows(sig_list)
          stats_data <- left_join(stats_data, sig_letters, by = c("Treatment", "Gene"))
          
          p <- ggplot(stats_data, aes(x = Treatment, y = mean, fill = Treatment)) +
            geom_bar(stat = "identity", alpha = 0.7) +
            geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
            facet_wrap(~Gene)
          
          if(any(!is.na(stats_data$letters))) {
            p <- p + geom_text(
              aes(y = mean + sem + max(sem)/10, label = letters),
              size = 4, hjust = -0.5
            )
          }
        }
        return(add_common_theme(p))
        
      } else if (figure_type == "Boxplot_KWT") {
        if(!multi) {
          stats_data <- plot_data %>%
            group_by(Treatment) %>%
            summarise(
              median = median(Expression),
              IQR = IQR(Expression),
              .groups = 'drop'
            )
          kw <- kruskal(plot_data$Expression, plot_data$Treatment, alpha = 0.05)
          stats_data$letters <- kw$groups[as.character(stats_data$Treatment), "groups"]
          has_sig_diff <- length(unique(stats_data$letters)) > 1
          
          p <- ggplot(stats_data, aes(x = Treatment, y = median, fill = Treatment)) +
            geom_bar(stat = "identity", alpha = 0.7) +
            geom_errorbar(aes(ymin = median - IQR/2, ymax = median + IQR/2), width = 0.2)
          if(has_sig_diff) {
            p <- p + geom_text(
              aes(y = median + IQR/2 + max(IQR)/10, label = letters),
              size = 4, hjust = -0.5
            )
          }
        } else {
          stats_data <- plot_data %>%
            group_by(Treatment, Gene) %>%
            summarise(
              median = median(Expression),
              IQR = IQR(Expression),
              .groups = 'drop'
            )
          sig_list <- lapply(unique(plot_data$Gene), function(g) {
            df_g <- plot_data %>% filter(Gene == g)
            kw <- kruskal(df_g$Expression, df_g$Treatment, alpha = 0.05)
            tibble(
              Treatment = rownames(kw$groups),
              letters = kw$groups[,"groups"],
              Gene = g
            )
          })
          sig_letters <- bind_rows(sig_list)
          stats_data <- left_join(stats_data, sig_letters, by = c("Treatment", "Gene"))
          
          p <- ggplot(stats_data, aes(x = Treatment, y = median, fill = Treatment)) +
            geom_bar(stat = "identity", alpha = 0.7) +
            geom_errorbar(aes(ymin = median - IQR/2, ymax = median + IQR/2), width = 0.2) +
            facet_wrap(~Gene)
          if(any(!is.na(stats_data$letters))) {
            p <- p + geom_text(
              aes(y = median + IQR/2 + max(IQR)/10, label = letters),
              size = 4, hjust = -0.5
            )
          }
        }
        return(add_common_theme(p))
        
      } else if (figure_type == "Boxplot with Dots") {
        p <- ggplot(plot_data, aes(x = Treatment, y = Expression, fill = Treatment)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(width = 0.2, size = 2, alpha = 0.7)
        if(multi) p <- p + facet_wrap(~Gene)
        return(add_common_theme(p))
        
      } else if (figure_type == "Violin Plot") {
        p <- ggplot(plot_data, aes(x = Treatment, y = Expression, fill = Treatment)) +
          geom_violin(trim = FALSE) +
          geom_jitter(width = 0.2, size = 2, alpha = 0.7)
        if(multi) p <- p + facet_wrap(~Gene)
        return(add_common_theme(p))
        
      } else if (figure_type == "Boxplot & Violin Plot") {
        p <- ggplot(plot_data, aes(x = Treatment, y = Expression, fill = Treatment)) +
          geom_violin(trim = FALSE, alpha = 0.5) +
          geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
          geom_jitter(width = 0.2, size = 2, alpha = 0.7)
        if(multi) p <- p + facet_wrap(~Gene)
        return(add_common_theme(p))
        
      } else if (figure_type == "Outlier Boxplot") {
        outlier_data <- detect_outliers()
        if(!multi) {
          outlier_data <- outlier_data %>% filter(Gene == selected_genes)
          p <- ggplot(plot_data, aes(x = Treatment, y = Expression, fill = Treatment)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
            geom_point(data = outlier_data, aes(x = Treatment, y = Value), 
                       color = "red", size = 3)
        } else {
          outlier_data <- outlier_data %>% filter(Gene %in% selected_genes)
          p <- ggplot(plot_data, aes(x = Treatment, y = Expression, fill = Treatment)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
            geom_point(data = outlier_data, aes(x = Treatment, y = Value), 
                       color = "red", size = 3) +
            facet_wrap(~Gene)
        }
        return(add_common_theme(p))
        
      } else if (figure_type == "Custom qPCR-style Plot") {
        req(input$selected_genes)
        gene <- input$selected_genes[1]
        
        df_expr <- mod_data()
        treatments <- attr(df_expr, "treatment")
        
        plot_data <- data.frame(
          Expression = as.numeric(df_expr[[gene]]),
          Treatment  = factor(treatments)
        )
        
        plot_data$Expression <- log10(plot_data$Expression + 1e-6)
        
        # Color palette
        n_treatments <- length(unique(plot_data$Treatment))
        colors <- c("steelblue", colorRampPalette(c("#ffcccc", "#ff0000"))(n_treatments - 1))
        
        p <- ggplot(plot_data, aes(x = Treatment, y = Expression, color = Treatment)) +
          geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
          stat_summary(fun = mean, geom = "point", size = 3, shape = 18) +
          stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2) +
          scale_color_manual(values = colors) +
          scale_fill_manual(values = colors) +
          theme_minimal() +
          labs(
            title = paste("Expression of", gene),
            y = expression(Log[10]~"Expression"),
            x = "Treatment"
          ) +
          theme(
            plot.title = element_text(size = input$title_size),
            axis.title = element_text(size = input$axis_title_size),
            axis.text  = element_text(size = input$axis_text_size)
          )
        
        if (nlevels(plot_data$Treatment) > 1) {
          stat_test <- rstatix::pairwise_wilcox_test(
            data = plot_data,
            formula = Expression ~ Treatment,
            p.adjust.method = "BH"
          ) %>%
            dplyr::filter(p.adj < 0.05) %>%
            rstatix::add_xy_position(x = "Treatment", data = plot_data, step.increase = 0.1)
          
          stat_test$label <- formatC(stat_test$p.adj, format = "e", digits = 2)
          
          if (nrow(stat_test) > 0) {
            p <- p + stat_pvalue_manual(stat_test, label = "label", tip.length = 0.01)
          }
        }
        
        return(p)
        
      } else if (figure_type == "Custom qPCR-style Plot (Parametric)") {
        req(input$selected_genes)
        gene <- input$selected_genes[1]
        
        df_expr <- mod_data()
        treatments <- attr(df_expr, "treatment")
        
        plot_data <- data.frame(
          Expression = as.numeric(df_expr[[gene]]),
          Treatment  = factor(treatments)
        )
        
        plot_data$Expression <- log10(plot_data$Expression + 1e-6)
        
        n_treatments <- length(unique(plot_data$Treatment))
        colors <- c("steelblue", colorRampPalette(c("#ffcccc", "#ff0000"))(n_treatments - 1))
        
        p <- ggplot(plot_data, aes(x = Treatment, y = Expression, color = Treatment)) +
          geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
          stat_summary(fun = mean, geom = "point", size = 3, shape = 18) +
          stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2) +
          scale_color_manual(values = colors) +
          scale_fill_manual(values = colors) +
          theme_minimal() +
          labs(
            title = paste("Expression of", gene),
            y = expression(Log[10]~"Expression"),
            x = "Treatment"
          ) +
          theme(
            plot.title = element_text(size = input$title_size),
            axis.title = element_text(size = input$axis_title_size),
            axis.text  = element_text(size = input$axis_text_size)
          )
        
        if (nlevels(plot_data$Treatment) > 1) {
          stat_test <- rstatix::pairwise_t_test(
            data = plot_data,
            formula = Expression ~ Treatment,
            p.adjust.method = "BH"
          ) %>%
            dplyr::filter(p.adj < 0.05) %>%
            rstatix::add_xy_position(x = "Treatment", data = plot_data, step.increase = 0.1)
          
          stat_test$label <- formatC(stat_test$p.adj, format = "e", digits = 2)
          
          if (nrow(stat_test) > 0) {
            p <- p + stat_pvalue_manual(stat_test, label = "label", tip.length = 0.01)
          }
        }
        
        return(p)
      } else if (figure_type == "Custom qPCR-style Plot (Raw)") {
        req(input$selected_genes)
        gene <- input$selected_genes[1]
        
        df_expr <- mod_data()
        treatments <- attr(df_expr, "treatment")
        
        plot_data <- data.frame(
          Expression = as.numeric(df_expr[[gene]]),
          Treatment  = factor(treatments)
        )
        
        # No log transform
        n_treatments <- length(unique(plot_data$Treatment))
        colors <- c("steelblue", colorRampPalette(c("#ffcccc", "#ff0000"))(n_treatments - 1))
        
        p <- ggplot(plot_data, aes(x = Treatment, y = Expression, color = Treatment)) +
          geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
          stat_summary(fun = mean, geom = "point", size = 3, shape = 18) +
          stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2) +
          scale_color_manual(values = colors) +
          scale_fill_manual(values = colors) +
          theme_minimal() +
          labs(
            title = paste("Raw Expression of", gene),
            y = "Expression (Fold Change)",
            x = "Treatment"
          ) +
          theme(
            plot.title = element_text(size = input$title_size),
            axis.title = element_text(size = input$axis_title_size),
            axis.text  = element_text(size = input$axis_text_size)
          )
        
        if (nlevels(plot_data$Treatment) > 1) {
          stat_test <- rstatix::pairwise_wilcox_test(
            data = plot_data,
            formula = Expression ~ Treatment,
            p.adjust.method = "BH"
          ) %>%
            dplyr::filter(p.adj < 0.05) %>%
            rstatix::add_xy_position(x = "Treatment", data = plot_data, step.increase = 0.1)
          
          stat_test$label <- formatC(stat_test$p.adj, format = "e", digits = 2)
          
          if (nrow(stat_test) > 0) {
            p <- p + stat_pvalue_manual(stat_test, label = "label", tip.length = 0.01)
          }
        }
        
        return(p)
        
      } else if (figure_type == "Custom qPCR-style Plot (Parametric, Raw)") {
        req(input$selected_genes)
        gene <- input$selected_genes[1]
        
        df_expr <- mod_data()
        treatments <- attr(df_expr, "treatment")
        
        plot_data <- data.frame(
          Expression = as.numeric(df_expr[[gene]]),
          Treatment  = factor(treatments)
        )
        
        # No log transform
        n_treatments <- length(unique(plot_data$Treatment))
        colors <- c("steelblue", colorRampPalette(c("#ffcccc", "#ff0000"))(n_treatments - 1))
        
        p <- ggplot(plot_data, aes(x = Treatment, y = Expression, color = Treatment)) +
          geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
          stat_summary(fun = mean, geom = "point", size = 3, shape = 18) +
          stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2) +
          scale_color_manual(values = colors) +
          scale_fill_manual(values = colors) +
          theme_minimal() +
          labs(
            title = paste("Raw Expression of", gene),
            y = "Expression (Fold Change)",
            x = "Treatment"
          ) +
          theme(
            plot.title = element_text(size = input$title_size),
            axis.title = element_text(size = input$axis_title_size),
            axis.text  = element_text(size = input$axis_text_size)
          )
        
        if (nlevels(plot_data$Treatment) > 1) {
          stat_test <- rstatix::pairwise_t_test(
            data = plot_data,
            formula = Expression ~ Treatment,
            p.adjust.method = "BH"
          ) %>%
            dplyr::filter(p.adj < 0.05) %>%
            rstatix::add_xy_position(x = "Treatment", data = plot_data, step.increase = 0.1)
          
          stat_test$label <- formatC(stat_test$p.adj, format = "e", digits = 2)
          
          if (nrow(stat_test) > 0) {
            p <- p + stat_pvalue_manual(stat_test, label = "label", tip.length = 0.01)
          }
        }
        
        return(p)
      }
      
      
    }
  }
  #########################
  generate_single_gene_plot <- function(gene, show_legend = TRUE, download_mode = TRUE, data = NULL) {
    df <- if (is.null(data)) mod_data() else data
    req(mod_data())
    replace_outliers_flag(TRUE)
    df <- mod_data()
    temp <- data.frame(Expression = as.numeric(df[[gene]]),
                       Treatment = factor(attr(df, "treatment")))
    
    default_title <- paste("Expression of", gene)
    final_title <- if (nchar(input$figure_title) > 0) input$figure_title else default_title
    ft <- input$figure_type
    
    # Inside generate_single_gene_plot:
    title_size      <- if (download_mode) 10 else input$title_size
    axis_title_size <- if (download_mode) 9 else input$axis_title_size
    axis_text_size  <- if (download_mode) 7 else input$axis_text_size
    dot_size        <- if (download_mode) 0.8 else 2
    
    
    p <- NULL
    
    if (ft == "Boxplot with Dots") {
      p <- ggplot(temp, aes(x = Treatment, y = Expression, fill = Treatment)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(aes(color = Treatment), width = 0.2, size = dot_size, alpha = 0.7) +
        scale_fill_manual(values = get_treatment_colors()) +
        scale_color_manual(values = get_treatment_colors()) +
        labs(title = final_title, x = input$x_label, y = input$y_label)
    } else if (ft == "Violin Plot") {
      p <- ggplot(temp, aes(x = Treatment, y = Expression, fill = Treatment)) +
        geom_violin(trim = FALSE) +
        geom_jitter(width = 0.2, size = dot_size, alpha = 0.7) +
        labs(title = final_title, x = input$x_label, y = input$y_label)
    } else if (ft == "Boxplot & Violin Plot") {
      p <- ggplot(temp, aes(x = Treatment, y = Expression, fill = Treatment)) +
        geom_violin(trim = FALSE, alpha = 0.5) +
        geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2, size = dot_size, alpha = 0.7) +
        labs(title = final_title, x = input$x_label, y = input$y_label)
    } else {
      p <- ggplot(temp, aes(x = Treatment, y = Expression, fill = Treatment)) +
        geom_boxplot() +
        labs(title = final_title, x = input$x_label, y = input$y_label)
    }
    
    p + theme_minimal() +
      theme(
        plot.title = element_text(size = title_size, margin = margin(b = 2)),
        axis.title = element_text(size = axis_title_size),
        axis.text  = element_text(size = axis_text_size),
        plot.background = element_rect(fill = input$bg_color),
        legend.position = if (show_legend) "right" else "none"
      )
  }
  
  
  
  
  ###### Render / Download ######
  output$main_plot <- renderPlot({
    p <- generate_plot()
    if (inherits(p, "gtable")) {
      grid::grid.newpage()
      grid::grid.draw(p)
    } else {
      print(p)
    }
  })  # PDF Download of All Gene Figures (2 per page)
  
  # PDF Download of All Gene Figures (2 per page)
  
  output$download_all_figures_pdf <- downloadHandler(
    filename = function(){
      paste0(gsub(" ", "_", input$figure_type), "_all_figures_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content = function(file){
      pdf(file, paper = "a4", width = 8.27, height = 11.69)
      
      all_genes <- colnames(mod_data())
      n <- length(all_genes)
      plots_per_page <- 12  # 3x4 layout
      
      ## ----- PAGE 1: Legend Only -----
      temp_data <- data.frame(
        Expression = rnorm(10), 
        Treatment  = factor(attr(mod_data(), "treatment")[1:10])
      )
      p_legend <- ggplot(temp_data, aes(x = Treatment, y = Expression, fill = Treatment)) +
        geom_boxplot() +
        scale_fill_manual(values = get_treatment_colors()) +
        theme_void() +
        theme(legend.position = "right")
      
      # Extract legend using cowplot
      legend_only <- cowplot::get_legend(p_legend)
      cowplot::plot_grid(legend_only)  # Add legend as full first page
      
      ## ----- PLOT PAGES -----
      for(i in seq(1, n, by = plots_per_page)) {
        end <- min(i + plots_per_page - 1, n)
        gene_subset <- all_genes[i:end]
        
        plots <- lapply(gene_subset, function(gene) {
          tryCatch({
            generate_single_gene_plot(gene, show_legend = FALSE)
          }, error = function(e) {
            message(sprintf("Error in gene %s: %s", gene, e$message))
            ggplot() + theme_void() +
              annotate("text", x = 0.5, y = 0.5, label = paste("Error in gene", gene), size = 5)
          })
        })
        
        gridExtra::grid.arrange(grobs = plots, ncol = 3, nrow = 4, top=NULL)
      }
      
      dev.off()
    }
  )
  
  output$download_selected_genes_pdf <- downloadHandler(
    filename = function(){
      paste0(gsub(" ", "_", input$figure_type), "_selected_genes_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf")
    },
    content = function(file){
      req(input$selected_genes)
      selected_genes <- input$selected_genes
      df <- mod_data()
      replace_outliers_flag(TRUE)
      
      pdf(file, paper = "a4", width = 8.27, height = 11.69)
      plots_per_page <- 12  # Layout: 3 x 4
      total <- length(selected_genes)
      
      for(i in seq(1, total, by = plots_per_page)) {
        genes_batch <- selected_genes[i:min(i + plots_per_page - 1, total)]
        
        plots <- lapply(genes_batch, function(gene) {
          tryCatch({
            generate_single_gene_plot(gene, show_legend = FALSE, download_mode = TRUE, data = df)
          }, error = function(e) {
            ggplot() + theme_void() +
              annotate("text", x = 0.5, y = 0.5, label = paste("Error:", gene), size = 5)
          })
        })
        
        gridExtra::grid.arrange(grobs = plots, ncol = 3, nrow = 4)
      }
      
      dev.off()
    }
  )
  
  
  output$download_excel_stats <- downloadHandler(
    filename = function() {
      paste0("expression_statistics_", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      library(openxlsx)
      wb <- createWorkbook()
      
      addWorksheet(wb, "Normality")
      writeData(wb, "Normality", normality_data())
      
      addWorksheet(wb, "Analysis_Type")
      writeData(wb, "Analysis_Type", analysis_type_data())
      
      
      # Just ANOVA p-values
      addWorksheet(wb, "ANOVA_KW_Pvalues")
      apv <- all_genes_pvalues()
      if (is.null(apv) || nrow(apv) == 0) apv <- data.frame(Gene=character(), ANOVA_p=numeric(), KWT_p=numeric())
      writeData(wb, "ANOVA_KW_Pvalues", apv)
      
      addWorksheet(wb, "Pairwise_All_Tests")
      pw <- pairwise_results()
      if (is.null(pw) || nrow(pw) == 0)
        pw <- data.frame(Gene=character(), Test=character(), group1=character(), group2=character(), p=numeric(), p.adj=numeric(), p.adj.signif=character())
      writeData(wb, "Pairwise_All_Tests", pw)
      
      addWorksheet(wb, "Pairwise_Tukey")
      pw_tukey <- pairwise_tukey_results()
      if (is.null(pw_tukey) || nrow(pw_tukey) == 0)
        pw_tukey <- data.frame(Gene=character(), Test=character(), group1=character(), group2=character(), p=numeric(), p.adj=numeric(), p.adj.signif=character())
      writeData(wb, "Pairwise_Tukey", pw_tukey)
      
      addWorksheet(wb, "All_Genes_Pvalues_FDR")
      apv_fdr <- all_genes_pvalues_fdr()
      if (is.null(apv_fdr) || nrow(apv_fdr) == 0)
        apv_fdr <- data.frame(Gene=character(), p_value=numeric(), method=character(),
                              adj_p_BH=numeric(), adj_p_Bonferroni=numeric())
      writeData(wb, "All_Genes_Pvalues_FDR", apv_fdr)
      
      addWorksheet(wb, "Significant_Genes")
      sig_genes <- find_significant_genes()
      writeData(wb, "Significant_Genes", sig_genes)
      
      
      saveWorkbook(wb, file, overwrite = TRUE)
    },
    contentType = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
  )
  
  
  output$pca_summary <- renderPrint({ req(pca_result()); summary(pca_result()) })
  
  output$data_preview <- renderDT({ req(mod_data()); datatable(mod_data(), options = list(scrollX = TRUE, pageLength = 5)) })
  
  output$sig_genes_table <- renderDT({
    req(mod_data())
    sig_genes_info <- find_significant_genes()
    sig_genes_info$P_Value <- format(sig_genes_info$P_Value, scientific = TRUE, digits = 3)
    datatable(sig_genes_info, options = list(pageLength = 5), caption = "Significant Genes (p < 0.05)")
  })
  
  output$download_figure <- downloadHandler(
    filename = function() { paste0(input$figure_type, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png") },
    content = function(file) {
      png(file, width = 2400, height = 2400, res = 300)
      grid.newpage()
      p <- generate_plot()
      if (inherits(p, "gtable")) {
        grid::grid.draw(p)  # ✅ This correctly renders heatmap objects
      } else {
        print(p)            # ✅ This works for ggplot2-based plots
      }
      dev.off()
    }
  )
  
  
  observeEvent(input$show_about, {
    showModal(modalDialog(
      title = "About: Gene Expression Analysis Shiny App",
      size = "l",
      easyClose = TRUE,
      footer = modalButton("Close"),
      tagList(
        tags$h4("What is this?"),
        tags$p("This program is an interactive platform for gene expression data analysis.
        Aslo can be used for any different type of data that have sample ID in the 1st column and Treatment group in the 2nd column and 3rd, 4th and so on different measurment for numeric data.
      You can perform outlier detection, imputation, statistical testing (ANOVA, Kruskal-Wallis, post-hoc tests), PCA, and generate publication-ready plots and summaries."),
        
        tags$h4("Input File Format"),
        tags$p("Upload your data as a CSV or TSV file with the following structure:"),
        tags$ul(
          tags$li("First column: Sample ID"),
          tags$li("Second column: Treatment or group (categorical variable)"),
          tags$li("Remaining columns: Gene expression values (one column per gene, numeric)")
        ),
        tags$p("Missing values are allowed and will be imputed as specified. Treatments should be named consistently."),
        
        tags$h4("Analyses Available"),
        tags$ul(
          tags$li("Outlier detection and replacement (per group or globally)"),
          tags$li("Imputation of missing values (mean, median, or zero)"),
          tags$li("Normality testing (Shapiro-Wilk)"),
          tags$li("Significant gene identification (ANOVA, Kruskal-Wallis, automatic based on normality, with optional FDR correction)"),
          tags$li("Post-hoc pairwise comparisons (Tukey, Dunn, t-test, Wilcoxon)"),
          tags$li("Principal Component Analysis (PCA) and biplots"),
          tags$li("Heatmaps (all genes or significant genes only)"),
          tags$li("Boxplots, violin plots, and custom qPCR-style plots")
        ),
        
        tags$h4("Outputs"),
        tags$ul(
          tags$li("Interactive plots in the browser"),
          tags$li("Downloadable plots (PNG, PDF for all genes)"),
          tags$li("Excel summary including normality, statistical test results, and significant gene lists"),
          tags$li("Tables of outliers, imputed values, and significant genes")
        ),
        
        tags$h4("Contact"),
        tags$p(HTML('Developed by Janan Gawra. For help, contact via <a href="https://www.linkedin.com/in/janangawra/" target="_blank">LinkedIn</a>.'))
      )
    ))
  })
  
}

shinyApp(ui, server)
