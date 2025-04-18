# library(EpiDISH)
# library(readxl)
# library(reticulate)
# library(ggplot2)
# library(pheatmap)
# library(tools)
# library(car) 
# library(caret)


load_data <- function() {
    source_python("Data_load.py")
    hannum_data <- py$load_pickle("hannum_data.pkl")
    arth_nor <- py$load_pickle("nor_arthritis.pickle")
    arth_dis <- py$load_pickle("dis_arthritis.pickle")
    return(list(hannum = as.data.frame(hannum_data),
                arth_nor = as.data.frame(arth_nor),
                arth_dis = as.data.frame(arth_dis)))
}
get_normalized_weights <- function(model_df, common_cpgs, clock_name = NULL) {
    model_weights_raw <- extract_model_weights(model_df)

    model_weights <- model_weights_raw[common_cpgs]
    model_weights <- model_weights[!is.na(model_weights)]

    if (length(model_weights) == 0 || all(model_weights == 0)) {
        if (!is.null(clock_name)) {
        warning(paste("no valid model weights found for", clock_name))
        }
        return(NULL)
    }

    model_weights <- model_weights / max(abs(model_weights), na.rm = TRUE)
    return(model_weights)
}
extract_model_weights <- function(model_df) {
    if (is.numeric(model_df)) {
        return(model_df)
    }

    possible_cols <- c("coef", "coefficient", "coefficients", "weight", "weights","CoefficientTraining")
    col_match <- which(tolower(colnames(model_df)) %in% possible_cols)

    if (length(col_match) == 0) {
        stop("no weight column found in model")
    }

    weight_col <- model_df[[col_match[1]]]
    names(weight_col) <- model_df[[1]]  # assumes first col is cpg ID
    return(weight_col)
}



detect_columns <- function(df) {
    cpg_id_column <- NULL
    
    if (any(grepl("^cg\\d+", rownames(df), ignore.case=TRUE))) {
        cpg_id_column <- "rownames"
    } else {
        for (col in colnames(df)) {
            if (any(grepl("^cg\\d+", df[[col]], ignore.case=TRUE))) {
                cpg_id_column <- col
                break
            }
        }
    }
    
    if (is.null(cpg_id_column)) {
        stop("no CpG column detected in dataset")
    }
    
    return(df[[cpg_id_column]])
}

compute_correlation <- function(hannum_data, arth_nor, arth_dis, model_cpgs) {
    
    hannum_selected <- hannum_data[, model_cpgs, drop=FALSE]
    arth_nor_selected <- arth_nor[, model_cpgs, drop=FALSE]
    arth_dis_selected <- arth_dis[, model_cpgs, drop=FALSE]
    

    hannum_selected <- hannum_selected[, apply(hannum_selected, 2, var, na.rm = TRUE) != 0]
    arth_nor_selected <- arth_nor_selected[, apply(arth_nor_selected, 2, var, na.rm = TRUE) != 0]
    arth_dis_selected <- arth_dis_selected[, apply(arth_dis_selected, 2, var, na.rm = TRUE) != 0]

    data(cent12CT450k.m)
    
    hannum_cell_fractions <- epidish(beta.m = t(hannum_data), ref.m = cent12CT450k.m, method = "RPC")
    
    arth_nor_cell_fractions <- epidish(beta.m = t(arth_nor), ref.m = cent12CT450k.m, method = "RPC")
    
    arth_dis_cell_fractions <- epidish(beta.m = t(arth_dis), ref.m = cent12CT450k.m, method = "RPC")

    
    cor_hannum <- cor(hannum_selected, hannum_cell_fractions$estF, use="complete.obs")
    cor_arth_nor <- cor(arth_nor_selected, arth_nor_cell_fractions$estF, use="complete.obs")
    cor_arth_dis <- cor(arth_dis_selected, arth_dis_cell_fractions$estF, use="complete.obs")
    

    return(list(hannum = cor_hannum,
                arth_nor = cor_arth_nor,
                arth_dis = cor_arth_dis))
}

library(pheatmap)
library(RColorBrewer)

plot_correlation <- function(cor_matrix, title) {
    output_file <- paste0("results/", gsub(" ", "_", title), "_heatmap.png")
    min_cor <- min(cor_matrix, na.rm = TRUE)
    max_cor <- max(cor_matrix, na.rm = TRUE)
    buffer <- (max_cor - min_cor) * 0.05
    adjusted_min <- min_cor - buffer
    adjusted_max <- max_cor + buffer

    heatmap_colors <- colorRampPalette(c("#2b8cbe", "#a6bddb", "#ffffbf", "#fdae61", "#d73027"))(50)
    pheatmap(t(cor_matrix),
         color = heatmap_colors,
         main = title,
         cluster_rows = FALSE, 
         cluster_cols = TRUE,  
         fontsize_row = 11,
         fontsize_col = 10,
         border_color = "grey", 
         cellwidth = 8,
         cellheight = 12,
         legend = TRUE,
         breaks = seq(-0.6, 0.6, length.out = 51), 
         filename = output_file) 
}




compute_vif <- function(data, model_cpgs, cell_fractions) {
    library(car)
    library(caret) 

    common_cpgs <- intersect(colnames(data), model_cpgs)

    if (length(common_cpgs) == 0) {
        stop("no common CpGs found")
    }

    print(paste("total CpGs in dataset:", length(colnames(data))))
    print(paste("total Model CpGs:", length(model_cpgs)))
    print(paste("common CpGs for VIF:", length(common_cpgs)))
    X <- data[, common_cpgs, drop = FALSE]
    print(paste("CpGs before removing zero variance:", ncol(X)))
    X <- X[, apply(X, 2, var, na.rm = TRUE) != 0, drop = FALSE]
    print(paste("CpGs after removing zero variance:", ncol(X)))

    X <- cbind(X, cell_fractions)
    X$Intercept <- 1  
    X_df <- as.data.frame(X)

    print("🔹 Final feature set (CpGs + Cell Fractions + Intercept):")
    print(colnames(X_df))

    print("checking for collinearity issues...")
    linear_combos <- findLinearCombos(X_df)
    
    if (!is.null(linear_combos$remove)) {
        print(paste("removing", length(linear_combos$remove), "collinear features:"))
        print(colnames(X_df)[linear_combos$remove])
        X_df <- X_df[, -linear_combos$remove, drop = FALSE]
    } else {
        print("no perfect collinearity detected.")
    }

    print("computing VIF values...")

    vif_values <- tryCatch({
        sapply(colnames(X_df), function(feature) {
            safe_feature <- paste0("`", feature, "`") 
            formula <- as.formula(paste(safe_feature, "~ ."))
            lm_model <- lm(formula, data = X_df)

            car::vif(lm_model)
        })
    }, error = function(e) {
        warning("VIF computation failed: ", e$message)
        return(NULL)
    })

    if (is.null(vif_values)) {
        stop("VIF computation failed")
    }

    print(paste("features in X_df (CpGs + Cell Fractions):", ncol(X_df)))
    print(paste("VIF Values Computed:", length(vif_values)))

    # Ensure feature count matches
    if (length(vif_values) != ncol(X_df)) {
        warning("feature count and VIF values do not match. some features might been dropped due to collinear")
    }

    # Step 8: Create DataFrame
    vif_data <- data.frame(Feature = names(vif_values), VIF = vif_values)

    return(vif_data)
}




# Plot VIF pie chart
plot_vif_pie <- function(vif_data, clock_name) {
    output_file <- paste0("results/", clock_name, "_VIF_pie.png")
    vif_data$Category <- cut(vif_data$VIF, 
                             breaks = c(0, 5, 10, 20, Inf), 
                             labels = c("VIF < 5", "VIF 5-10", "VIF 10-20", "VIF > 20"),
                             right = FALSE)
    
    vif_counts <- as.data.frame(table(vif_data$Category))
    colnames(vif_counts) <- c("Category", "Count")

    pie_chart <- ggplot(vif_counts, aes(x = "", y = Count, fill = Category)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar(theta = "y") +
        theme_minimal() +
        scale_fill_manual(values = c("VIF < 5" = "#2b8cbe", 
                                     "VIF 5-10" = "#a6bddb", 
                                     "VIF 10-20" = "#fdae61", 
                                     "VIF > 20" = "#d73027")) +
        labs(title = paste("VIF Proportion -", clock_name), fill = "VIF Range") +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank())

    ggsave(output_file, pie_chart, width = 6, height = 6)
    
    print(paste("VIF Pie saved:", output_file))
    
    return(pie_chart)
}





plot_model_level_summary <- function(summary_df) {
  # summary_df must have columns: model_name, percent_vif5, percent_disagree, MAE or MSE
  p1 <- ggplot(summary_df, aes(x = percent_vif5, y = percent_disagree, label = model_name)) +
    geom_point(size = 3, color = "darkblue") +
    geom_text(vjust = 1.5, hjust = 1, size = 3) +
    theme_minimal() +
    labs(title = "Model Scatter: VIF vs Feature Disagreement",
         x = "% Features with VIF > 5",
         y = "% Features with Weight Sign Disagreement")

  ggsave("results/model_vif_vs_disagreement.png", p1, width = 6, height = 5)

  return(p1)
}


plot_feature_quadrant <- function(agreement_df, clock_name, save_path = "trend_disagree_scatterplots") {
    quad_counts <- list(
    top_left     = sum(agreement_df$Univ > 0 & agreement_df$Model < 0),
    bottom_left  = sum(agreement_df$Univ < 0 & agreement_df$Model < 0),
    top_right    = sum(agreement_df$Univ > 0 & agreement_df$Model > 0),
    bottom_right = sum(agreement_df$Univ < 0 & agreement_df$Model > 0)
    )
 
    # Axis ranges
    xlim_range <- c(-1.1, 1.1)
    ylim_range <- c(-1.1, 1.1)

    # Plot
    p <- ggplot(agreement_df, aes(x = Model, y = Univ, color = Disagree)) +
        geom_point(size = 1.8, alpha = 0.8, shape = 16) +
        scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2")) +
        geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
        geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
        scale_x_continuous(breaks = seq(-1, 1, 0.5), limits = xlim_range) +
        scale_y_continuous(breaks = seq(-1, 1, 0.5), limits = ylim_range) +
        annotate("text", x = -0.9, y =  0.95, label = quad_counts$top_left, color = "#D55E00", size = 4.2, fontface = "bold") +
        annotate("text", x = -0.9, y = -0.95, label = quad_counts$bottom_left, color = "#0072B2", size = 4.2, fontface = "bold") +
        annotate("text", x =  0.9, y =  0.95, label = quad_counts$top_right, color = "#0072B2", size = 4.2, fontface = "bold") +
        annotate("text", x =  0.9, y = -0.95, label = quad_counts$bottom_right, color = "#D55E00", size = 4.2, fontface = "bold") +
        labs(x = "Model Assigned Coefficient", y = "Pearson Correlation with Age") +
        theme_minimal(base_size = 11) +
        theme(
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 9),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.3),
        legend.position = "none",
        plot.margin = margin(6, 6, 6, 6)
        )

    # Save path
    file_name <- file.path(save_path, paste0(clock_name, "_quadrant_plot.png"))
    ggsave(file_name, p, width = 3.2, height = 2.4, dpi = 600)

    message(paste("✅ Saved quadrant plot for", clock_name, "→", file_name))
    return(p)
}
