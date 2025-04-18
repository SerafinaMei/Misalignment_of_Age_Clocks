library(ggplot2)
library(car)
library(dplyr)
library(tidyr)
library(stringr)

compute_vif_increase <- function(hannum_selected, cell_fractions, clock_name) {
    initial_vif <- vif_with_features(py$pd$DataFrame(hannum_selected))
    initial_vif <- as.data.frame(initial_vif)
    colnames(initial_vif) <- c("Feature", "VIF")
    initial_vif$VIF[!is.finite(initial_vif$VIF)] <- NA
    
    results <- data.frame()
    
    for (cell_type in colnames(cell_fractions)) {
        combined_data <- cbind(hannum_selected, cell_fractions[, cell_type, drop=FALSE])
        vif_results <- vif_with_features(py$pd$DataFrame(combined_data))
        vif_results <- as.data.frame(vif_results)
        colnames(vif_results) <- c("Feature", "VIF")
        vif_results$VIF[!is.finite(vif_results$VIF)] <- NA
        
        merged <- merge(initial_vif, vif_results, by = "Feature", suffixes = c("_initial", "_new"))
        merged$Increase <- merged$VIF_new - merged$VIF_initial
        
        total_cpgs <- nrow(merged)
        cpgs_affected <- sum(merged$Increase > 0, na.rm = TRUE)
        pct_any_increase <- cpgs_affected / total_cpgs
        
        if (cpgs_affected > 0) {
            pct_0_0.25 <- sum(merged$Increase > 0 & merged$Increase <= 0.25, na.rm = TRUE) / cpgs_affected
            pct_0.25_0.5 <- sum(merged$Increase > 0.25 & merged$Increase <= 0.5, na.rm = TRUE) / cpgs_affected
            pct_0.5_1 <- sum(merged$Increase > 0.5 & merged$Increase <= 1, na.rm = TRUE) / cpgs_affected
            pct_1_2 <- sum(merged$Increase > 1 & merged$Increase <= 2, na.rm = TRUE) / cpgs_affected
            pct_above_2 <- sum(merged$Increase > 2, na.rm = TRUE) / cpgs_affected
        } else {
            pct_0_0.25 <- 0
            pct_0.25_0.5 <- 0
            pct_0.5_1 <- 0
            pct_1_2 <- 0
            pct_above_2 <- 0
        }
        
        results <- rbind(results, data.frame(Cell_Type = cell_type, Clock = clock_name,
                                             Total_CpGs = total_cpgs,
                                             CpGs_Affected = cpgs_affected,
                                             Total_Pct_Any_Increase = pct_any_increase,
                                             Pct_0_0.25 = pct_0_0.25,
                                             Pct_0.25_0.5 = pct_0.25_0.5,
                                             Pct_0.5_1 = pct_0.5_1,
                                             Pct_1_2 = pct_1_2,
                                             Pct_Above_2 = pct_above_2))
    }
    
    print(results[, c("Cell_Type", "Total_CpGs", "CpGs_Affected")])
    
    results_long <- pivot_longer(results, cols = c("Pct_0_0.25", "Pct_0.25_0.5", "Pct_0.5_1", "Pct_1_2", "Pct_Above_2"),
                                 names_to = "Threshold", values_to = "Proportion")
    
    ggplot(results_long, aes(x = Cell_Type, y = Total_Pct_Any_Increase * Proportion, fill = Threshold)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_y_continuous(labels = scales::percent) +
        theme_minimal() +
        labs(title = paste("Proportion of CpGs with Increased VIF -", clock_name),
             x = "Cell Type", y = "Total % of CpGs Affected", fill = "VIF Increase") +
        coord_flip()
}


plot_vif_histograms <- function(initial_vif, vif_results, clock_name) {
    merged <- merge(initial_vif, vif_results, by = "Feature", suffixes = c("_initial", "_new"))
    merged$Increase <- merged$VIF_new - merged$VIF_initial
    
    ggplot(merged, aes(x = Increase)) +
        geom_histogram(bins = 30, fill = "blue", alpha = 0.6) +
        facet_wrap(~ Feature, scales = "free_y") +
        theme_minimal() +
        labs(title = paste("Histogram of VIF Increases -", clock_name),
             x = "VIF Increase", y = "Count")
}
