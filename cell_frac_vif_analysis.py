import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor
import matplotlib.pyplot as plt

def compute_category_wise_vif(cpg_matrix, cell_fractions):
    results = []
    X_cpgs = cpg_matrix.copy()
    X_cpgs.insert(0, 'Intercept', 1)
    vif_baseline = [variance_inflation_factor(X_cpgs.values, i) for i in range(X_cpgs.shape[1])]

    results.append({
        "CellType": "Baseline (CpGs Only)",
        "Mean_VIF": np.mean(vif_baseline[1:]),  
        "Max_VIF": np.max(vif_baseline[1:]),
        "VIF_List": vif_baseline[1:],
        "VIF_Values": vif_baseline[1:],  
        "CpGs": X_cpgs.columns[1:]
    })
    
    for cell_type in cell_fractions.columns:
        X_combined = cpg_matrix.copy()
        X_combined[cell_type] = cell_fractions[cell_type]
        X_combined.insert(0, 'Intercept', 1)
        
        vif_values = [variance_inflation_factor(X_combined.values, i) for i in range(X_combined.shape[1])]
        mean_vif = np.mean(vif_values[1:]) 
        max_vif = np.max(vif_values[1:])

        results.append({
            "CellType": cell_type,
            "Mean_VIF": mean_vif,
            "Max_VIF": max_vif,
            "VIF_List": vif_values[1:],
            "VIF_Values": vif_values[1:],  
            "CpGs": X_combined.columns[1:]
        })
    
    return pd.DataFrame(results)

def run_analysis(cpg_data, cell_fraction_data):
    vif_results = compute_category_wise_vif(cpg_data, cell_fraction_data)
    return vif_results
