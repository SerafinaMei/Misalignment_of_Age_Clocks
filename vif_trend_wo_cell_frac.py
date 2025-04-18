import numpy as np
import pandas as pd
from statsmodels.stats.outliers_influence import variance_inflation_factor


def vif_with_features(df, cpg):
    X = df[df.columns.intersection(cpg)].copy()
    X.insert(0, 'Intercept', 1)
    vif_data = pd.DataFrame()
    vif_data["Feature"] = X.columns
    vif_data["VIF"] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
    return vif_data

def compute_same_trend(correlation_df, model_df, feature_col, correlation_col, model_index_col, model_value_col):
    model_index = model_df.index if model_index_col == "Index" else model_df[model_index_col]
    correlation_values = correlation_df.set_index(feature_col).loc[
        correlation_df.set_index(feature_col).index.intersection(model_index)
    ][correlation_col]
    model_filtered = model_df[~model_df.index.str.lower().str.contains("intercept", na=False)] if model_index_col == "Index" else \
                     model_df[~model_df[model_index_col].str.lower().str.contains("intercept", na=False)]
    model_values = model_filtered[model_value_col] if model_index_col == "Index" else model_filtered.set_index(model_index_col)[model_value_col]
    correlation_values = pd.to_numeric(correlation_values, errors='coerce')
    model_values = pd.to_numeric(model_values, errors='coerce')
    return ((correlation_values * model_values) > 0).sum()
