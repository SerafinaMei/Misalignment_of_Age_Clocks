import re
import pandas as pd
import numpy as np

def detect_columns(df):
    """
    Automatically detects the column names for CpG IDs and Coefficients (or weights).
    If the CpG ID column is already in the index, return "Index" instead.
    """
    cpg_id_column = None
    coefficient_column = None

    if df.index.astype(str).str.match(r"^(?!.*intercept).*cg\d+", case=False, na=False).sum() > len(df) * 0.8:
        cpg_id_column = "Index"

    if cpg_id_column is None:
        for col in df.columns[:2]:
            if df[col].astype(str).str.match(r"^(?!.*intercept).*cg\d+", case=False, na=False).sum() > 0:
                cpg_id_column = col
                break

    for col in df.columns:
        if re.search(r"coef(ficient)?|weight(s)?", col, re.IGNORECASE):
            coefficient_column = col
            break

    return cpg_id_column, coefficient_column

def inverse_F(x, adult_age=20):
    """
    Computes the inverse function transformation used for DNAmAge calculation.
    """
    return np.where(
        x < 0,
        (1 + adult_age) * np.exp(x) - 1,
        (1 + adult_age) * x + adult_age
    )

def count_minus_intercept(filtered_vif):
    """
    Counts rows excluding those containing 'intercept'.
    """
    filtered_vif = filtered_vif[~filtered_vif.apply(lambda row: row.astype(str).str.contains("intercept", case=False, na=False)).any(axis=1)]
    return filtered_vif.shape[0]
