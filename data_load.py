import pickle
import numpy as np
import pandas as pd

def load_pickle(file_path):
    with open(file_path, 'rb') as f:
        data = pickle.load(f)
    if isinstance(data, np.ndarray):
        return data.tolist()
    return data
