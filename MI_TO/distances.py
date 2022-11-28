"""
Module to create custom distance function among cell AF profiles.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import numba
from joblib import Parallel, delayed, cpu_count

# 
