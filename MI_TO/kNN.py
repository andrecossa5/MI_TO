"""
Module for kNN utils. Import kNN from Cellula by default, but implements utilities for creating 
kNN masked affinity matrices from full affinity matrices. 
"""

import numpy as np
import pandas as pd
import scanpy as sc
import anndata

from Cellula.preprocessing._neighbors import kNN_graph

