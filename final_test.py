import numpy as np
import array_api_compat.numpy as xp

print("--- 1. Testing Iterator Ingestion ---")
percentiles_iter = iter([25.0, 50.0])
_, res_iter = xp.asarray(percentiles_iter), None # Simulating what ingest_array does
# Note: In NumPy, asarray(iterator) creates an object array
res_iter = xp.asarray(percentiles_iter)
print(f"Result Dtype: {res_iter.dtype}")
try:
    print(f"Can we compare? {res_iter < 0}")
except Exception as e:
    print(f"Comparison Failed: {e}")

print("\n--- 2. Testing Float Truncation ---")
# Input matrix of integers (counts)
matrix = xp.asarray([10, 20], dtype=xp.int64)
# Initializing pval_mat using the integer dtype (what happens without your fix)
pval_mat = xp.empty((2,), dtype=matrix.dtype)
p_val = 0.045
pval_mat[0] = p_val
print(f"Original p-value: {p_val}")
print(f"Stored in Integer Matrix: {pval_mat[0]}")
if pval_mat[0] == 0:
    print("BUG CONFIRMED: Value truncated to 0!")
