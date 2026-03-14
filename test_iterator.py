import numpy as np
import array_api_compat.numpy as xp

# 1. Simulate an input table of integers (like raw counts)
matrix = xp.asarray([[10, 20], [30, 40]], dtype=xp.int64)
print(f"Input Matrix Dtype: {matrix.dtype}")

# 2. Simulate what happened without the fix: 
# Initializing pval_mat using the matrix's integer dtype
pval_mat_int = xp.empty((2, 2), dtype=matrix.dtype)

# 3. A typical float p-value from a statistical test
test_p_value = 0.045

# 4. Try to assign the float to the integer matrix
pval_mat_int[0, 1] = test_p_value

print(f"Assigned p-value: {test_p_value}")
print(f"Value stored in pval_mat (int dtype): {pval_mat_int[0, 1]}")

if pval_mat_int[0, 1] == 0:
    print("Verification: The float was truncated to 0!")
