# The default type
alias DEFAULT_TYPE = DType.float64
# The default tolerance used for `is_close` comparisons, matrix invertability checks, etc.
alias DEFAULT_TOL = 5e-14
# The default zero threshold for sparse matrices and statevectors. Values with norm less 
# than or equal to the zero threshold are considered zero.
alias DEFAULT_ZERO_THRESHOD = 1e-20

# IEEE 754
# Machine epsilon for Float64
alias EPS64 = 2.220446049250313e-16
# Smallest positive normalized Float64
alias REALMIN64 = 2.2250738585072014e-308
# Largest positive Float64
alias REALMAX64 = 1.7976931348623157e308
# Machine epsilon for Float32
alias EPS32 = 1.1920929e-07
# Smallest positive normalized Float32
alias REALMIN32 = 1.1754944e-38
# Largest positive Float32
alias REALMAX32 = 3.4028235e+38
# Machine epsilon for Float16
alias EPS16 = 0.00097656
# Smallest positive normalized Float16
alias REALMIN16 = 6.10352e-05
# Largest positive Float16
alias REALMAX16 = 65504.0