# The default type to use for the quantum circuit library. Does not apply to cplx.
alias DEFAULT_TYPE = DType.float64
# The default tolerance used for `is_close` comparisons, matrix invertability checks, etc.
alias DEFAULT_TOL = 2.22044604925e-16  # Machine epsilon for 64 bit floats
# The default zero threshold for sparse matrices and statevectors. Values with norm less 
# than or equal to the zero threshold are considered zero.
alias DEFAULT_ZERO_THRESHOD = 1e-20
