"""MMS case names and fields, independent of the symbolic input generator."""

EVOLVED = ("Ni", "Nj", "NVi", "NVj", "Pe", "Vort", "Eion")
RECOVERED = ("Pi", "pi", "pj", "Tion", "phi", "psi", "chi")
PARALLEL_VARIANTS = ("parallel_current", "parallel_current_wall")
VARIANTS = ("baseline", "skew", "viscosity", "curvature", *PARALLEL_VARIANTS)
