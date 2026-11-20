"""Runtime fields for the single 3D slab MMS, independent of SymPy."""

EVOLVED = ("Ni", "Nj", "NVi", "NVj", "Pe", "Vort", "Eion")
RECOVERED = ("Pi", "pi", "pj", "Tion", "phi", "psi", "chi")
PARALLEL_VARIANTS = ("parallel_current_wall",)
VARIANTS = PARALLEL_VARIANTS
