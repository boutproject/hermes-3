#!/usr/bin/env python3
from job_functions import run_neutral_mixed_manufactured_solutions_test
from perpendicular_laplacian import x, y, z, div_a_grad_perp_f_symbolic
from perpendicular_laplacian import div_par_k_grad_par_f_symbolic
from perpendicular_laplacian import div_par_f_symbolic
from perpendicular_laplacian import grad_par_f_symbolic
from sympy import sin, cos

# specify symbolic inputs
# contravariant metric coeffs
g11 = 1.1 + 0.16*x*cos(y)
g22 = 0.9 + 0.09*x*cos(y)
g33 = 1.2 + 0.2*x*cos(y)
g12 = 0.0
g23 = 0.5 + 0.15*x*cos(y)
g13 = 0.0
# f and a
f = (x**2)*sin(y)*sin(z)
a = 1.0 + 0.1*x**3*sin(2*y)*sin(2*z)
# Div . ( a Grad_perp f )
div_a_grad_perp_f = div_a_grad_perp_f_symbolic(g11, g12, g13, g22, g23, g33, a, f)
# Div . ( a Grad_par f)
div_par_k_grad_par_f = div_par_k_grad_par_f_symbolic(g11, g12, g13, g22, g23, g33, a, f)
# Div . ( vec(b) f)
div_par_f = div_par_f_symbolic(g11, g12, g13, g22, g23, g33, f)
# vec(b) . Grad f
grad_par_f = grad_par_f_symbolic(g11, g12, g13, g22, g23, g33, f)

test_input = {
    "ntest" : 1, 
    "ngrid" : 20,
    "a_string": str(a),
    "f_string": str(f),
    "g11_string": str(g11),
    "g22_string": str(g22),
    "g33_string": str(g33),
    "g12_string": str(g12),
    "g13_string": str(g13),
    "g23_string": str(g23),
    "test_dir" : "neutral_mixed",
    "interactive_plots" : False
}

run_neutral_mixed_manufactured_solutions_test(test_input)
