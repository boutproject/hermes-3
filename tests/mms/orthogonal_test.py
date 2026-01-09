#!/usr/bin/env python3
from job_functions import run_manufactured_solutions_test
from boutdata.mms import x, y, z, Div_par, Grad_par
from perpendicular_laplacian import Div_a_Grad_perp_f, Div_par_k_Grad_par_f, GeneralMetric
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
# the metric object
metric = GeneralMetric(
            g11=g11, g12=g12, g13=g13,
            g22=g22, g23=g23, g33=g33)
# Div . ( a Grad_perp f )
div_a_grad_perp_f = Div_a_Grad_perp_f(a, f, metric=metric)
# Div . ( a Grad_par f)
div_par_k_grad_par_f = Div_par_k_Grad_par_f(a, f, metric=metric)
# Div . ( vec(b) f)
div_par_f = Div_par(f, metric=metric)
# vec(b) . Grad f
grad_par_f = Grad_par(f, metric=metric)

test_input = {
    "ntest" : 3,
    "ngrid" : 20,
    # list of list of ["name", symbolic function, expected convergence order]
    "differential_operator_list": [["FV::Div_a_Grad_perp(a, f)", str(div_a_grad_perp_f), 2],
                                   ["Div_a_Grad_perp_nonorthog(a, f)", str(div_a_grad_perp_f), 2],
                                   ["Div_a_Grad_perp_flows(a, f)", str(div_a_grad_perp_f), 2],
                                   ["Div_par_K_Grad_par_mod(a, f)", str(div_par_k_grad_par_f), 2],
                                   ["Div_par(f)", str(div_par_f), 2],
                                   ["FV::Div_par_fvv(f, v, wave_speed)", str(div_par_f), 1],
                                   ["FV::Div_par_mod(f, v, wave_speed)", str(div_par_f), 1],
                                   ["Grad_par(f)", str(grad_par_f), 2],
                                   ],
    "a_string": str(a),
    "f_string": str(f),
    "g11_string": str(g11),
    "g22_string": str(g22),
    "g33_string": str(g33),
    "g12_string": str(g12),
    "g13_string": str(g13),
    "g23_string": str(g23),
    "test_dir" : "orthogonal",
    "interactive_plots" : False
}

success, output_message = run_manufactured_solutions_test(test_input)

print(output_message)
if success:
    exit(0)
else:
    exit(1)
