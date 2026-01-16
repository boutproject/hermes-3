#!/usr/bin/env python3
from job_functions import run_neutral_mixed_manufactured_solutions_test
from boutdata.mms import x, y, z, Div_par, Grad_par
from perpendicular_laplacian import Div_a_Grad_perp_f, Div_par_k_Grad_par_f, GeneralMetric
from sympy import sin, cos, log

# specify symbolic inputs
# contravariant metric coeffs
# g11 = 1.1 + 0.16*x*cos(y)
# g22 = 0.9 + 0.09*x*cos(y)
# g33 = 1.2 + 0.2*x*cos(y)
# g12 = 0.0
# g23 = 0.5 + 0.15*x*cos(y)
# g13 = 0.0
# we can safely specify the simplest metric
# here knowing that the symbolic operators
# are used in orthogonal_test.py and nonorthogonal_test.py
# to test more complex metrics.
g11 = 1.0
g22 = 1.0
g33 = 1.0
g12 = 0.0
g23 = 0.0
g13 = 0.0
# the metric object
metric = GeneralMetric(
            g11=g11, g12=g12, g13=g13,
            g22=g22, g23=g23, g33=g33)
# extract covariant metric coeffs and Jacobian
g_11 = metric.g_11
g_12 = metric.g_12
g_13 = metric.g_13
g_22 = metric.g_22
g_23 = metric.g_23
g_33 = metric.g_33
J = metric.J

def neutral_mixed_equations(evolve_momentum=False,
                            conservation_test=False,
                            # expected N in discretisation error ~ (grid_spacing)^N
                            expected_convergence_order=2.0):
    pfac = 0.2
    nfac = 0.2
    # mass
    AA = 2.0
    # collision frequency
    nu_cfreq = 1.0
    # Manufactured solutions
    Nn = 0.5 + nfac*(0.5 + x**3 - 1.5*x**2)*sin(y)
    Pn = 1.0 + pfac*(0.5 + x**3 - 1.5*x**2)*sin(2*y)
    if evolve_momentum:
        # choose function that is both zero on boundary in x
        # and has no derivative
        NVn = (x**2 - 2*x**3 + x**4)*sin(y)
        #NVn = 0.0 + x*(x-1)*cos(y)
    else:
        NVn = 0.0

    # Derived quantities
    Vn = NVn/(Nn*AA)
    Tn = Pn/Nn
    logPn = log(Pn)
    Dn = (Tn/AA)/nu_cfreq
    Kn = (5.0/2.0)*Nn*Dn
    Etan = (2.0/5.0)*AA*Kn
    neutral_conduction = True
    neutral_viscosity = True

    # time evolution equations without sources
    # N.B.
    # Div . ( a Grad_perp f )
    #div_a_grad_perp_f = div_a_grad_perp_f_symbolic(g11, g12, g13, g22, g23, g33, a, f)
    # Div . ( a Grad_par f)
    #div_par_k_grad_par_f = div_par_k_grad_par_f_symbolic(g11, g12, g13, g22, g23, g33, a, f)
    # Div . ( vec(b) f)
    #div_par_f = div_par_f_symbolic(g11, g12, g13, g22, g23, g33, f)
    # vec(b) . Grad f
    #grad_par_f = grad_par_f_symbolic(g11, g12, g13, g22, g23, g33, f)

    ddt_Nn = (  -Div_par(Nn*Vn, metric=metric)
                -Div_a_Grad_perp_f(-Dn*Nn, logPn, metric=metric)
                )

    ddt_Pn = (  -Div_par(Pn*Vn, metric=metric)
                -(5.0/3.0)*Div_a_Grad_perp_f(-Dn*Pn, logPn, metric=metric)
                -(2.0/3.0)*Pn*Div_par(Vn, metric=metric)
                )

    ddt_NVn = ( - AA*Div_par(Nn*Vn*Vn, metric=metric)
                + AA*Div_a_Grad_perp_f(Nn*Vn*Dn, logPn, metric=metric)
                - Grad_par(Pn, metric=metric)
                )
    # if neutral conduction included in test
    if neutral_conduction:
        ddt_Pn = (ddt_Pn +
                +(2.0/3.0)*Div_a_Grad_perp_f(Kn, Tn, metric=metric)
                +(2.0/3.0)*Div_par_k_Grad_par_f(Kn, Tn, metric=metric)
                )
    # if neutral viscosity included in test
    if neutral_viscosity:
        viscosity_source = (Div_a_Grad_perp_f(Etan, Vn, metric=metric)
                        + Div_par_k_Grad_par_f(Etan, Vn, metric=metric))
        ddt_NVn = (ddt_NVn +
                    viscosity_source
                    )
        ddt_Pn = (ddt_Pn -
                (2.0/3.0)*Vn*viscosity_source
                )

    # Sources for steady state, set to exactly cancel
    # the RHS operators if aiming to test the definition
    # of the RHS, or set to zero if we just want to test
    # that the RHS is conservative.
    if conservation_test:
        # we just want to make ddt and test conservation properties
        # so no sources are required here
        source_Nn = 0.0
        source_Pn = 0.0
        source_NVn = 0.0
    else:
        # we want to cancel the numerically calculated ddt with an
        # analytical source, so set the sources
        source_Nn = -ddt_Nn
        source_Pn = -ddt_Pn
        source_NVn = -ddt_NVn

    test_input = {
        "ntest" : 3,
        "ngrid" : 10,
        "collision_frequency" : nu_cfreq,
        "g11_string": str(g11),
        "g22_string": str(g22),
        "g33_string": str(g33),
        "g12_string": str(g12),
        "g13_string": str(g13),
        "g23_string": str(g23),
        "g_11_string": str(g_11),
        "g_22_string": str(g_22),
        "g_33_string": str(g_33),
        "g_12_string": str(g_12),
        "g_13_string": str(g_13),
        "g_23_string": str(g_23),
        "J_string": str(J),
        "mass" : AA,
        "Nd_string" : str(Nn),
        "Pd_string" : str(Pn),
        "NVd_string" : str(NVn),
        "source_Nd_string" : str(source_Nn),
        "source_Pd_string" : str(source_Pn),
        "source_NVd_string" : str(source_NVn),
        "neutral_conduction" : neutral_conduction,
        "neutral_viscosity" : neutral_viscosity,
        "evolve_momentum" : evolve_momentum,
        "test_dir" : "neutral_mixed",
        "sub_test_dir" : "evolve_momentum_"+str(evolve_momentum)+"_conservation_test_"+str(conservation_test),
        "interactive_plots" : False,
        "conservation_test" : conservation_test,
        "expected_convergence_order" : expected_convergence_order,
    }
    return test_input

global_success = True

test_input = neutral_mixed_equations(evolve_momentum=False,
                                    conservation_test=False,
                                    expected_convergence_order=2.0)
success, output_message_0 = run_neutral_mixed_manufactured_solutions_test(test_input)
global_success = success and global_success

test_input = neutral_mixed_equations(evolve_momentum=False,
                                    conservation_test=True,
                                    expected_convergence_order=2.0)
success, output_message_1 = run_neutral_mixed_manufactured_solutions_test(test_input)
global_success = success and global_success

test_input = neutral_mixed_equations(evolve_momentum=True,
                                    conservation_test=False,
                                    expected_convergence_order=1.0)
success, output_message_2 = run_neutral_mixed_manufactured_solutions_test(test_input)
global_success = success and global_success

test_input = neutral_mixed_equations(evolve_momentum=True,
                                    conservation_test=True,
                                    expected_convergence_order=0.8)
success, output_message_3 = run_neutral_mixed_manufactured_solutions_test(test_input)
global_success = success and global_success

print(output_message_0)
print(output_message_1)
print(output_message_2)
print(output_message_3)

if global_success:
    exit(0)
else:
    exit(1)
