# 1D heat conduction

Solves electron heat conduction in 1D on a non-uniform grid as
used in the 1D-recycling and other examples.

The main purpose of this is to test preconditioning methods.

With CVODE solver, no preconditioning:

    ./hermes-3 -d examples/1D-conduction/ solver:use_precon=false

    2.000e+03       1523       5.23e-01    89.3    0.0    0.3    4.7    5.6
    4.000e+03       1225       4.27e-01    89.9    0.0    0.3    3.9    5.9
    6.000e+03       1057       3.64e-01    89.0    0.0    0.3    4.6    6.1
    8.000e+03       3563       1.18e+00    93.6    0.0    0.3    1.4    4.7

CVODE solver, with preconditioning:

     ./hermes-3 -d examples/1D-conduction/ solver:use_precon=true

    2.000e+03        209       1.18e-01    57.6    0.0    0.2   20.1   22.1
    4.000e+03         85       5.84e-02    47.5    0.0    0.2   28.0   24.3
    6.000e+03         81       5.81e-02    46.0    0.0    0.2   28.2   25.7
    8.000e+03        130       7.73e-02    54.5    0.0    0.2   20.9   24.4
