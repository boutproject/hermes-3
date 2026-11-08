# 2D interchange energy regression

The legacy pressure formulation checks the signs of electron, ion and total
pressure changes in a curvature-driven slab with viscosity.

```sh
ctest --test-dir build -R '^2D-energy$' --output-on-failure
```

See the [Sphinx description](https://hermes3.readthedocs.io/en/latest/developer.html#subsec-2d-energy-test)
for the energy and entropy exchanges and the limits of these checks.
