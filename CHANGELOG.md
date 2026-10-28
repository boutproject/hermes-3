# Changelog

## Unreleased

### Changed

- The `diamagnetic_drift` gradient form always uses second-order central
  (`C2`) derivatives, independent of the mesh's default derivative methods.
  With `radial_boundary = closed` or `closed_average`, it removes the
  implied fluxes through the inner and outer radial faces before any core
  averaging. The `open` setting retains those fluxes.
