# Changelog

## Unreleased

- Use Neumann radial ion-pressure and pressure-potential guards with a conducting,
  Z-uniform electrostatic wall potential. The generalised potential includes the
  local wall pressure potential.
- Carry polarisation wall work, particles, momentum and thermal/parallel kinetic
  energy with local face weights. Identify this contract with ion-energy version
  1 and conservative drift version 1; retain the auxiliary chi approximation.
- Remove analysis support for unreleased ion-energy versions and untagged data;
  require the explicit ion-energy version-1 tag.
- Carry ion viscous drift work through polarisation energy transport;
  retain `DivEviscous` as a diagnostic and preserve electron heating transfer.
- Default to ten auxiliary polarisation corrections.
- Add generated-restart gauge, wall-energy and slow fixed-band wall-amplitude
  regressions, exact frozen-state wall convergence and varying-fraction
  multispecies identities. Add a slow D = 0.1 publication convergence regression;
  update native-output analysis with local wall transport reconstruction.
