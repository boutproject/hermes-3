# Changelog

## Unreleased

- Make `evolve_vorticity_ion_energy` polarisation work consistent with the
  pressure potential and perpendicular kinetic energy. Use uniform radial-wall
  particle, momentum and energy weights formed from separate wall-line means.
- Carry ion viscous drift work through polarisation energy transport;
  retain `DivEviscous` as a diagnostic and preserve electron heating transfer.
- Default to ten auxiliary polarisation corrections and identify the transport
  contract with component metadata version 4.
- Add generated-restart gauge, wall-energy and slow fixed-band wall-amplitude
  regressions; update native-output flux and spectral analysis.
