# 1D low-density sources integration test

This is an extremely simple 1D example designed to test the low_n_sources component (the LowNSources and PseudoReaction classes).
The only top-level components are the species "t", "t+", "d", "d+" and "e"; and low_n_sources.
Density, pressure and momentum are evolved for all species with no-flow boundaries.
All species densities are started at values that trigger the low-n sources, restoring each to its respective threshold (10x the species-specific density floor, by default.)
