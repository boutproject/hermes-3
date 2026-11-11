# 1D Sheath Boundaries

This is a low-resolution version of the 1D-sheath example. This tests
the `sheath_boundary_simple` component under a few different
cases. This is purely a golden answer regression test.

The golden answers can be updated with:

```console
$ ./test_sheath.py --snaptol-update-all
```

This test requires `pytest` and `snaptol`: install these with `uv sync
--only-dev`.
