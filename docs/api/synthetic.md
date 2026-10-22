# `virgil_vlti.synthetic`

Simulated GRAVITY products with a known truth: `write_product` writes a noise-free product whose V² and closure phases are functions you give, and `add_systematics` adds noise, gains and closure-phase offsets of known size. The [tutorials](../calibration.md) and the test suite use them.

::: virgil_vlti.synthetic
    options:
      members:
        - write_product
        - add_systematics
