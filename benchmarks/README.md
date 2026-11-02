# Included benchmarks

These truth tables are the complete benchmark subset used by the documented
experiments, copied from [`changmg/DALTA`](https://github.com/changmg/DALTA)
commit `f8458f7` so that no separate DALTA checkout is required.

- `numerics/`: the six functions used in Table 1 at 9 inputs and Figure 4 at
  16 inputs (`cos`, `tan`, `exp`, `ln`, `erf`, and `denoise`);
- `axbench-verilog/adder_8.txt`: the 16-input, 9-output Brent-Kung adder;
- `axbench-verilog/kin_16.txt`, `kin_inv.txt`, and `mult_8.txt`: the
  Forwardk2j, Inversek2j, and Multiplier truth tables from Figure 4.

Each line contains an integer input address followed by its exact integer
output value. The solver infers the input and output widths from the table.

The referenced upstream DALTA snapshot does not contain a project-wide license
file. The datasets retain their upstream attribution and are included here for
experiment reproducibility; consult the DALTA authors about reuse or
redistribution beyond this artifact.
