# Included benchmarks

These truth tables are the complete benchmark subset used by the documented
experiments, copied from [`changmg/DALTA`](https://github.com/changmg/DALTA)
commit `f8458f7` so that no separate DALTA checkout is required.

- `numerics/`: the six 9-input functions used in Table 1 (`cos`, `tan`,
  `exp`, `ln`, `erf`, and `denoise`);
- `axbench-verilog/adder_8.txt`: the 16-input, 9-output Brent-Kung adder used
  for the Figure 4 run.

Each line contains an integer input address followed by its exact integer
output value. The solver infers the input and output widths from the table.

The upstream DALTA snapshot does not contain a project-wide license file.
Confirm redistribution permission with the DALTA authors before publishing
these copied datasets in a public repository.
