# Figure 4 reproduction

All ten 16-input benchmarks from Figure 4 were run in joint mode with the
`paper` preset, `n_free=7`, `n_bound=9`, `P=1000`, `R=5`, `seed=1`, and one
bSB trajectory per core problem. For 16-input truth tables, the dynamic-stop
settings are `t=s=10` and `epsilon=1e-8`.

For `n>10`, the solver now exactly mirrors DALTA's `DrawVarOrd` partition
selection: it starts from the natural free set `{0,...,6}`, repeatedly
shuffles all input variables, takes and sorts the first seven variables, and
allows duplicate partitions. A direct comparison confirmed that DALTA's
`boost::random::mt19937` and this solver's `std::mt19937` generate the same
1,000 free sets, in the same order, for seed 1.

## Results

| Benchmark | This reproduction MED | Figure 4 DALTA MED | MED ratio | Time (s) | Figure 4 DALTA time (s) |
|---|---:|---:|---:|---:|---:|
| cos(x) | 9.701981 | 10.18 | 0.9530 | 1563.09 | 4477.86 |
| tan(x) | 2.967453 | 3.14 | 0.9450 | 1481.60 | 3389.72 |
| exp(x) | 10.687546 | 10.30 | 1.0376 | 1601.86 | 3865.39 |
| ln(x) | 10.075836 | 10.65 | 0.9461 | 1487.36 | 5312.16 |
| erf(x) | 10.063354 | 15.01 | 0.6704 | 1808.24 | 3983.88 |
| denoise(x) | 9.858566 | 33.94 | 0.2905 | 1594.30 | 4503.74 |
| Brent-Kung | 0.970703 | 1.11 | 0.8745 | 846.06 | 1886.98 |
| Forwardk2j | 904.525009 | 1168.18 | 0.7743 | 1738.42 | 8261.07 |
| Inversek2j | 359.461868 | 379.68 | 0.9467 | 2007.20 | 6970.90 |
| Multiplier | 426.158341 | 481.47 | 0.8851 | 2028.32 | 8506.03 |

The mean of the ten per-case MED ratios is **0.8323**, so this reproduction
has a 16.8% smaller MED than the Figure 4 DALTA values on average. It improves
MED on nine of ten cases. The paper reports an average ratio of approximately
0.89 (11% smaller MED) and improvements on six of ten cases.

The mean per-case runtime ratio is 0.3474, corresponding to 2.88x relative to
the DALTA times printed in the paper. This runtime comparison is indicative,
not controlled: these runs used an Apple M3 Pro and three concurrent jobs,
whereas the paper used a 1.9 GHz AMD Ryzen 7 5800U.

The complete commands, stdout, timings, and ratios are in `results.json`.
Run or resume all cases with:

```sh
python3 experiments/run_figure4.py --jobs 3
```

## Why the earlier MEDs were anomalously large

The earlier implementation enumerated all `C(16,7)=11,440` partitions and
sampled 1,000 unique combinations uniformly. That is not DALTA's strategy and
changes both the candidate partitions and the bSB random-number consumption.
For ordered numeric LUTs, DALTA's first natural partition preserves the
low-bit/high-bit structure of the fixed-point input and is especially strong.

For example, on `tan(x)`:

| Configuration | MED |
|---|---:|
| Old unique-combination sampling, `P=1000, R=5` | 13.436905 |
| Correct DALTA draw, natural partition only, `P=1, R=5` | 3.359024 |
| Correct DALTA draw, `P=1000, R=1` | 3.094467 |
| Correct DALTA draw, full `P=1000, R=5` | 2.967453 |

Joint decomposition processes output bits from MSB to LSB. A poor partition
can therefore introduce a high-weight error early, which changes the `D` and
`q` coefficients seen by every later core problem and amplifies the final
MED. The old partial results are retained for diagnosis in
`results-legacy-combination-sampling.json`; they are not reproduction data.
