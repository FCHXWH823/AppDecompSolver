# Figure 4 reproduction

The first reproduced 16-input example is Brent-Kung. The included DALTA
benchmark file `benchmarks/axbench-verilog/adder_8.txt` enumerates the addition
of two 8-bit inputs, giving 16 input bits and 9 output bits, as described in
the paper.

The full paper-scale `paper` run uses joint mode, a 7-bit free set, a 9-bit
bound set, 1000 randomly selected partitions, five rounds, and seed 1. For a
16-input truth table, the solver automatically uses the paper's dynamic-stop
settings `t=s=10` and `epsilon=1e-8`.

```sh
./build/ising-decomp \
  --input benchmarks/axbench-verilog/adder_8.txt \
  --mode joint --n-free 7 --rounds 5 --partitions 1000 \
  --seed 1 --max-steps 2000 --preset paper
```

Measured result:

| Method | MED | Time (s) |
|---|---:|---:|
| This reproduction (`paper`) | 0.839844 | 784.51 |
| DALTA value printed in Figure 4 | 1.11 | 1886.98 |

Relative to the DALTA numbers printed in Figure 4, the reproduced result has
MED ratio 0.7566, runtime ratio 0.4157, and a 2.41x speedup. Figure 4 presents
the proposed method as ratios rather than printing its absolute MED/runtime,
so this comparison does not claim an exact absolute match to an unprinted
paper value.

The machine used here is an Apple M3 Pro; the paper used a 1.9 GHz AMD Ryzen 7
5800U. Runtime comparisons are therefore indicative rather than controlled.
The machine-readable record is in `brent-kung-paper.json`.
