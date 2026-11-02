# Reproduction status

## Environment

- macOS / Apple M3 Pro
- AppleClang 15
- Eigen 3.4.0
- Gurobi 12.0.2 (the papers used Gurobi 9.1.1)

## DALTA baseline

The minimal adapted DALTA snapshot in `third_party/dalta-baseline/` builds on
macOS and exposes all previously hard-coded solver paths:

- `joint-heuristic`: DALTA's alternating row-pattern/row-type heuristic;
- `joint-ilp`: DALTA's joint Gurobi model;
- `separate-legacy`: the released separate Gurobi formulation;
- `separate-ilp`: a portable, explicitly binary reformulation of the same
  separate objective.

The full six-benchmark joint heuristic run is stored in
`results/table1/results.json`. With seed 1 it obtains average MED 3.790 versus
3.61 in Table 1. Individual values are close except for the stochastic `ln`
and `denoise` runs. Average runtime is not expected to match across CPUs.

The joint ILP path has been executed end-to-end with one partition and one
round (cos: MED 2.57227, 3.72 s). The paper's full configuration explores 126
partitions for five rounds and reports the 3600-second timeout for every
benchmark. It is therefore opt-in in the experiment driver.

The released separate formulation is not practically compatible with Gurobi
12: even one fixed-partition model can take more than 30 seconds, whereas the
paper reports roughly 140--260 seconds for the complete benchmark under
Gurobi 9.1. The portable formulation is retained to separate this solver
version issue from the mathematical model.

A complete 126-partition `separate-ilp` cos run with a 0.2-second limit per
model and 10% MIP gap obtains MED 10.3984 in 333.1 seconds. Table 1 reports
11.64 in 258.37 seconds. This confirms the full path and approximately matches
the reported quality despite the different Gurobi version. Its command and
raw output are stored in `results/table1/dalta-separate-cos.json`.

## Ising solver

`src/` is a clean-room C++/Eigen implementation containing:

1. column-based approximate disjoint decomposition;
2. separate and joint second-order Ising coefficient construction;
3. ballistic simulated bifurcation with symplectic Euler updates and
   inelastic walls;
4. the paper's energy-variance stopping rule (`t=s=20` for 9 inputs and
   `t=s=10` for 16 inputs, epsilon `1e-8`);
5. exact conditional updates of the column-type vector;
6. MSB-to-LSB optimization repeated for five rounds;
7. exhaustive partitions for 9 inputs and bounded random exploration for
   larger functions.

### Paper preset

The `paper` preset uses the paper's free-row orientation,
64-spin Ising core, separate ER and joint MED formulations, bSB dynamics,
20/20 energy-variance stop, exact conditional column-type feedback, all 126
partitions, and five MSB-to-LSB rounds. It uses exactly one complete seeded bSB
trajectory per core and disables the enhanced portfolio, MED-based agent
selection, separate-mode ER slack, and all hybrid refinements.

The target paper does not state how the linear Ising field is scheduled inside
bSB. This preset makes one documented reconstruction choice per formulation:
`eta(t)=0.7*sqrt(p(t))` for separate ER and `eta(t)=0.329*c0*p(t)` for joint
MED. Each mode still executes exactly one schedule, rather than selecting
among a portfolio. The integration step, pressure slope, initialization,
coupling normalization, 2000-step safety cap, and seed are also necessarily
inferred because the paper does not report them.

| Benchmark | `paper` separate | Reported separate | `paper` joint | Reported joint |
|---|---:|---:|---:|---:|
| cos | 8.174 | 8.33 | 2.404 | 2.50 |
| tan | 7.105 | 10.45 | 2.494 | 2.50 |
| exp | 6.105 | 7.07 | 2.688 | 2.66 |
| ln | 7.180 | 6.57 | 2.508 | 2.72 |
| erf | 7.381 | 4.61 | 2.447 | 1.90 |
| denoise | 8.576 | 9.69 | 2.543 | 2.80 |
| average | **7.420** | **7.83** | **2.514** | **2.51** |

Average measured runtime in the latest run is 0.26 seconds for separate and
1.32 seconds for joint. The complete commands and raw outputs are stored in
`results/table1/paper-results.json`; the field calibration is recorded in
`results/table1/FIELD_SWEEP.md`.

The remaining per-function gaps are caused by schedule sensitivity and, in
separate mode, objective mismatch. Separate mode minimizes unweighted bit ER
but Table 1 reports numeric MED; additional ER-minimizing restarts actually
increased the six-function average MED from 7.420 to 8.098, 8.367, and 8.798
for two, four, and eight restarts. Seeds 1--10, smaller initialization,
coupling scales, pressure ramps, integration caps, tie policies, and five
separate rounds were also checked without removing the `erf` gap. A different
single field schedule can obtain `erf` MED 4.160, but degrades the other five
functions, so there is no single tested calibration that covers both basins.

### Paper-plus preset

`paper-plus` preserves the `paper` preset and adds only a fixed
single-seed field portfolio: two complete trajectories for separate and three
for joint, run in parallel. The complete LUT with the lowest MED is retained.
It does not use multiple seeds, ER slack, local refinement, simulated
annealing, or the Gurobi oracle. Final MED selection is not reported by the
paper, so this preset is explicitly an enhancement rather than strict
reproduction.

| Benchmark | `paper-plus` separate | Reported separate | `paper-plus` joint | Reported joint |
|---|---:|---:|---:|---:|
| cos | 8.174 | 8.33 | 2.404 | 2.50 |
| tan | 7.105 | 10.45 | 2.494 | 2.50 |
| exp | 6.105 | 7.07 | 2.584 | 2.66 |
| ln | 7.180 | 6.57 | 2.508 | 2.72 |
| erf | 4.160 | 4.61 | 1.783 | 1.90 |
| denoise | 8.576 | 9.69 | 2.543 | 2.80 |
| average | **6.883** | **7.83** | **2.386** | **2.51** |

Five of six separate results and all six joint results are no worse than the
paper. Average measured runtime in the recorded run is 0.23 seconds for
separate and 1.91 seconds for joint. Full records are stored in
`results/table1/paper-plus-results.json`.

## Figure 4: Brent-Kung

The full 16-input Brent-Kung example has also been run with `paper` in joint
mode using `n-free=7`, 1000 partitions, five rounds, and seed 1. It obtains MED
0.839844 in 784.51 seconds. Figure 4 prints DALTA MED 1.11 and runtime 1886.98
seconds for this benchmark, giving reproduced ratios 0.7566 and 0.4157,
respectively. The exact command and machine-readable record are documented in
`results/figure4/`.

## Commands

Build the Ising solver and the optional DALTA baseline in one build tree:

```sh
cmake -S . -B build \
  -DBUILD_DALTA_BASELINE=ON \
  -DGUROBI_HOME=/Library/gurobi1202/macos_universal2 \
  -DEigen3_DIR=/opt/homebrew/share/eigen3/cmake
cmake --build build -j
```

Run the normal self-contained Ising Table-1 suite:

```sh
python3 experiments/run_table1.py
```

Run the single-trajectory `paper` suite:

```sh
python3 experiments/run_table1.py \
  --methods ising_separate_paper ising_joint_paper \
  --output results/table1/paper-results.json
```

Run the minimal field-portfolio enhancement:

```sh
python3 experiments/run_table1.py \
  --methods ising_separate_paper_plus ising_joint_paper_plus \
  --output results/table1/paper-plus-results.json
```

Opt into a long Gurobi run:

```sh
python3 experiments/run_table1.py \
  --methods dalta_joint_ilp --timeout 3700
```

Run the bounded separate-ILP check used above:

```sh
python3 experiments/run_table1.py \
  --methods dalta_separate_ilp --benchmarks cos \
  --mip-time 0.2 --mip-gap 0.1 --timeout 600 \
  --output results/table1/dalta-separate-cos-rerun.json
```
