# Efficient Approximate Decomposition Solver using Ising Model (DAC 2024)

> **Abstract.** Computing with memory is an energy-efficient computing
> approach. It pre-computes a function and stores its values in a lookup table
> (LUT), which can be retrieved at runtime. Approximate Boolean decomposition
> reduces the LUT size for implementing complex functions, but it takes a long
> time to find a decomposition with a minimal error. In this work, to address
> this issue, we propose an efficient Ising model-based approximate Boolean
> decomposition solver. First, a new column-based approximate disjoint
> decomposition method is proposed to fit the Ising model. Then, it is adapted
> to the Ising model-based optimization solver. Moreover, two improvement
> techniques are developed for an efficient search of the approximate disjoint
> decomposition when using simulated bifurcation to solve the Ising model.
> Experimental results show that compared to the state-of-the-art work, our
> approach achieves an 11% smaller mean error distance with an average 1.16x
> speedup when approximately decomposing 16-input Boolean functions.

This repository is a reproducible code artifact for the DAC 2024 paper by
Weihua Xiao, Tingting Zhang, Xingyue Qian, Jie Han, and Weikang Qian,
*Efficient Approximate Decomposition Solver using Ising Model*. The primary
solver is our clean-room C++/Eigen reproduction, not an official source release
from the paper's authors.

The repository is self-contained for the experiments documented here. The
required DALTA truth tables are under `benchmarks/`; users do **not** need to
clone or download the DALTA repository. A minimal adapted DALTA comparison
snapshot is included under `third_party/dalta-baseline/` and is only needed
when rebuilding the optional baseline.

## Implemented paper components

- column-based decomposition with two column patterns and one column-type vector;
- the second-order Ising formulations for separate (ER) and joint (MED) modes;
- ballistic simulated bifurcation with symplectic Euler integration;
- variance-based dynamic stopping;
- exact conditional updates of the column-type vector during energy sampling;
- DALTA's outer MSB-to-LSB, multi-round optimization and partition exploration.

## Repository layout

- `src/`: Ising decomposition solver;
- `benchmarks/`: the six Table 1 truth tables and the Figure 4 Brent-Kung case;
- `experiments/`: reproducible experiment driver;
- `results/`: recorded Table 1 and Figure 4 outputs;
- `third_party/dalta-baseline/`: optional adapted DALTA baseline snapshot;
- `TECHNIQUES.md` and `REPRODUCTION.md`: technique and reproduction details.

## Installation

The main bSB solver is CPU-only and requires:

- CMake 3.16 or newer;
- a C++20 compiler (GCC 10+, Clang 12+, or AppleClang 14+);
- Eigen 3.3 or newer;
- Python 3.9+ only when running the experiment driver.

Gurobi and Boost are **not** required for the normal Ising solver.

### Ubuntu/Debian

Install the build dependencies:

```sh
sudo apt update
sudo apt install -y build-essential cmake libeigen3-dev python3 git
```

Clone and build:

```sh
git clone https://github.com/FCHXWH823/AppDecompSolver.git
cd AppDecompSolver
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

### macOS

Install the Apple command-line tools and Homebrew dependencies:

```sh
xcode-select --install
brew install cmake eigen python
```

Then clone and build:

```sh
git clone https://github.com/FCHXWH823/AppDecompSolver.git
cd AppDecompSolver
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DEigen3_DIR="$(brew --prefix eigen)/share/eigen3/cmake"
cmake --build build --parallel
```

### Verify the installation

Run a fast one-partition smoke test from the repository root:

```sh
./build/ising-decomp \
  --input benchmarks/numerics/cos-i-9-o-10-a-0-b-1_57079632-acc.txt \
  --mode joint --n-free 4 --rounds 1 --partitions 1 \
  --seed 1 --max-steps 2000 --preset paper --quiet
```

A successful run ends with a line beginning with `RESULT mode=joint`.

### Optional DALTA/Gurobi baseline

This step is needed only to rerun the DALTA heuristic or ILP baselines. Install
Boost and a licensed Gurobi distribution first. On Ubuntu, Boost can be
installed with `sudo apt install libboost-all-dev`; on macOS, use
`brew install boost`. Set `GUROBI_HOME` to the platform directory that contains
Gurobi's `include/` and `lib/` directories, then reconfigure the same build:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_DALTA_BASELINE=ON \
  -DGUROBI_HOME=/path/to/gurobi/platform-directory
cmake --build build --parallel
./build/dalta-baseline -h
```

The vendored baseline source and all required truth tables are already in this
repository; a separate DALTA clone is not needed.

### Troubleshooting

- If CMake cannot find Eigen, pass the directory containing
  `Eigen3Config.cmake` through `-DEigen3_DIR=/path/to/eigen3/cmake`.
- Run all example and experiment commands from the repository root because the
  documented benchmark paths are relative to it.
- If Gurobi is not detected, verify that both
  `$GUROBI_HOME/include/gurobi_c++.h` and the Gurobi libraries under
  `$GUROBI_HOME/lib` exist. The bSB solver can still be built by leaving
  `BUILD_DALTA_BASELINE=OFF`.
- The recorded runtimes are machine-dependent; MED should be used for the
  primary quality comparison.

## Quick start

Example matching the 9-input configuration in Table 1:

```sh
./build/ising-decomp \
  --input benchmarks/numerics/cos-i-9-o-10-a-0-b-1_57079632-acc.txt \
  --mode joint --n-free 4 --rounds 5 --partitions 1000 --seed 1
```

Two reproducibility presets are provided:

- `--preset paper` (default) follows the paper technique stack closely: free
  variables index rows, one ballistic trajectory is used per core, dynamic
  stopping is 20/20, and the exact conditional type intervention is enabled.
  It uses one fixed independent external-field calibration per formulation to
  fill a parameter omitted by the paper: `0.7*sqrt(p)` for separate and
  `0.329*c0*p` for joint. It has no portfolio, ER slack, or hybrid refinement.
- `--preset paper-plus` is the smallest quality-oriented extension of
  the `paper` preset. It runs two fixed field trajectories for separate and
  three for joint, all with seed 1, and selects the final LUT by MED. It still
  has no ER slack, local search, annealing, or exact-solver oracle.

Example:

```sh
./build/ising-decomp \
  --input benchmarks/numerics/cos-i-9-o-10-a-0-b-1_57079632-acc.txt \
  --mode joint --n-free 4 --rounds 5 --partitions 1000 \
  --preset paper-plus
```

Run the full six-function `paper` comparison:

```sh
python3 experiments/run_table1.py \
  --methods ising_separate_paper ising_joint_paper \
  --output results/table1/paper-results.json
```

Run the minimally enhanced `paper-plus` comparison:

```sh
python3 experiments/run_table1.py \
  --methods ising_separate_paper_plus ising_joint_paper_plus \
  --output results/table1/paper-plus-results.json
```

The active-technique comparison between the target paper and both presets is
listed in `TECHNIQUES.md`.

## Best tested single-trajectory configuration per case

The configurations below are manual, per-benchmark calibrations of the
`paper` solver. Every command runs exactly one complete bSB trajectory with
seed 1; no portfolio or post-run selection is used. They are the best retained
configurations found in this reproduction study, not parameters published by
the paper and not one universal configuration. Separate uses one pass because
its per-bit ER problems are independent; joint uses the paper's five rounds.

| Mode | Benchmark | Field schedule | Scale | Normalization | MED | Paper MED |
|---|---|---|---:|---|---:|---:|
| separate | cos | sqrt | 0.7 | full | 8.174 | 8.33 |
| separate | tan | sqrt | 0.7 | full | 7.105 | 10.45 |
| separate | exp | sqrt | 0.7 | full | 6.105 | 7.07 |
| separate | ln | coupled-linear | 0.5 | pairwise | 6.230 | 6.57 |
| separate | erf | coupled-linear | 1.0 | full | 4.160 | 4.61 |
| separate | denoise | sqrt | 0.7 | full | 8.576 | 9.69 |
| joint | cos | coupled-linear | 0.329 | full | 2.404 | 2.50 |
| joint | tan | coupled-linear | 0.329 | full | 2.494 | 2.50 |
| joint | exp | linear | 1.15 | full | 2.584 | 2.66 |
| joint | ln | coupled-linear | 0.329 | full | 2.508 | 2.72 |
| joint | erf | coupled | 1.0 | pairwise | 1.783 | 1.90 |
| joint | denoise | coupled-linear | 0.329 | full | 2.543 | 2.80 |

### General-case commands

For separate `cos`, `tan`, `exp`, and `denoise`, use the following general
configuration and replace only the input path:

```sh
./build/ising-decomp \
  --input benchmarks/numerics/cos-i-9-o-10-a-0-b-1_57079632-acc.txt \
  --mode separate --n-free 4 --rounds 1 --partitions 1000 --seed 1 \
  --max-steps 2000 --preset paper --field-schedule sqrt \
  --field-scale 0.7 --normalization full --quiet
```

For joint `cos`, `tan`, `ln`, and `denoise`, use this general configuration
and replace only the input path:

```sh
./build/ising-decomp \
  --input benchmarks/numerics/cos-i-9-o-10-a-0-b-1_57079632-acc.txt \
  --mode joint --n-free 4 --rounds 5 --partitions 1000 --seed 1 \
  --max-steps 2000 --preset paper --field-schedule coupled-linear \
  --field-scale 0.329 --normalization full --quiet
```

### Four special-case commands

Separate `ln`:

```sh
./build/ising-decomp \
  --input benchmarks/numerics/ln-i-9-o-10-a-1-b-10-acc.txt \
  --mode separate --n-free 4 --rounds 1 --partitions 1000 --seed 1 \
  --max-steps 2000 --preset paper --field-schedule coupled-linear \
  --field-scale 0.5 --normalization pairwise --quiet
```

Separate `erf`:

```sh
./build/ising-decomp \
  --input benchmarks/numerics/erf-i-9-o-10-a-0-b-3-acc.txt \
  --mode separate --n-free 4 --rounds 1 --partitions 1000 --seed 1 \
  --max-steps 2000 --preset paper --field-schedule coupled-linear \
  --field-scale 1.0 --normalization full --quiet
```

Joint `exp`:

```sh
./build/ising-decomp \
  --input benchmarks/numerics/exp-i-9-o-10-a-0-b-3-acc.txt \
  --mode joint --n-free 4 --rounds 5 --partitions 1000 --seed 1 \
  --max-steps 2000 --preset paper --field-schedule linear \
  --field-scale 1.15 --normalization full --quiet
```

Joint `erf`:

```sh
./build/ising-decomp \
  --input benchmarks/numerics/erf-i-9-o-10-a-0-b-3-acc.txt \
  --mode joint --n-free 4 --rounds 5 --partitions 1000 --seed 1 \
  --max-steps 2000 --preset paper --field-schedule coupled \
  --field-scale 1.0 --normalization pairwise --quiet
```

The defaults use the standard bSB values `time_step=0.1`,
`pressure_slope=0.01`, `a0=1`, uniform initialization in `(-1,1)`, and an
RMS-based coupling normalization calibrated by `coupling_scale=0.3` against
Table 1. Every inferred parameter is exposed on the command line.
Position and momentum initialization can be controlled independently with
`--initial-scale` and `--momentum-scale` for sensitivity experiments.
`--restarts` controls full bSB trajectories. `--refine-restarts` adds cheap
random starts of the exact alternating pattern/type update from Theorem 3;
`--refine-sweeps` bounds each of these local refinements.
`--swap-restarts` enables a stronger hybrid refinement that reoptimizes both
column patterns after trial column-type swaps; `--swap-sweeps` bounds it.
`--anneal-restarts`, `--anneal-steps`, and `--anneal-temperature` enable a
Metropolis refinement over the two column-pattern vectors while keeping the
column-type vector conditionally optimal.
`--row-sa-restarts`, `--row-sa-steps`, and `--row-sa-temperature` provide a
dual row-pattern simulated-annealing refinement over the same decomposition
feasible set; no additional LUT capacity is introduced.
For biased Ising instances, `--field-schedule` controls the external-field
coefficient independently from the pairwise coupling (`coupled`,
`coupled-linear`, `linear`, or `sqrt`), and `--field-scale` scales it.
`--dynamics ballistic` is the paper mode; `--dynamics discrete` switches the
interaction MAC to sign-valued positions for a higher-accuracy dSB comparison.
`--normalization full` includes linear fields in the bSB RMS normalization;
`pairwise` retains the zero-field normalization for comparison.
In separate mode, `--separate-slack` optionally permits a small normalized ER
slack when using current MED as a deterministic secondary objective.

When Gurobi is found at build time, `--core-solver gurobi --oracle-time 0.1`
enables a bounded exact column-core oracle. This is a diagnostic/calibration
path; the default and reproduced solver remain `--core-solver bsb`.

## Figure 4: 16-input joint mode

The paper repeats its MSB-to-LSB joint optimization for `R` rounds and sets
`R=5` in the experimental configuration. For `n>10`, it also states that the
proposed solver uses DALTA's random input-partition strategy. This
implementation mirrors DALTA's `DrawVarOrd`: the first candidate is the
natural free set `{0,...,n_free-1}`, later candidates are obtained by
repeatedly shuffling the input variables, and duplicate partitions are
allowed. This detail materially affects the 16-input results; sampling unique
combinations uniformly is not the same experiment.

Run all ten Figure 4 cases, with atomic per-case checkpointing, using:

```sh
python3 experiments/run_figure4.py --jobs 3
```

The full 16-input Brent-Kung command is:

```sh
./build/ising-decomp \
  --input benchmarks/axbench-verilog/adder_8.txt \
  --mode joint \
  --n-free 7 \
  --rounds 5 \
  --partitions 1000 \
  --seed 1 \
  --max-steps 2000 \
  --preset paper
```

The completed records and comparison with the DALTA values printed in Figure 4
are stored under `results/figure4/`.
