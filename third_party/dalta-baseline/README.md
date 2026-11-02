# DALTA baseline snapshot

This directory contains the minimum DALTA solver sources needed by the
optional Table 1 baseline. They were taken from
[`changmg/DALTA`](https://github.com/changmg/DALTA) at commit `f8458f7` and
adapted for this reproduction:

- machine-specific Gurobi paths were replaced by CMake discovery;
- separate/joint and heuristic/ILP modes were exposed through `--solver`;
- the separate ILP received an algebraically equivalent formulation for
  newer Gurobi versions;
- MED and runtime are printed in a machine-readable form.

The Ising solver does not compile or link these files. They are provided only
for baseline comparison and are built when `BUILD_DALTA_BASELINE=ON`.

The upstream DALTA snapshot does not contain a project-wide license file.
Before redistributing this directory in a public repository, obtain or confirm
redistribution permission from the DALTA authors. The embedded
`cmdline.hpp` retains its original BSD license notice.
