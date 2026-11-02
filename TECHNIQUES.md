# Technique inventory

All presets use `n-free=4`, `n-bound=5`, a 16-by-32 free-row Boolean matrix,
and the same 64-spin column core for the 9-input Table-1 experiments.

| Technique | Target paper | `paper` | `paper-plus` |
|---|---:|---:|---:|
| Column-based disjoint decomposition | yes | yes | yes |
| Separate ER and joint MED Ising formulations | yes | yes | yes |
| Ballistic simulated bifurcation | yes | yes | yes |
| Euler integration | yes | yes | yes |
| Symplectic-Euler ordering and inelastic walls | bSB reference | yes | yes |
| Energy-variance dynamic stopping | yes | yes | yes |
| Exact conditional update and feedback of `T` | yes | yes | yes |
| Exhaustive 126 partitions for 9 inputs | yes | yes | yes |
| Five MSB-to-LSB joint rounds | yes | yes | yes |
| Independent field `eta(t)*h` | not specified | one per mode | 2 separate / 3 joint |
| `dt`, pressure slope, integration cap | not reported | inferred | inferred |
| Random position/momentum initialization | not reported | one trajectory | one seed, 2/3 trajectories |
| Coupling/field normalization | not reported | full RMS | full plus joint pairwise |
| Multiple fixed random seeds | not reported | no | no |
| Final-LUT MED selection across agents | not reported | no | yes |
| ER slack with MED secondary objective | no | no | no |
| Parallel execution of agents | not reported | no | yes |

The target paper includes bias terms in its Ising Hamiltonian but does not
specify how they enter bSB dynamics. Related SB work by authors in the same
group uses `c0*J*x + eta(t)*h`; `paper` therefore uses one fixed
calibration for each objective: `eta(t)=0.7*sqrt(p(t))` for separate ER and
`eta(t)=0.329*c0*p(t)` for joint MED. There is no runtime schedule selection.
These are explicit reconstruction assumptions, not a claim that the target
paper used those schedules. Their scales, initialization, coupling
normalization, integration cap, and random seed likewise remain inferred.

`paper` is therefore best read as "paper-technique faithful": it uses
all algorithmic ideas stated in the target paper and only one deterministic
choice for each omitted numerical detail. `paper-plus` isolates the
smallest tested enhancement that addresses schedule sensitivity: a fixed,
single-seed field portfolio and final MED selection. It is the only retained
quality-oriented extension, and those two additions are deliberately outside
the reported target algorithm.

The following diagnostic features are implemented but disabled in both
presets: auxiliary-spin bias embedding, dSB dynamics, alternating local
refinement, type-swap search, Metropolis column-pattern search, row-pattern SA,
and the optional Gurobi core oracle.
