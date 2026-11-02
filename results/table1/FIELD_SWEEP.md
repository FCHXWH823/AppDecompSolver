# Independent-field sweep

The `paper` calibration varies only the independent external field
`eta(t)*h`. Every run retains one bSB trajectory, seed 1, full normalization,
the free-row 64-spin core, and disables portfolio selection and ER slack.

Initial six-function averages were:

| Schedule | Scale | Separate MED | Joint MED |
|---|---:|---:|---:|
| `linear` | 0.50 | 7.469 | 3.201 |
| `linear` | 0.70 | 7.434 | 2.729 |
| `linear` | 0.90 | 8.056 | **2.687** |
| `linear` | 1.00 | 8.030 | 2.707 |
| `linear` | 1.20 | 8.030 | 2.895 |
| `sqrt` | 0.50 | 7.470 | 2.972 |
| `sqrt` | 0.70 | **7.420** | 2.711 |
| `sqrt` | 1.00 | 8.024 | 2.747 |

A finer joint-only linear scan at scales 0.75, 0.80, 0.85, 0.95, 1.05, and
1.10 did not improve on scale 0.90. A subsequent coupled-linear scan found a
better joint calibration at scale 0.329. Its six-function average MED is
2.5140, while the mean of the paper's six unrounded table entries is 2.5133.
The preset therefore makes a fixed, mode-dependent choice rather than running
a portfolio: `0.7*sqrt(p(t))` for separate ER and `0.329*c0*p(t)` for joint
MED.

## Gap diagnosis and paper-plus

The strict separate gap is not caused by too few bSB restarts. With the same
configuration, one, two, four, and eight ER-selected restarts produced average
MED 7.420, 8.098, 8.367, and 8.798, respectively. Better ER search can worsen
MED because separate mode does not weight output bits by numeric significance.
Scans over seeds 1--10 across all functions, seeds 1--100 for separate `ln`,
initialization scales, coupling scales, pressure ramps, integration caps, tie
policies, and five separate rounds did not remove the strict single-setting
gap. The best `ln` seed in that diagnostic sweep was seed 85 with MED 6.881;
the retained presets both use seed 1.

A second, fixed coupled-linear field trajectory reaches a different separate
basin and reduces `erf` from 7.381 to 4.160. For joint, fixed linear and
coupled/pairwise trajectories cover the remaining `exp` and `erf` basins.
`paper-plus` runs these candidates with the same seed and retains the
complete LUT with minimum MED. This final selection is a documented extension,
not a technique claimed by the target paper.
