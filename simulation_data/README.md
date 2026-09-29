# Simulation data used in the manuscript

These comma-separated tables contain the run-level summaries behind the manuscript. They were
produced by `analisar_lote.py` from the runs in `../Lotes/` on 2026-09-28. The series contain E1: 49, E2: 168, E2chi: 30, E3: 42, E3chi: 15,
E4: 280 and E4v: 35 runs (619 total).

- `E1.csv`: Laplace-law radius series and tension fits.
- `E2.csv`: static-drop spurious-current map over relaxation time and nominal tension. The
  F1 and F2 rows at $\tau=0.6$, input $\sigma=0.1$ contain converged 100,000-step follow-up
  values in place of their original 38,000-step values.
- `E2_superseded.csv`: those two original nonconverged rows, retained for audit.
- `E2chi.csv`: fifteen-member higher-moment scan at two relaxation times.
- `E3.csv`, `E3chi.csv`: straight-interface orientation and chi scans.
- `E4.csv`: fixed-horizon stability map.
- `E4v.csv`: translating-drop run summaries. The `umax_rel_media` and `soma_Fx_med` columns
  in the original batch-analysis CSV are `nan` because those quantities were not written to
  `resumo.dat`; use `E4v_derived.csv` for them.
- `E4v_derived.csv`: velocity disturbance relative to domain-mean velocity, net force and
  population diagnostics extracted from the time histories for the 35 published operators/runs.
- `E4v/<operator>/v0=<speed>/galilean.dat`: the 35 translating-drop time histories used to
  verify mean-relative velocities, phase slip and domain-mean momentum change.

`simulation_data/control_runs.csv` records the five targeted follow-up simulations and the
two original, nonconverged F1/F2 runs. The underlying time histories are included under
`simulation_data/control_timeseries/`.

The manuscript's tables and figures are produced from these files by the scripts in
`../paper_scripts/`:
`python3 tabelas.py ../simulation_data ../simulation_data/E4v_derived.json` and
`python3 figuras.py ../simulation_data ../simulation_data/E4v_derived.json`.

Units are lattice units. `sigma_nominal` is the model input, `sigma_med` is the Laplace-law
measurement, and `umax_sigma`/`Ca_s` use the measured tension for normalization. Numerical
stability denotes only the stated finite run horizon and velocity criterion; negative color
populations can occur even when total populations remain positive. The `E4v_derived.csv`
column `net_force_x_per_step` is the domain sum of the force at each step, averaged over the
second half of the run. The source code is in the root of this repository and the per-run
outputs, with their input files, are in `../Lotes/`.
