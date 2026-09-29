# Interfacial-tension operators in colour-gradient lattice Boltzmann models

Code, simulation results and analysis scripts for the manuscript

> *Interfacial-Tension Operators in Color-Gradient Lattice Boltzmann Models* (in preparation).

The study compares seven ways of imposing interfacial tension in a colour-gradient lattice
Boltzmann model. Everything else is kept identical, so that only the tension operator changes
between runs:

- D3Q19 lattice with one node in *z*;
- single-relaxation-time (BGK) collision with Guo forcing;
- Latva-Kokko–Rothman recolouring (β = 0.8);
- the same isotropic gradient stencil, evaluated with mediators.

## Operators

| Paper label | Code (`operador.txt`) | Description |
|---|---|---|
| PG | `P1` | even perturbation, weighted Gunstensen form (λ = 0) |
| PL, RP | `P3 -0.4` | Liu–Valocchi–Kang family, χ = −2/5 (Reis–Phillips on the D3Q19 slab) |
| PL, LKR | `P3 0.5` | same family, χ = 1/2 (weighted Latva-Kokko–Rothman) |
| PL, LVK | `P3 2` | same family, χ = 2 (Liu–Valocchi–Kang) |
| F1 | `F1` | continuum-surface force, F = (σ/2) κ ∇ρᴺ |
| F2 | `F2` | divergence of the capillary stress tensor |
| F3 | `F3 1.569988` | chemical-potential force μ∇ρᴺ (the number is the interface width ξ) |

Any other value of χ can be given to `P3`. All formulas are in `Operadores_tensao.cpp`.

## Repository layout

| Path | Content |
|---|---|
| `Operadores_tensao.cpp` | the tension operators, the common collision step and the interface Fourier modes |
| `Imm_OPS.cpp` | drop program: static or translating cylindrical drop in a periodic box |
| `Imm_INCL.cpp` | straight interface inclined at angle θ = atan(q/p) in a periodic box |
| `Galilean_diagnostics.cpp` | diagnostics in the frame of the drop (mean flow, slip, negative populations) |
| `verifica_operadores.cpp` | analytical checks of the operators (series E0); output in `E0_resultado.txt` |
| `SimBoltz_Functions/` | the subset of the SimBoltz library used by these programs (lattice, BGK, recolouring, mediators, I/O) |
| `rodar_lote.py` | builds and runs the batches E1–E4v in parallel on the local machine |
| `analisar_lote.py` | collects the batch results into `Lotes/*.csv` and `Lotes/RESULTADOS.md` |
| `roda_caso.sh`, `ger_caixa.py` | run a single case; generate the geometry file `meio.vtk` |
| `Lotes/` | per-run outputs of the 619 simulations of the paper (see below) |
| `simulation_data/` | the run-level tables used in the manuscript, including the follow-up and control runs |
| `paper_scripts/` | `tabelas.py` and `figuras.py`, which produce the tables and figures of the manuscript |

Comments and variable names in the source code are in Portuguese.

## Requirements

- A C++17 compiler with OpenMP (tested with g++ 13).
- Python 3. The batch and analysis scripts use only the standard library; `paper_scripts/` needs
  `numpy` and `matplotlib`.

## Building

```bash
make              # Imm_OPS, Imm_INCL and verifica_operadores
make verifica     # runs the operator checks and writes E0_resultado.txt
```

The default flags use `-O3 -march=native`. A different compiler or set of flags changes the
results only at the level of round-off; bitwise reproduction requires the same compiler and flags.

## Running a single case

Each case is a folder containing these input files:

- `data_in.txt`: the usual SimBoltz input. `fat_R_B` is the amplitude A for the perturbations and
  the tension σ for the forces.
- `operador.txt`: one line, for example `P3 0.5` or `F2`.
- `velocity.txt`: the initial translation velocity v0.
- `meio.vtk`: the geometry, created by `ger_caixa.py`.
- `inclinacao.txt`: `p q`, only for `Imm_INCL`.

For the perturbations, the tension is σ = (4/9) A τ.

```bash
./roda_caso.sh teste/LKR "P3 0.5" 0.225 1.0 0.8 30000 100 0
```

## Reproducing the batches

```bash
python3 rodar_lote.py todos --procs 6 --lista   # list the cases and the estimated time
python3 rodar_lote.py E2 --procs 6              # run one batch; interrupted runs resume
python3 analisar_lote.py                        # rebuild Lotes/*.csv and Lotes/RESULTADOS.md
```

| Batch | Content | Runs |
|---|---|---|
| E1 | Laplace law, R = 10–50 | 49 |
| E2 | static drop, τ × σ map | 168 |
| E2chi | fine scan of χ at τ = 0.6 and 1 | 30 |
| E3 | straight interface at six orientations | 42 |
| E3chi | χ scan at θ = 26.6° | 15 |
| E4 | stability map, τ = 0.51–2, σ = 0.1–2 | 280 |
| E4v | translating drop, v0 = 0.05–0.25 | 35 |

In `Lotes/`, each case keeps these outputs:

- `resumo.dat`: final state and second-half averages;
- `laplace.dat`: time series of radius, pressures, tension and velocity;
- `galilean.dat`: diagnostics in the drop frame;
- `iso.dat`: Fourier modes of the interface;
- `inclinada.dat` and `perfil.dat`: the corresponding outputs of the inclined-interface program.

The final field files (`*.vtk`) are not included because of their size, and neither is `meio.vtk`,
which `ger_caixa.py` regenerates.

## Data used in the manuscript

`simulation_data/` contains the tables from which the manuscript's tables and figures are built.
Its `README.md` describes each file. Two points differ from the raw `Lotes/` tables:

- The F1 and F2 runs at τ = 0.6, σ = 0.1 of series E2 were not converged after 38 000 steps.
  In `simulation_data/E2.csv` they are replaced by 100 000-step runs; `Lotes/` keeps the originals.
- The E4v runs were made with a build of `Imm_OPS` that did not write the net force. The net force
  and the velocity relative to the mean flow were derived from the time histories and are stored in
  `simulation_data/E4v_derived.csv` and `.json`.

To regenerate the tables and figures:

```bash
cd paper_scripts
python3 tabelas.py ../simulation_data ../simulation_data/E4v_derived.json
python3 figuras.py ../simulation_data ../simulation_data/E4v_derived.json
```

## Citation

If you use this code or data, please cite the manuscript (reference to be added on publication).
See `CITATION.cff`.

## License

Copyright (C) 2026 the authors (see `CITATION.cff`).

This program is free software: you can redistribute it and/or modify it under the terms of the
GNU General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version. It is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See `LICENSE` for the full text.
