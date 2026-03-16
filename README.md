<div align="center">

# WFLO-GGA

**A Geometry-guided Genetic Algorithm for Integrated Offshore Wind Farm**
**Layout and Electrical Cable Routing Optimization**

*Submitted to Applied Energy*

![MATLAB](https://img.shields.io/badge/MATLAB-R2018a%2B-blue?logo=mathworks)
![License](https://img.shields.io/badge/License-Academic-green)
![Algorithms](https://img.shields.io/badge/Algorithms-10-orange)
![Wind%20Farms](https://img.shields.io/badge/Benchmark%20Sites-8-teal)

</div>

---

## Overview

This repository provides a unified MATLAB framework for benchmarking metaheuristic algorithms on the **integrated wind farm layout and cable routing optimization** problem. All algorithms minimize **Levelized Cost of Energy (LCOE)** over a 25-year horizon using the Jensen wake model, evaluated on 8 real offshore wind farm sites with three interchangeable cable routing strategies.

---

## Benchmark Sites

<p align="center">
  <img src="figures/benchmark_layouts.png" alt="Benchmark wind farm boundaries and layouts" width="90%">
  <br><em>Figure 1. Boundaries and candidate turbine positions for the 8 benchmark wind farm sites</em>
</p>

<p align="center">
  <img src="figures/benchmark_windrose.png" alt="Benchmark wind roses" width="90%">
  <br><em>Figure 2. Wind roses for the 8 benchmark sites</em>
</p>

---

## Sample Result

<p align="center">
  <img src="figures/final_layout.png" alt="Optimized turbine layout and cable routing" width="55%">
  <br><em>Figure 3. Example optimized turbine layout and cable routing produced by GGA</em>
</p>

---

## Highlights

| | |
|---|---|
| **Algorithms** | 10 metaheuristic methods (genetic, swarm, DE, learning-based) |
| **Cable routing** | 3 strategies: Balance-Sector (BSR), MST, Sweep |
| **Benchmark** | 8 real offshore wind farms across 4 countries |
| **Reproducibility** | Fixed random seeds, full generation history, config snapshots |
| **Parallelism** | Optional multi-run parallel execution via `parfor` |

---

## Repository Structure

```
WFLO-GGA/
├── main.m                      # Main entry point & experiment orchestration
├── alg/                        # Algorithm implementations
│   ├── GGA.m                   # ★ Geometry-guided GA (proposed)
│   ├── GA.m                    # Classical Genetic Algorithm
│   ├── AGA.m                   # Adaptive Genetic Algorithm
│   ├── BPSO.m                  # Binary Particle Swarm Optimization
│   ├── AGPSO.m                 # Adaptive Genetic PSO
│   ├── BDE.m                   # Ranking/Half-Split Differential Evolution
│   ├── SaOFGDE.m               # Self-adaptive Fractional-order GDE
│   ├── DOLSSA.m                # Opposition-based Sparrow Search
│   ├── RLPS_TLBO.m             # RL Phase-Selection TLBO
│   └── EJAYA.m                 # Enhanced Jaya Algorithm
├── utils/
│   ├── load_problem_poisson.m  # Problem loader (Poisson-disk sampling)
│   ├── evaluate.m              # LCOE / AEP / CF evaluation
│   ├── cr_sector.m             # Balance-Sector Routing (BSR)
│   ├── cr_mst.m                # Minimum Spanning Tree routing
│   ├── cr_sweep.m              # Sweep-line routing
│   ├── load_layout.m           # Lat/lon to Cartesian projection
│   └── unique_fix.m            # Chromosome constraint repair
├── figures/                    # README figures
└── data/
    └── OWF8.qgz                # Compressed wind farm dataset (8 sites)
```

---

## Algorithms

| Algorithm | Category | Key Mechanism |
|:----------|:---------|:--------------|
| **GGA** ★ | Genetic | Geometry-guided crossover via angular sector partitioning |
| GA | Genetic | Tournament selection · single-point crossover · elitism |
| AGA | Genetic | Adaptive mutation and crossover rates |
| BPSO | Swarm | Binary PSO with top-M velocity ranking |
| AGPSO | Hybrid | GA + PSO with stagnation recovery |
| BDE | Differential Evolution | Ranking-based pBest / pWorst mutation |
| SaOFGDE | Differential Evolution | Fractional-order historical memory · adaptive CR |
| DOLSSA | Swarm | Opposition-based learning Sparrow Search |
| RLPS_TLBO | Learning-based | Q-learning phase selection in TLBO |
| EJAYA | Learning-based | Enhanced Jaya with attraction / repulsion dynamics |

> All algorithms share an identical function signature and evaluation pipeline for fair comparison.

---

## Benchmark Wind Farms

| # | Site | Country | Turbines | Scale |
|:-:|:-----|:-------:|:--------:|:-----:|
| 1 | Zhuhai Guishan Hai | China | 20 | Small |
| 2 | Egmond aan Zee | Netherlands | 36 | Small |
| 3 | Shanghai Lingang | China | 80 | Large |
| 4 | Prinses Amaliawindpark | Netherlands | 60 | Medium |
| 5 | Nysted | Denmark | 72 | Medium |
| 6 | Sheringham Shoal | UK | 88 | Large |
| 7 | Rodsand II | Denmark | 90 | Large |
| 8 | London Array | UK | 175 | Extra Large |

> Turbine candidate positions are generated via **Poisson-disk sampling** constrained to each site boundary.

---

## Quick Start

### Requirements

- MATLAB R2018a or later
- No external toolboxes required for core functionality
- Parallel Computing Toolbox (optional, for multi-run parallelism)

### Steps

**1. Clone the repository**

```bash
git clone https://github.com/zbh0528/WFLO-GGA.git
cd WFLO-GGA
```

**2. Open MATLAB and add to path**

```matlab
addpath(genpath(pwd));
```

**3. Configure experiment in `main.m`**

```matlab
% Select cable routing strategy
cfg.routing_fn = @cr_sector;   % Balance-Sector Routing (default)
% cfg.routing_fn = @cr_mst;    % Minimum Spanning Tree
% cfg.routing_fn = @cr_sweep;  % Sweep-line

% Select algorithms to benchmark
cfg.algorithms = { @GGA, @GA, @BDE };
cfg.algonames  = { 'GGA', 'GA', 'BDE' };

% Select wind farm cases
cfg.case_list = { 'Denmark_Nysted', 'UK_London_Array' };

% Adjust run parameters
cfg.popsize  = 30;
cfg.max_it   = 500;
cfg.runTime  = 10;
```

**4. Run**

```matlab
main
```

Results are saved automatically to `results/results_YYYYMMDD_HHMMSS/`.

---

## Configuration Reference

| Parameter | Default | Description |
|:----------|:-------:|:------------|
| `cfg.popsize` | `30` | Population size |
| `cfg.max_it` | `100` | Maximum iterations per run |
| `cfg.runTime` | `1` | Number of independent runs |
| `cfg.base_seed` | `20260316` | Base random seed |
| `cfg.enable_parallel` | `false` | Enable `parfor` parallel runs |
| `cfg.routing_fn` | `@cr_sector` | Cable routing function handle |

---

## Output Structure

```
results/
└── results_YYYYMMDD_HHMMSS/
    ├── global_run_summary.csv        <- aggregated LCOE across all cases
    ├── metadata/
    │   ├── config_snapshot.mat       <- full configuration
    │   └── environment_info.mat      <- MATLAB version, hardware
    └── {case_name}/
        └── {algorithm_name}/
            ├── run_summary.csv
            ├── run_N_summary.mat     <- best fitness, success flag
            └── {algorithm}_runN.mat <- full generation history
```

Each `.mat` file stores the complete `generations` struct array with population, fitness, AEP, CF, cable topology, and per-generation timing.

---

## Citation

If you use this code or benchmark dataset, please cite:

```bibtex
@article{wflo_gga_2026,
  title   = {A Geometry-guided Genetic Algorithm for Integrated Offshore Wind Farm
             Layout and Electrical Cable Routing Optimization},
  journal = {Applied Energy},
  year    = {2026},
  note    = {Under review}
}
```

---

## License

Released for academic research purposes. Please contact the authors prior to any commercial use of the benchmark dataset or source code.
