# WFLO-GGA

**A Geometry-guided Genetic Algorithm for Integrated Offshore Wind Farm Layout and Electrical Cable Routing Optimization**

This repository provides the MATLAB implementation for the paper submitted to *Applied Energy*. The framework benchmarks 10 metaheuristic algorithms on 8 real offshore wind farm sites, with three interchangeable cable routing strategies and a unified LCOE-based objective.

---

## Features

- **10 algorithms** - GA, GGA, AGA, BPSO, AGPSO, BDE, SaOFGDE, DOLSSA, RLPS_TLBO, EJAYA
- **3 cable routing strategies** - Sector-based (BSR), Minimum Spanning Tree, Sweep
- **8 real wind farm benchmark sites** across Denmark, Netherlands, UK, and China
- Unified LCOE objective with Jensen wake model and 25-year economic model
- Fully reproducible: deterministic seeding, config snapshots, generation history saved per run
- Optional parallel execution via MATLAB Parallel Computing Toolbox

---

## Repository Structure

```
WFLO-GGA/
├── main.m                      # Main entry point
├── alg/                        # Algorithm implementations
│   ├── GGA.m                   # Geometry-guided Genetic Algorithm (proposed)
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
│   ├── load_problem_poisson.m  # Problem loader (Poisson-disk candidate generation)
│   ├── evaluate.m              # LCOE / AEP / CF evaluation
│   ├── cr_sector.m             # Balance-Sector Routing (BSR)
│   ├── cr_mst.m                # Minimum Spanning Tree routing
│   ├── cr_sweep.m              # Sweep-line routing
│   ├── load_layout.m           # Lat/lon to Cartesian coordinate conversion
│   └── unique_fix.m            # Index constraint repair
└── data/
    └── OWF8.qgz                # Compressed wind farm dataset
```

---

## Algorithms

| Algorithm | Category | Key Mechanism |
|-----------|----------|---------------|
| **GGA** | Genetic | Geometry-guided crossover using angular sectors |
| GA | Genetic | Tournament selection, single-point crossover, elitism |
| AGA | Genetic | Adaptive mutation/crossover rates |
| BPSO | Swarm | Binary PSO with top-M velocity ranking |
| AGPSO | Hybrid | GA + PSO with stagnation recovery |
| BDE | Differential Evolution | Ranking-based pBest/pWorst mutation |
| SaOFGDE | Differential Evolution | Fractional-order historical memory + adaptive CR |
| DOLSSA | Swarm | Opposition-based learning Sparrow Search |
| RLPS_TLBO | Learning-based | Q-learning phase selection in TLBO |
| EJAYA | Learning-based | Enhanced Jaya with attraction/repulsion dynamics |

All algorithms share an identical function signature and evaluation pipeline, enabling controlled comparison.

---

## Benchmark Wind Farms

| # | Site | Country | Turbines |
|---|------|---------|----------|
| 1 | Zhuhai Guishan Hai | China | 20 |
| 2 | Egmond aan Zee | Netherlands | 36 |
| 3 | Shanghai Lingang | China | 80 |
| 4 | Prinses Amaliawindpark | Netherlands | 60 |
| 5 | Nysted | Denmark | 72 |
| 6 | Sheringham Shoal | UK | 88 |
| 7 | Rodsand II | Denmark | 90 |
| 8 | London Array | UK | 175 |

Candidate turbine positions are generated via Poisson-disk sampling constrained to each site geographic boundary.

---

## Quick Start

**Requirements**: MATLAB R2018a or later. No external toolboxes required for core functionality; Parallel Computing Toolbox is optional.

```matlab
% Open MATLAB, set project root as working directory, then run:
main
```

**Select routing strategy** in `main.m`:

```matlab
cfg.routing_fn = @cr_sector;   % Balance-Sector Routing (default)
% cfg.routing_fn = @cr_mst;    % Minimum Spanning Tree
% cfg.routing_fn = @cr_sweep;  % Sweep-line
```

**Select algorithms** in `main.m`:

```matlab
cfg.algorithms = { @GGA, @GA, @BDE };
cfg.algonames  = { 'GGA', 'GA', 'BDE' };
```

---

## Configuration

Key parameters in `main.m`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `cfg.popsize` | 30 | Population size |
| `cfg.max_it` | 100 | Maximum iterations per run |
| `cfg.runTime` | 1 | Number of independent runs |
| `cfg.base_seed` | 20260316 | Base random seed for reproducibility |
| `cfg.enable_parallel` | false | Enable parallel runs via parfor |
| `cfg.routing_fn` | @cr_sector | Cable routing strategy |

---

## Output Structure

Each experiment session generates a timestamped directory:

```
results/
└── results_YYYYMMDD_HHMMSS/
    ├── global_run_summary.csv
    ├── metadata/
    │   ├── config_snapshot.mat
    │   └── environment_info.mat
    └── {case_name}/
        └── {algorithm_name}/
            ├── run_summary.csv
            ├── run_1_summary.mat
            └── {algorithm}_run01.mat
```

Each .mat file stores the complete generations cell array, including population, fitness, AEP, CF, cable configurations, and per-generation timing.

---

## Citation

If you use this code or benchmark dataset in your research, please cite:


---

## License

This repository is released for academic research purposes. Please contact the authors before using the benchmark dataset or code in commercial applications.
