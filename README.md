<div align="center">

# WFLO-GGA

**An Engineering Optimization Platform for Integrated Offshore Wind Farm**
**Layout and Electrical Cable Routing**

*Submitted to Applied Energy*

![MATLAB](https://img.shields.io/badge/MATLAB-R2018a%2B-blue?logo=mathworks)
![License](https://img.shields.io/badge/License-Academic-green)
![Algorithms](https://img.shields.io/badge/Algorithms-10-orange)
![Wind%20Farms](https://img.shields.io/badge/Benchmark%20Sites-8-teal)

</div>

---

## Overview

Offshore wind farm development involves two tightly coupled engineering decisions: **where to place turbines** within a geographically constrained site, and **how to connect them** via submarine cable networks. Optimizing these decisions jointly—rather than sequentially—can meaningfully reduce the Levelized Cost of Energy (LCOE), but requires a platform that faithfully models the physical, electrical, and economic realities of real projects.

This repository provides such a platform. It couples a physics-based wake and power model with a detailed CAPEX/OPEX cost model and three electrical cable routing strategies, grounded in a benchmark dataset of **8 real offshore wind farm sites**. Ten metaheuristic optimization algorithms are implemented as interchangeable solvers within this shared evaluation environment, enabling controlled, reproducible comparisons under consistent engineering assumptions.

---

## Benchmark Sites

<p align="center">
  <img src="figures/benchmark_layouts.png" alt="Benchmark wind farm boundaries and candidate layouts" width="90%">
  <br><em>Figure 1. Geographic boundaries and Poisson-disk-sampled candidate turbine positions for the 8 benchmark sites</em>
</p>

<p align="center">
  <img src="figures/benchmark_windrose.png" alt="Wind roses for the 8 benchmark sites" width="90%">
  <br><em>Figure 2. Directional wind speed distributions for the 8 benchmark sites</em>
</p>

---

## Optimized Layout Example

<p align="center">
  <img src="figures/final_layout.png" alt="Optimized turbine layout and cable routing" width="55%">
  <br><em>Figure 3. Optimized turbine placement and Balance-Sector cable routing produced by GGA on a representative site</em>
</p>

---

## Platform Architecture

The platform is organized as a layered evaluation pipeline. Optimization algorithms interact only with the top-level interface; all engineering models are encapsulated in the evaluation stack.

```
┌─────────────────────────────────────────────────────┐
│              Optimization Algorithm Layer           │
│   GGA · GA · AGA · BPSO · AGPSO · BDE · SaOFGDE     │
│              DOLSSA · RLPS_TLBO · EJAYA             │
└────────────────────────┬────────────────────────────┘
                         │  candidate layout (turbine indices)
                         ▼
┌─────────────────────────────────────────────────────┐
│               Cable Routing Layer                   │
│   cr_sector (BSR)  ·  cr_mst  ·  cr_sweep           │
│   → cable topology, segment lengths, current loads  │
└────────────────────────┬────────────────────────────┘
                         │  cable struct
                         ▼
┌─────────────────────────────────────────────────────┐
│               Objective Evaluation Layer            │
│                    evaluate.m                       │
│   Wake model → AEP → CAPEX + OPEX → LCOE            │
└─────────────────────────────────────────────────────┘
```

Each algorithm receives a single scalar fitness value (LCOE) and has no direct access to the physical or electrical models, ensuring all comparisons reflect algorithmic performance under identical engineering conditions.

---

## Physical & Economic Model

### Wake Model

Turbine interactions are computed using the **Jensen (top-hat) wake model**:

- Thrust coefficient: *C*_T = 0.8
- Wake decay constant: *k* = 0.04 (offshore)
- Multiple wake superposition via quadratic velocity deficit addition

### Power Model

Turbine output is interpolated from a **manufacturer power curve** (cut-in 3 m/s, rated 4.2 MW at 13 m/s, cut-out 25 m/s) accounting for wake-reduced inflow speed at each turbine for each wind direction sector.

### Cost Model (25-year LCOE)

| Component | Model |
|:----------|:------|
| Turbine CAPEX | Unit cost × rated capacity |
| Foundation & installation | Uniform per-turbine cost (depth-independent) |
| Inter-array cables | Length × unit cost, typed by current load (3 grades) |
| Export cable | Fixed per-site |
| O&M | Annual percentage of CAPEX |
| Decommissioning | Terminal lump sum |
| Discount rate | 5% |

> Foundation costs are modeled as spatially uniform to ensure cross-site algorithmic comparability. The framework accepts depth-dependent cost functions when bathymetric data are available.

### Annual Energy Production

AEP integrates over the site wind rose (16 directional sectors × wind speed bins), applying the wake-reduced power curve at each turbine:

```
AEP = Σ_d Σ_v  f(d,v) · P_wake(layout, d, v) · 8760
```

---

## Electrical Cable Routing

Three routing strategies are provided as drop-in modules, all subject to the same **current-capacity constraints** (33 kV inter-array voltage, three cable grades):

| Strategy | Description | Characteristic |
|:---------|:------------|:---------------|
| **Balance-Sector (BSR)** | Partitions turbines into angular sectors relative to the substation; builds a capacity-constrained MST per sector | Spatial locality, low imbalance |
| **MST** | Global minimum spanning tree with capacity-constraint repair | Minimum total cable length |
| **Sweep** | Sweep-line sector assignment with sequential cable loading | Fast, geometry-aware |

Routing is called once per candidate layout evaluation and returns the full cable topology, segment lengths, cable grades, and total cable CAPEX.

---

## Benchmark Dataset

Eight real offshore wind farms are included, spanning four countries and a turbine count range of 20–175:

| # | Site | Country | Turbines | Characteristic |
|:-:|:-----|:-------:|:--------:|:---------------|
| 1 | Zhuhai Guishan Hai | China | 34 | Irregular coastal boundary |
| 2 | Egmond aan Zee | Netherlands | 36 | Compact rectangular layout |
| 3 | Shanghai Lingang | China | 56 | Complex multi-polygon boundary |
| 4 | Prinses Amaliawindpark | Netherlands | 61 | Regular offshore grid |
| 5 | Nysted | Denmark | 72 | Large rectangular array |
| 6 | Sheringham Shoal | UK | 88 | Irregular offshore boundary |
| 7 | Rodsand II | Denmark | 90 | Curved elongated boundary |
| 8 | London Array | UK | 175 | Large-scale complex site |

**Site data includes**: geographic boundary polygons, directional wind speed distributions (16 sectors), and site-specific electrical parameters. Turbine candidate positions are generated via **Poisson-disk sampling** to ensure minimum spacing compliance while providing uniform spatial coverage.

---

## Algorithm Suite

Ten metaheuristic algorithms are implemented as interchangeable solvers. All use the same chromosome encoding (a vector of *M* distinct candidate-point indices), the same evaluation pipeline, and the same initial population when provided.

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

---

## Repository Structure

```
WFLO-GGA/
├── main.m                      # Experiment orchestration entry point
├── alg/                        # Optimization algorithm implementations
│   ├── GGA.m                   # ★ Proposed method
│   ├── GA.m / AGA.m            # Genetic algorithm variants
│   ├── BPSO.m / AGPSO.m        # Swarm intelligence variants
│   ├── BDE.m / SaOFGDE.m       # Differential evolution variants
│   ├── DOLSSA.m                # Sparrow search variant
│   ├── RLPS_TLBO.m             # Teaching-learning variant
│   └── EJAYA.m                 # Jaya variant
├── utils/
│   ├── evaluate.m              # LCOE / AEP / CF evaluation (core model)
│   ├── cr_sector.m             # Balance-Sector Routing (BSR)
│   ├── cr_mst.m                # Minimum Spanning Tree routing
│   ├── cr_sweep.m              # Sweep-line routing
│   ├── load_problem_poisson.m  # Site loader + Poisson-disk sampling
│   ├── load_layout.m           # Lat/lon to Cartesian projection
│   └── unique_fix.m            # Chromosome feasibility repair
├── figures/                    # README figures
└── data/
    └── OWF8.qgz                # Compressed benchmark dataset (8 sites)
```

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
% Choose cable routing strategy
cfg.routing_fn = @cr_sector;   % Balance-Sector Routing (default)
% cfg.routing_fn = @cr_mst;    % Minimum Spanning Tree
% cfg.routing_fn = @cr_sweep;  % Sweep-line

% Choose algorithms
cfg.algorithms = { @GGA, @GA, @BDE };
cfg.algonames  = { 'GGA', 'GA', 'BDE' };

% Choose wind farm sites
cfg.case_list = { 'Denmark_Nysted', 'UK_London_Array' };

% Set run parameters
cfg.popsize  = 30;    % population size
cfg.max_it   = 500;   % iterations per run
cfg.runTime  = 10;    % independent runs per case
cfg.base_seed = 42;   % random seed for reproducibility
```

**4. Run**

```matlab
main
```

Results are saved to `results/results_YYYYMMDD_HHMMSS/`.

---

## Output Structure

```
results/
└── results_YYYYMMDD_HHMMSS/
    ├── global_run_summary.csv        <- LCOE across all cases and algorithms
    ├── metadata/
    │   ├── config_snapshot.mat       <- full experiment configuration
    │   └── environment_info.mat      <- MATLAB version, hardware info
    └── {site_name}/
        └── {algorithm}/
            ├── run_summary.csv
            ├── run_N_summary.mat     <- best LCOE, AEP, CF per run
            └── {algorithm}_runN.mat <- complete generation history:
                                        population, fitness, AEP, CF,
                                        cable topology, timing breakdown
```

---

## Citation

If you use this platform, dataset, or any algorithm implementation in your research, please cite:


---

## License

Released for academic research purposes. Please contact the authors prior to any commercial use of the benchmark dataset or source code.
