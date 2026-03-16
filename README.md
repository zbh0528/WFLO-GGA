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

## Highlights

| | |
|---|---|
| **Algorithms** | 10 metaheuristic methods (genetic, swarm, DE, learning-based) |
| **Cable routing** | 3 strategies: Balance-Sector (BSR), MST, Sweep |
| **Benchmark** | 8 real offshore wind farms across 4 countries |
| **Reproducibility** | Fixed random seeds, full generation history, config snapshots |
| **Parallelism** | Optional multi-run parallel execution via  |

---

## Repository Structure



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
| 1 | Zhuhai Guishan Hai | 🇨🇳 China | 20 | Small |
| 2 | Egmond aan Zee | 🇳🇱 Netherlands | 36 | Small |
| 3 | Shanghai Lingang | 🇨🇳 China | 80 | Large |
| 4 | Prinses Amaliawindpark | 🇳🇱 Netherlands | 60 | Medium |
| 5 | Nysted | 🇩🇰 Denmark | 72 | Medium |
| 6 | Sheringham Shoal | 🇬🇧 UK | 88 | Large |
| 7 | Rodsand II | 🇩🇰 Denmark | 90 | Large |
| 8 | London Array | 🇬🇧 UK | 175 | Extra Large |

> Turbine candidate positions are generated via **Poisson-disk sampling** constrained to each site boundary.

---

## Quick Start

### 1. Clone and open in MATLAB



### 2. Select routing strategy



### 3. Select algorithms



---

## Configuration Reference

| Parameter | Default | Description |
|:----------|:-------:|:------------|
|  |  | Population size |
|  |  | Maximum iterations per run |
|  |  | Number of independent runs |
|  |  | Base random seed |
|  |  | Enable  parallel runs |
|  |  | Cable routing function handle |

---

## Output Structure



Each  file stores the complete  struct array with population, fitness, AEP, CF, cable topology, and per-generation timing.

---

## Citation

If you use this code or benchmark dataset, please cite:

Need exactly one file argument

---

## License

Released for academic research purposes. Please contact the authors prior to any commercial use of the benchmark dataset or source code.
