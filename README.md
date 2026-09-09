# Code for "Neglecting movement-disease feedback biases wildlife epidemic forecasts"

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

### Contents

- [Overview](#overview)
- [Algorithm Description](#algorithm-description)
- [Repo Contents](#repo-contents)
- [System Requirements](#system-requirements)
- [Installation Guide](#installation-guide)
- [Demo](#demo)
- [Instructions for Use](#instructions-for-use)
- [License](#license)
- [Citation](#citation)

# Overview

This repository contains the MATLAB simulation code accompanying the
manuscript *"Neglecting movement-disease feedback biases wildlife epidemic
forecasts."* The code implements a stochastic, individual-based SIR/SI
epidemic model on a 2-D grid in which agent movement (diffusion + drift
toward a home range, optionally with group cohesion) is coupled to
infection status, so infected and susceptible agents move differently and
deposit/respond to an environmental pathogen signal left on the grid. Two
model variants are provided: a periodic-boundary ("torus") model and a
hard-edge model with a group-cohesion term.

Each run generates its own simulated dataset (agent trajectories and SIR
counts); no external or proprietary data is required.

# Algorithm Description

A complete, detailed description of the code's functionality is provided
in the manuscript's **Methods section**: the movement model (Eq. 1,
Euler-Maruyama discretization of the drift-diffusion process), the
group-cohesion extension (Eq. 3b), the environmental pathogen
decay/transmission rule, and the initialization, burn-in, and replicate
procedure are all specified there in full mathematical/procedural detail.

# Repo Contents

- [Simulation_Framework_Main.m](Simulation_Framework_Main.m): torus
  (periodic-boundary) model, independent agents/groups.
- [Simulation_Framework_Cohesion.m](Simulation_Framework_Cohesion.m):
  hard-edge model with a group-cohesion term pulling each agent toward its
  group's centroid.
- [LICENSE](LICENSE): MIT license.

# System Requirements

## Hardware Requirements

A standard desktop or laptop computer is sufficient; no GPU or other
non-standard hardware is required.

## Software Requirements

### OS Requirements

The code is platform-independent (no OS-specific calls). Tested on:

- Windows 11, MATLAB R2023a

### Dependencies

MATLAB only, using exclusively base functions (`rng`, `rand`, `randn`,
`mod`, `mean`, `std`, `cell`, `fullfile`, `mkdir`, `exist`, `csvwrite`,
`squeeze`). **No additional MATLAB toolboxes are required.**

# Installation Guide

Nothing to compile — this is interpreted MATLAB source code.

1. Download or `git clone` this repository.
2. Open MATLAB and add the repository folder to the path.

Typical install time: a few seconds (download only; no build step).

# Demo

Both scripts ship with a small, illustrative parameter set (not the full
manuscript scale) so a first run finishes quickly and confirms the
install works. The random-number generator is seeded
(`rng(42, 'multFibonacci')`) for reproducible output. See
[Instructions for Use](#instructions-for-use) below for the exact
parameters used to produce the manuscript's results.

Before running, open the script (`Simulation_Framework_Main.m` or
`Simulation_Framework_Cohesion.m`) and:

1. **Set `output_folder`** (near the top of `Simulation_Framework()`) to a
   valid, writable path on your system — the shipped value is a
   placeholder and must be changed before running.
2. **Optionally adjust the epidemiological/movement parameters** in the
   same block if you want a different test run.

Then run from the MATLAB command line, with the repository folder on the
path:

```matlab
Simulation_Framework_Main
Simulation_Framework_Cohesion
```

Each takes no arguments.

**Expected runtime** (MATLAB R2023a, Windows 11, standard desktop, shipped
illustrative parameters):

| Script | Runtime |
|---|---|
| `Simulation_Framework_Main.m` | 251.16 s |
| `Simulation_Framework_Cohesion.m` | 207.78 s |

**Expected output:** each run creates a subfolder under the
`output_folder` you set, named after the parameter combination used,
containing six CSV files:

| File | Contents |
|---|---|
| `MeanSirOUT.csv` | Mean S/I/R counts per timestep, averaged across repeats. |
| `StdSirOUT.csv` | Standard deviation across repeats of S/I/R counts per timestep. |
| `Peak.csv` | Per repeat: [peak simultaneously-infected count, timestep of peak]. |
| `Coordinates.csv` | Per timestep, per agent (repeat 1 only): [timestep, agent ID, X, Y, SIR status]. |
| `Raw_Last_Step.csv` | Per repeat: final-timestep [S, I, R] counts. |
| `SIR_Status.csv` | Final SIR status of every agent, in every repeat. |

# Instructions for Use

To run either model with your own parameters:

1. Open `Simulation_Framework_Main.m` or `Simulation_Framework_Cohesion.m`.
2. Edit the parameter block at the top of `Simulation_Framework()` (board
   size, agents/groups, disease and movement parameters, and — for the
   cohesion script — `group_cohesion`).
3. Set `output_folder` to a writable path (e.g. `fullfile(pwd, 'output')`).
4. Run the script. CSVs are written under a subfolder of `output_folder`
   named after the parameter combination (see table above).

To reproduce the manuscript's results (Fig. 3 baseline configuration),
set:

| Variable | Manuscript value | Where set |
|---|---|---|
| `board_size_x`, `board_size_y` | 100, 100 | top of `Simulation_Framework()` |
| `num_agents` | 50 | top of `Simulation_Framework()` |
| `num_infected_agents` × `num_infected_groups` | 2 initial infections total | top of `Simulation_Framework()` |
| `num_groups` | swept: 1, 2, 5, 25 | top of `Simulation_Framework()` |
| `simulation_time` | 2,500+ (run continues until no infected hosts remain; the code itself runs a fixed step count, so set this long enough to reach extinction) | top of `Simulation_Framework()` |
| `burn_in_steps` | 1,000 | **hardcoded inside `simulation()`** — edit that line directly, not the top-level parameter block |
| `num_repeats` | 750 (Fig. 3; other figures use different counts — see the manuscript) | top of `Simulation_Framework()` |
| `pathogen_decay`, `diffusion`, `drift`, `infection_rate`, `recovery_rate` | 0.25, 3, 0.2, 0.1, 0.005 | top of `Simulation_Framework()` — already match the shipped defaults |
| `group_cohesion` (Cohesion script only) | swept: 0, 0.5, 0.75, 1 | top of `Simulation_Framework()` |

For any other figure, use the parameter values given in that figure's
caption in the manuscript.

Two implementation details worth knowing when reusing this code:
- `sir_mode` is a documentation flag; recovery is actually controlled by
  `recovery_rate` (set it to `0` for SI dynamics).


# License

This project is covered under the MIT License — see [LICENSE](LICENSE).

# Citation

If you use this code, please cite:

*"Neglecting movement-disease feedback biases wildlife epidemic
forecasts."*