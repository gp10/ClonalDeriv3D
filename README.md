# ClonalDeriv3D :: Quantitative clonal behavior in in vitro 3D culture
Computational methods to analyze clone dynamics in primary organotypic culture of esophageal epithelium explants.

This repository contains code used in the following manuscript:
  > Herms A, Fernandez-Antoran D, Alcolea MP, et al. (2022) Epithelioids: Ultra long-term self-maintaining epithelial tissue equivalents reveal drivers of clonal expansion in normal esophageal epithelium. _submitted_

### Overview
Code in this repository is structured around two different purposes:
- Statistical fit of the experimental clonal data under the assumption of a single progenitor (SP) model paradigm.
- Lattice-based simulations of in vitro cell competition given realistic SP-model parameter values.

To fit experimental clonal data, we rely on statistical theory on neutral competition dynamics, which stablishes that the average clone size should increase linearly with time, while the clone density should decline hyperbolically with time according to the SP-model (Piedrafita et al, 2020). We implement a least-squares minimization algorithm that attempts to fit both types of data all-at-once with a given set of parameter values. Optimum theoretical parameter fits are displayed overlaid on the experimental data trends for the average clone size, clone density and total labelled area.

Given optimum parameter values, we can run realistic simulations of clone dynamics through a 2D lattice-based implementation of the SP-model, a platform first introduced in Colom et al (2020). Simulation results can be overlaid on experimental data too.

### Main scripts
- **Analysis-OptimumFit.m** : implements a least square minimization algorithm to find optimum SP-model parameters fitting both the experimental average clone size and clone density values at different time points. Obtained goodness-of-fit (GoF) values are summarized as a heatmap for the different parameter values tested in the pre-defined grid search. (Complete search takes ~20s on a standard 16GB RAM computer).
  + Optimum SP-model parameter values are saved in `paramFit_optim.mat`) file.
- **Analysis-AvgCloneSize-plot.m** : represents the experimental average clone size over time and allows to overlay the results from the optimum theoretical SP-model fit by loading `paramFit_optim.mat`. Alternatively, it provides the independent, best fit to the average clone size, saved in `paramFit_avg.mat` file.
- **Analysis-ClonalPersis-plot.m** : represents the experimental number of surviving clones over time and allows to overlay the results from the optimum theoretical SP-model fit by loading `paramFit_optim.mat`. Alternatively, it provides the independent, best fit to the average clone size, saved in `paramFit_decay.mat` file.
- **Analysis-TotalLabelledArea-plot.m** : represents the total labelled area (or sum of all clone sizes) over time and allows to overlay the results from the optimum theoretical SP-model fit by loading `paramFit_optim.mat`.
- **Analysis-Simul-2Dgrid.m** : main script to run 2D simulations of clone competition in a 2D lattice and plot the results overlaid on those of the experimental data (average clone size, clone density and total labelled area individual scripts are called). Simulation outcome is shown accompanied by plausible intervals computed from permutation-made subsets that represent limited sample sizes equivalent to those found in the experimental data.
  + Simulation outcome is stored in the `Simulation_data` folder as a `.mat` file.
- **plot_Simul_2Dmovie.m** : represents 2D views of the simulated basal, progenitor cell lattice highlighting simulated clone areas at different time points, and allows to save a video of the simulated clone dynamics.

### Dependencies
- **ClonalData-invitro.mat** : Contains the experimental clonal dataset.
- **Simul-2Dgrid-SPdynamics.m** : Simulator of SP-model clone competition dynamics in a 2D lattice grid. Used in all scripts for model simulations.
- freezeColors.m : Required for color layout in plot_Simul_2Dmovie.m.
- `Simulation_data` folder : contains results from 2D lattice simulations of clone dynamics.

### Requirements
- Matlab (built & tested on version R2020b)
- No dependencies needed
- Runs on all OS
