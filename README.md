# Robust MIMO Control for Microgrid Frequency Stabilisation

MATLAB/Simulink implementation of a robust multivariable controller for frequency regulation in an islanded microgrid with battery energy storage.

This repository contains the simulation files developed for a Robust Control course project at Sharif University of Technology.

## Overview

Islanded microgrids are sensitive to load variations and fluctuations in renewable generation because they do not benefit from the frequency support of a large utility grid. In this project, the microgrid combines:

- a conventional diesel generator,
- a wind turbine generator,
- a battery-based storage system,
- uncertain plant dynamics,
- load and renewable-generation disturbances,
- measurement noise.

The control objective is to reduce frequency deviations while coordinating the slower conventional generator with the faster battery storage system.

## Control problem

The plant is represented by a linearized MIMO model around the nominal operating point. Wind-power fluctuations and load variations are treated as disturbance inputs, while the generator and battery power commands are the control inputs.

Model uncertainty is included in the diesel-generator and rotating-mass dynamics. The uncertain plant is written in a structured form suitable for robust-control analysis.

The main performance objectives are:

- maintain small frequency deviations,
- preserve stability under model uncertainty,
- attenuate load and wind-power disturbances,
- limit excessive battery usage,
- avoid unnecessarily large control effort.

## Method

The controller design follows a structured robust-control approach.

### 1. Nominal and uncertain plant model

The microgrid model is linearized and represented in state-space form. Dynamic uncertainties are introduced through weighting functions and structured perturbation blocks.

### 2. Performance weighting

Frequency deviation, battery state of charge, and control effort are weighted to encode the desired closed-loop behaviour.

The weighting strategy reflects the different roles of the actuators:

- the conventional generator mainly handles slower and larger load variations,
- the battery responds to faster transients,
- frequency deviation is strongly penalized over the relevant bandwidth.

### 3. Baseline \(H_\infty\) design

An \(H_\infty\) controller is first designed for the nominal weighted plant. Structured singular-value analysis is then used to assess robust stability and robust performance.

The baseline controller does not satisfy the required robustness conditions for the full uncertainty set.

### 4. \(\mu\)-synthesis with D-K iteration

A robust MIMO controller is then synthesized using D-K iteration. The procedure alternates between:

1. \(\mu\)-analysis of the closed-loop system,
2. frequency-dependent D-scaling,
3. controller redesign using the scaled plant.

After two D-K iterations, the resulting controller satisfies the robust-stability and robust-performance requirements for the modeled uncertainty set.

## Simulation study

The robust controller is evaluated under a worst-case perturbed plant with:

- strong load-power variations,
- wind-generation fluctuations,
- measurement noise,
- simultaneous model uncertainty.

The report compares the D-K-iteration controller with the baseline \(H_\infty\) controller.

In the considered worst-case simulation, the \(H_\infty\) design exhibits a frequency drop of about **8 Hz**, whereas the robust \(\mu\)-synthesis controller keeps the peak frequency deviation to roughly **0.2 Hz**.

These simulations illustrate the advantage of explicitly accounting for structured model uncertainty during controller synthesis.

## Repository contents

The repository currently contains the MATLAB/Simulink files used for the course-project simulations, including the plant model, uncertainty representation, controller synthesis, and closed-loop robustness tests.

The original project report describes the modelling assumptions, weighting functions, controller design procedure, and simulation results in more detail.

## Software

The project was developed in MATLAB/Simulink and uses functionality from the Robust Control Toolbox.

Relevant methods include:

- state-space modelling,
- linearization of Simulink models,
- \(H_\infty\) control,
- structured singular-value (\(\mu\)) analysis,
- linear fractional transformations,
- D-K iteration,
- robust stability and robust performance analysis.

## Background

The project follows the robust-control framework commonly used for uncertain multivariable systems and is motivated by frequency regulation in islanded microgrids with renewable generation and battery storage.

Two key references used in the original report are:

1. H. Bevrani, M. R. Feizi, and S. Ataee, “Robust Frequency Control in an Islanded Microgrid: \(H_\infty\) and \(\mu\)-Synthesis Approaches,” *IEEE Transactions on Smart Grid*, 2016.
2. Y. Han, P. M. Young, A. Jain, and D. Zimmerle, “Robust Control for Microgrid Frequency Deviation Reduction With Attached Storage System,” *IEEE Transactions on Smart Grid*, 2014.

## Author

**Houman Asgari**  
Sharif University of Technology
