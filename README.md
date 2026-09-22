# Robust MIMO Control for Microgrid Frequency Stabilisation

The project develops a robust multivariable control system for frequency regulation in an islanded microgrid containing a conventional generator, a wind turbine, and battery storage. The generalized MIMO plant includes multiple disturbance, measurement, control, and performance channels. The controller jointly coordinates the slower conventional generator and the faster battery while penalizing frequency deviation, battery use, and control effort. The model also includes load and wind-power disturbances, measurement noise, and structured uncertainty in the plant dynamics.
A nominal $H_\infty$ MIMO controller is first synthesized and evaluated through structured singular-value analysis. Since it does not satisfy the required robust-stability and robust-performance conditions, the final controller is obtained through $\mu$-synthesis using D-K iteration.

**Author:** Houman Asgari

## Control problem

After linearization around the nominal operating point, the microgrid is represented as an uncertain MIMO system. The two main control inputs act on the conventional generator and the battery, while the main regulated output is the frequency deviation.

The generalized plant can be written in the standard form

$$
\begin{aligned}
\dot{x} &= Ax + B_1 w + B_2 u,\\
z &= C_1 x + D_{11}w + D_{12}u,\\
y &= C_2 x + D_{21}w,
\end{aligned}
$$

where $w$ collects external disturbances and measurement noise, $u$ contains the control inputs, $z$ contains the weighted performance signals, and $y$ is available to the controller.

Model uncertainty is represented by structured multiplicative perturbations collected in a block-diagonal matrix $\Delta$ and connected to the nominal plant through a linear fractional transformation (LFT).

The design objective is to find a stabilizing controller $K$ that

- suppresses frequency deviations caused by load and wind-power fluctuations;
- remains stable for the modeled plant uncertainty;
- satisfies robust-performance requirements;
- coordinates the slow conventional generator with the faster battery storage system;
- limits excessive battery use and control effort.

Frequency-dependent weights are used to shape this behavior. The generator is assigned mainly to slower variations, while the battery is used for faster transients.

## From $H_\infty$ control to $\mu$-synthesis

A nominal $H_\infty$ controller is first obtained by minimizing the worst-case gain from disturbances to weighted performance outputs,

$$
\min_K \left\|T_{wz}(K)\right\|_\infty .
$$

The resulting controller is then evaluated with structured singular-value analysis. For the uncertainty structure $\Delta$, robust performance requires the structured singular value to remain below one over frequency:

$$
\sup_\omega \mu_\Delta\!\left(M(j\omega)\right) < 1 .
$$

The nominal $H_\infty$ design used in the project does not satisfy this condition, motivating a structured robust-control design.

## D-K iteration

Direct minimization of $\mu$ is difficult. D-K iteration instead alternates between two simpler optimization steps.

For a closed-loop interconnection $M$, the structured singular value is upper-bounded by

$$
\mu_\Delta(M)
\le
\inf_D
\bar{\sigma}\!\left(DMD^{-1}\right),
$$

where $D$ is a scaling matrix compatible with the uncertainty structure.

The algorithm proceeds as follows:

### 1. K-step

For a fixed scaling $D$, solve a scaled $H_\infty$ synthesis problem:

$$
\min_K \|D M(K)D^{-1}\|_\infty.
$$

### 2. D-step

With the controller fixed, perform frequency-by-frequency $\mu$-analysis and compute $D(j\omega)$ scalings that tighten the upper bound on $\mu$.

### 3. Fit the scaling

Approximate the frequency-dependent scaling by a stable low-order transfer function so that it can be included in the next synthesis step.

### 4. Repeat

Repeat the K-step and D-step until the robust-performance bound is satisfactory.

In this project, the D-scalings were fitted with third-order transfer functions. Two D-K iterations were sufficient to satisfy the nominal-performance, robust-stability, and robust-performance conditions for the modeled uncertainty set.

## Simulation result

The final controller is tested on a worst-case perturbed plant with load variations, wind-power fluctuations, measurement noise, and model uncertainty.

In the reported simulation:

- the baseline $H_\infty$ controller allows a frequency drop of about **8 Hz**;
- the D-K-iteration controller keeps the peak frequency deviation to about **0.2 Hz**.

The result illustrates why nominal disturbance attenuation alone is not sufficient when structured plant uncertainty is significant.

## MATLAB/Simulink implementation

The repository contains the simulation and controller-design files used in the project. The implementation uses MATLAB/Simulink and the Robust Control Toolbox for

- plant linearization and state-space modelling;
- uncertain LTI and LFT interconnections;
- $H_\infty$ synthesis;
- structured singular-value analysis;
- D-K iteration;
- worst-case perturbation studies.

## Scope and references

This repository is a course-project implementation and numerical study of robust frequency control for an uncertain islanded microgrid.

The project follows the robust-control framework described in:

1. H. Bevrani, M. R. Feizi, and S. Ataee, “Robust Frequency Control in an Islanded Microgrid: $H_\infty$ and $\mu$-Synthesis Approaches,” *IEEE Transactions on Smart Grid*, 2016.
2. Y. Han, P. M. Young, A. Jain, and D. Zimmerle, “Robust Control for Microgrid Frequency Deviation Reduction With Attached Storage System,” *IEEE Transactions on Smart Grid*, 2014.
3. J. C. Doyle, “Analysis of Feedback Systems with Structured Uncertainty,” 1982.

$$
u =
\begin{bmatrix}
u_g\\
u_{\mathrm{batt}}
\end{bmatrix}.
$$

The main regulated output is the frequency deviation $\Delta f$. Battery state of charge and control effort are also penalized in the generalized plant.

The linearized open-loop model is obtained from the Simulink model using `linmod`. The resulting state-space realization

$$
\dot{x}=Ax+Bu,\qquad
y=Cx+Du
$$

has order 12.

## Uncertainty and disturbances

Two plant components are treated as uncertain:

- rotating-mass/load dynamics: $50\%$ uncertainty;
- diesel-engine dynamics: $40\%$ uncertainty.

The uncertainties are represented as multiplicative SISO perturbations and collected into a block-diagonal structured uncertainty matrix

$$
\Delta = \mathrm{diag}(\Delta_1,\Delta_2,\ldots).
$$

The uncertain plant is connected to $\Delta$ through a linear fractional transformation (LFT). This produces the generalized interconnection used for structured singular-value analysis and $\mu$-synthesis.

The disturbance channels include

$$
w =
\begin{bmatrix}
\Delta P_{\mathrm{load}}\\
\Delta P_{\mathrm{wind}}\\
n
\end{bmatrix},
$$

where $n$ denotes measurement noise.

The design therefore addresses two different issues:

1. **disturbance rejection**, mainly against load and renewable-generation fluctuations;
2. **robustness to structured model uncertainty**, mainly in the diesel-generator and rotating-mass dynamics.

## Frequency-dependent performance shaping

The generalized plant uses frequency-dependent weighting functions to encode the desired division of control effort between the conventional generator and the battery.

The generator-control penalty is

$$
W_{cg}(s) = \frac{0.2s+0.1}{100s+0.1}
$$

and the battery-control penalty is

$$
W_{cb}(s) = \frac{s+10^{-4}}{5s+1}
$$

These weights reflect the different actuator bandwidths. The conventional generator is intended to supply slower and larger power variations. The battery is intended to suppress faster transients.

Additional weights reported in the project are


$$
W_{w2}(s)
=\frac{s+5\times 10^{-4}}{s+10^{-5}},
$$


$$
W_{be}(s)
=\frac{20s+100}{s+0.001},
$$

and


$$
W_{se}(s)
=\frac{50s+0.001}{0.5s+0.1}.
$$

Together, the weighting functions penalize frequency deviation, battery usage, and excessive control action over the frequency ranges relevant to each signal.

## Baseline $H_\infty$ design

For a generalized plant $P$ and controller $K$, let

$$
T_{wz}(s;K)
$$

denote the closed-loop transfer matrix from exogenous inputs $w$ to weighted performance outputs $z$.

The nominal $H_\infty$ synthesis problem can be written as

$$
\min_{K\ \mathrm{stabilizing}}
\left\|T_{wz}(s;K)\right\|_\infty.
$$

The resulting controller is designed for the weighted nominal plant. It is then tested against the structured uncertainty using $\mu$-analysis.

Let $\mathcal{U}$ denote the set of admissible structured perturbations. The structured singular value is

$$
\mu_{\mathcal{U}}(M)
=\frac{1}{
\displaystyle
\min_{\Delta\in\mathcal{U}:\,\det(I-M\Delta)=0}
\bar{\sigma}(\Delta)
}.
$$

If no admissible perturbation makes $I-M\Delta$ singular, then
$\mu_{\mathcal{U}}(M)=0$.


A standard robust-performance test is

$$
\sup_{\omega}
\mu_{\boldsymbol{\Delta}}
\left(M(j\omega)\right)<1.
$$

The  $H_\infty$ controller in this project does not satisfy the required $\mu$-based robustness condition. The analysis therefore indicates that an admissible perturbation can violate robust stability or robust performance.

## μ-Synthesis and D-K Iteration

The robust-controller design seeks a controller that minimizes the worst-case structured singular value,

$$
\min_K
\sup_{\omega}
\mu_{\boldsymbol{\Delta}}
\left(M(K,j\omega)\right),
$$

where $M(K,s)$ is the closed-loop interconnection seen by the structured uncertainty blocks.

Direct optimization of $\mu$ with respect to $K$ is difficult. D-K iteration replaces it with alternating controller synthesis and scaling steps. The key upper bound is

$$
\mu_{\boldsymbol{\Delta}}(M)
\leq
\inf_{D\in\mathcal{D}}
\bar{\sigma}
\left(
DMD^{-1}
\right),
$$

where $D$ belongs to a set of scaling matrices that commute with the uncertainty structure.

### K-step

For a fixed scaling $D^{(k)}$, solve a scaled $H_\infty$ problem:


$$
K^{(k+1)}
=\arg\min_K
\left\|
D^{(k)}
M(K)
\left(D^{(k)}\right)^{-1}
\right\|_\infty.
$$

This step finds a controller for the current approximation of the structured robust-performance objective.

### D-step

With $K^{(k+1)}$ fixed, perform frequency-by-frequency $\mu$-analysis and compute scaling matrices that reduce the upper bound,

$$
D^{(k+1)}(j\omega)
\approx
\arg\min_{D\in\mathcal{D}}
\bar{\sigma}
\left[
D
M(K^{(k+1)},j\omega)
D^{-1}
\right].
$$

The resulting $D(j\omega)$ is frequency dependent. To use it in the next synthesis step, the project fits the scaling response with a third-order rational transfer function.

The procedure is then repeated:

1. start with the current controller or scaling;
2. perform $\mu$-analysis over frequency;
3. extract the $D$-scaling matrices;
4. fit the frequency-dependent scaling with a third-order transfer function;
5. solve the scaled $H_\infty$ problem for a new controller;
6. evaluate robust stability and robust performance;
7. repeat until the $\mu$ bound is acceptable.

In the reported implementation, two D-K iterations were sufficient to satisfy the nominal-performance, robust-stability, and robust-performance conditions for the modeled uncertainty set.

The final robust controller has order 12, equal to the order of the design interconnection.

## Interpretation of the MIMO design

The controller acts on the diesel generator and battery simultaneously.

At low frequencies, the conventional generator is preferred because it supplies the bulk of the power and handles slow load variations. At higher frequencies, the battery is preferred because of its faster response.

The frequency-dependent penalties therefore encourage a closed-loop allocation in which

$$
\text{slow power imbalance}
\longrightarrow
\text{conventional generator},
$$

while

$$
\text{fast power imbalance}
\longrightarrow
\text{battery storage}.
$$

This allocation is not imposed by switching logic. It emerges from the weighted MIMO synthesis.

## Worst-case perturbation study

The final $\mu$-synthesis controller is compared with the baseline $H_\infty$ controller under a worst-case perturbed plant.

The test includes:

- the modeled plant uncertainties;
- large load-power variations;
- wind-power fluctuations;
- Gaussian measurement noise.

In the reported simulation, the baseline $H_\infty$ controller allows the system frequency to drop by approximately

$
8\,\mathrm{Hz},
$

under the worst-case perturbation.

With the D-K-iteration controller, the peak frequency variation remains around

$
0.2\,\mathrm{Hz}.
$

The $\mu$-analysis also verifies robust stability and robust performance for the uncertainty model used in the project.

## MATLAB/Simulink implementation

The repository contains the simulation and controller-design files used for the project. The implementation uses MATLAB/Simulink and the Robust Control Toolbox for tasks including

- Simulink plant modelling;
- linearization with `linmod`;
- state-space modelling;
- uncertain LTI models;
- LFT interconnections;
- $H_\infty$ synthesis;
- structured singular-value analysis;
- D-scaling;
- D-K iteration;
- worst-case perturbation simulation.

The original course report does not specify a canonical execution order or document the individual simulation filenames, so no file-by-file run sequence is imposed here.

## Scope

This repository is a course-project implementation and study of robust frequency control for an uncertain islanded microgrid. The contribution is the modelling, MATLAB/Simulink implementation, weighting design, controller synthesis, robustness analysis, and worst-case simulation study.

The robust-control framework follows the literature cited in the original report, in particular:

1. H. Bevrani, M. R. Feizi, and S. Ataee, *Robust Frequency Control in an Islanded Microgrid: $H_\infty$ and $\mu$-Synthesis Approaches*, IEEE Transactions on Smart Grid, 2016.
2. Y. Han, P. M. Young, A. Jain, and D. Zimmerle, *Robust Control for Microgrid Frequency Deviation Reduction With Attached Storage System*, IEEE Transactions on Smart Grid, 2014.
3. J. Doyle, *Analysis of Feedback Systems with Structured Uncertainty*, 1982.
4. S. Skogestad and I. Postlethwaite, *Multivariable Feedback Control: Analysis and Design*.

## Author

**Houman Asgari**  
Sharif University of Technology
