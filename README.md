# Stochasticity in multi-layered Quorum Sensing systems

This repository contains MATLAB scripts and functions for the analysis of the paper:
"Stochasticity in multi-layered Quorum Sensing systems".

It describes _Burkholderia thailandensis_ quorum sensing network, that has three  LuxI/LuxR-type circuits: BtaI1/BtaR1 (QS-1), BtaI2/BtaR2 (QS-2), and BtaI3/BtaR3 (QS-3). They synthesise and respond to distinct signalling molecules: HC8, HC10 and C8, and are hierarchically organised and interconnected, resulting in sequential activation and coordinated regulation of gene expression.

The workflow compares a **Complete Model** (QS-1, QS-2 and QS-3) against a **Reduced Model** (QS-3 removed) in terms of the second moments, to assess the stochastic fluctuations simulated via a linearised Langevin equation.

The main steps of this analysis are:
- Assessing the role of kinetic parameters in shaping fluctuations: parameter scan over the inhibition and activation rates $k_{in}$ and $k_{att}$
- Validate analytical predictions using the fitted data and the set of parameters found.

---

## How to cite

paper citation...

---

## File Overview

### Data

| File | Description |
|------|-------------|
| `experimentalData.dat` | Experimental time-series data for HC8, HC10, and C8 (from https://doi.org/10.1128/jb.00727-17). |
| `parLong.dat` | Fitted model parameters. |

---

### Functions

#### `sistemODE.m`
Defines the deterministic ODE system for the **Complete Model** (11 species).

`dydt = sistemaODE(~, y, params)`

- `y` — state vector of length 11: [R1, R2, R3, M, I, R1\*, R2\*, R3\*, HC8, HC10, C8]
- `params` — parameter vector of length 10 (see parameter table below)
- Returns `dydt`, the time derivative of the state vector.

---

#### `sistemODE_12.m`
Defines the deterministic ODE system for the **Reduced Model** (9 species, R3 and R3\* absent).

`dydt = sistemaODE_12(~, y, params)`

- `y` — state vector of length 9: [R1, R2, M, I, R1\*, R2\*, HC8, HC10, C8]
- `params` — parameter vector of length 10
- Returns `dydt` of length 9.

---

#### `calcMat_M_B.m`
Computes the matrix **M** and **B** for the **Complete Model** around the equilibrium point.

`[M, B] = calcMat_M_B(params, Eq)`

- `params` — parameter vector of length 10
- `Eq` — equilibrium concentration vector of length 11
- Returns `M` (11×11) and `B` (11×11, symmetric).

---

#### `calcMat_M_B_12.m`
Same as above for the **Reduced Model**. Builds the full 11×11 matrices internally, then removes the rows and columns corresponding to R3 and R3\*.

`[M, B] = calcMat_M_B_12(params, Eq)`

- `params` — parameter vector of length 10
- `Eq` — equilibrium concentration vector of length 9
- Returns `M` (9×9) and `B` (9×9, symmetric).

---

#### `matSecMoms.m`
Builds the linear system for the steady-state second moments from M and B.

`[bigM, bigB] = matSecMoms(M, B, n)`

- `M` — n×n matrix
- `B` — n×n matrix
- `n` — number of species
- Returns `bigM` (nMoments×nMoments) and `bigB` (nMoments×1), where `nMoments = n*(n+1)/2`
- The steady-state second moments are obtained as `secMom = -inv(bigM) * bigB`

---

#### `solveLangevinLinear.m`
Implements the linearised Langevin equation using the **Euler-Maruyama** algorithm:

```
d(csi) = M*csi*dt + G*dW
```

where G is the Cholesky factor of B (satisfying B = G\*G').

`[t, sol] = solveLangevinLinear(M, G, nSteps, dt, startPoint, n)`

- `M` - (n×n) drift matrix
- `G` - (n×n) noise matrix, Cholesky factor of B
- `nSteps` - total number of time steps
- `dt` - integration time step
- `startPoint` - (n×1) initial condition (typically zeros for fluctuations around equilibrium)
- `n` - system dimension
- Returns `t` (1×nSteps time vector) and `sol` (nSteps×n trajectory matrix)


---

## Analysis and plot Scripts

#### `plot_det_sys_exp_data.m`
Plots the deterministic solution of the Complete Model against experimental data for the three measured species: HC8, HC10, and C8.

**Output:** A 3-panel figure with error bars from `experimentalData.dat` (Figure 2).

**Depends on:** `sistemODE.m`, `experimentalData.dat`, `parLong.dat`

---

#### `scan_KinKatt_CompVSReduc.m`
Performs a 2D parameter scan over a grid of kin and katt values for two alpha3 values, computing standard deviations and convergence plots for both CM and RM.

**What it does:**
1. Loops over `alpha3Vet` (2 values), `k_ins` (11 values), and `k_atts` (18 values).
2. For each combination, solves the ODE system, computes second moments, and
   extracts the standard deviation of each species (sqrt of diagonal second moments).
3. Optionally saves convergence plots as PDF files.

**Output files:** 
- `scan_secMom_a3_{alpha3}_kin_{kin}_katt_{katt}_CM/RM.dat`
- `scan_eqPoints_a3_{alpha3}_kin_{kin}_katt_{katt}_CM/RM.dat`
- convergence plots PDF

**Depends on:** `sistemODE.m`, `sistemODE_12.m`, `calcMat_M_B.m`,
`calcMat_M_B_12.m`, `matSecMoms.m`

---

#### `plot_scan_intensity.m`
Visualises the ratio of CM to RM second moments across the kin–katt parameter grid as scatter heatmaps.

**What it does:**
1. Loads scan results for both alpha3 values and computes the ratio `secmoms_CM / secmoms_RM` for the 9 shared species, for each `k_ins` and `k_atts` combination.
2. Produces a 2×3 figure showing the last 3 species (HC8, HC10, C8) for both
   alpha3 values (one row per alpha3, Figure 3).
3. Optionally produces a 3×3 figure for all 9 species for a single alpha3.

**Depends on:** output of `scan_KinKatt_CompVSRidott.m`

---

#### `langevin_sims.m`
Computes and saves equilibrium points and steady-state second moments for all combinations of parameter type (FIT/MIN) and alpha3 values.
Runs `nAverages` independent Langevin simulations for the same combinations of parameter type and alpha3, for both CM and RM.

**What it does:**
1. Loops over `typeVet = {'FIT', 'MIN'}` and `alpha3Vet = [0.03, 0.3]`
2. For each combination, solves the deterministic system (CM and RM) to find the equilibrium points, then computes the second moments analyticallycomputes M and B, and derives G via Cholesky decomposition
3. Runs `nAverages = 25` independent trajectories using `solveLangevinLinear` and do a subsample every `skip = 1e5` steps before saving

**Key parameters:**
- `dt = 1.25e-4` — integration time step
- `nSteps = 9.6e8` — total steps per trajectory
- `skip = 1e5` — subsampling interval

**Output files:** 
- `{type}_eq_Points_a3_{alpha3}_CM/RM.dat`
- `{type}_secMom_a3_{alpha3}_CM/RM.dat`
- `{type}_{alpha3}/nAv_CM/fileOut_nAV_{i}_CM/RM.dat`

**Depends on:** `sistemODE.m`, `sistemODE_12.m`, `calcMat_M_B.m`,
`calcMat_M_B_12.m`, `solveLangevinLinear.m`

---

#### `plot_histograms.m`
Loads Langevin simulation output and theoretical second moments for all 4 cases (FIT/MIN × alpha3=0.03/0.3), computes normalised mean histograms, and overlays theoretical Gaussian curves from the linearised second moments.

**What it does:**
1. For each case, loads all `nAverages` trajectory files for CM and RM
2. Finds a common bin range across averages, then recomputes all histograms on the same edges and averages them
3. Normalises to area = 1 and plots alongside theoretical Gaussians

**Output:** One 3×3 figure per case (9 subplots, one per species in the RM), with CM and RM bars and theoretical curves overlaid.

**Depends on:** output of `langevin_sims.m`

---

## Suggested Workflow

```
1. plot_det_sys_exp_data.m  →  check deterministic fit to data
2. scan_KinKatt_CompVSRidott.m  →  run kin/katt parameter scan
3. plot_scan_intensity.m    →  visualise scan results
4. langevin_sims.m          →  compute equilibrium points and second moments and run stochastic simulations
5. plot_histograms.m        →  compare simulation histograms to theory

```
---

## Parameter Vector

All scripts share the same 10-element `params` vector:

| Index | Symbol | Description |
|-------|--------|-------------|
| 1 | beta | Basic transcription rate |
| 2 | mu | Molecule degradation rate |
| 3 | kon (primary) | Primary association rate constant |
| 4 | koff (primary) | Primary dissociation rate constant |
| 5 | kon (secondary) | Secondary association rate constant |
| 6 | koff (secondary) | Secondary dissociation rate constant |
| 7 | kin | Inhibition rate constant |
| 8 | katt | Activation rate constant |
| 9 | alpha2 | R2 production rate |
| 10 | alpha3 | R3 production rate |

---

## Species Ordering

| Index | CM species | RM species |
|-------|-----------|-----------|
| 1 | R1 | R1 |
| 2 | R2 | R2 |
| 3 | R3 | *(absent)* |
| 4 | M | M |
| 5 | I | I |
| 6 | R1\* | R1\* |
| 7 | R2\* | R2\* |
| 8 | R3\* | *(absent)* |
| 9 | HC8 | HC8 |
| 10 | HC10 | HC10 |
| 11 | C8 | C8 |

