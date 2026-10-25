# CRN_v8 – ChemPlasKin for reactor networks, with an analytic Jacobian

ChemPlasKin solves coupled gas–plasma kinetics in 0D reactors:
- the chemistry comes from Cantera;
- the electron energy distribution function (EEDF) comes from the CppBOLOS Boltzmann solver;
- the integrator is CVODE.

`CRN_v8` is `CRN_v7` with the option of an **analytic Jacobian** for CVODE (section 7). The Jacobian is built by `ChemPlasReactor::jacobian` and passed to CVODES through a new integrator (`cvodeJac.h`), instead of being computed by finite differences, which costs N right-hand-side evaluations per Jacobian.

On the CRN_serialPlasma case (full run to steady state, W 1000, k 10, tol 1e-3) the run takes 230 s of CPU time against 299 s for v7 (**−23%**). The reactors freeze at the same pulses (R1 ~8200, R2 ~35 100), and the steady states agree within 4·10⁻⁵.

`CRN_v7` is `CRN_v6` with a **steady-state criterion for the network**. A reactor that reaches steady state, together with the reactors upstream of it, is no longer integrated: its outflow stays at the average of its last block of pulses. The run stops when the whole network is at steady state, and a frozen reactor is integrated again if its inflows change.

On the CRN_serialPlasma case, R2 (residence time ~2400 pulses) is not at steady state after 7500 pulses. v7 brings it there: with the `autoplay` settings R1 freezes at pulse ~8200, R2 at pulse ~35 100, in about 300 s of CPU time. At steady state the outlet NO is 35% and H₂O 34% higher than at 7500 pulses.

New in v7:

| Feature | Description |
|---|---|
| Steady-state monitor (`steadyState.h`) | A sample (pulse average) every W/k pulses, compared with the sample W pulses before; relative criterion projected over the residence time; one comparison per group (section 6) |
| Recirculation | The reactors of a loop of the flow graph form one group with a common counter: they converge together |
| Mass balance | A reactor with a net inflow or outflow is never at steady state |
| Freezing and release | A reactor at steady state is no longer integrated. It is integrated again if a flow rate changes, or the state of an inflow from a reservoir or a non-frozen reactor changes |
| Ep not reset in reactors without discharge | An offset replaces the reset: Ep does not enter the other equations, so CVODE is not restarted at every pulse when the inflows do not change. This makes the phase in which a downstream reactor advances alone fast |
| Exact update times | startTime + n·`flowUpdateInterval` (round-off used to build up) |
| Last row at runTime | At the end of the run every reactor, frozen or not, writes a row at the final time |

With `steadyStateW 0` (default) the monitor is off: R1 is bit-identical to v6, and R2 differs by at most 10⁻⁴ because of Ep.

v6 integrates the reactors independently up to the flow update:
- CVODE is restarted only when the flow terms change;
- reactors without discharge take a single integration;
- the outflow averages use the trapezoidal rule on the CVODE internal steps.

v6 is built on v5 (EEDF table), and v5 is v4 with the table: everything else is summarized below.

v4 extends the code to a **chemical reactor network** (CRN). Perfectly stirred reactors (PSR) are connected by mass flow controllers. The flow rates are computed at every update by a Python script of the case. One or more reactors can have a nanosecond repetitively pulsed (NRP) discharge.

Compared with the original version (CRN_v3), v4, and therefore all later versions, adds:

| Feature | Description |
|---|---|
| Boltzmann state per reactor | Each reactor has its own EEDF and rate constants. The EEDF is recomputed only when E/N, T or the mixture change beyond the tolerances |
| CVODE restarted only when needed | Only after a new EEDF, new flows, the start of a pulse or the temperature jump of the heat-loss model. Option `reinitEveryStep` restores the v3 behaviour |
| `flowUpdateInterval` | The flows are updated at fixed intervals. Between two updates each downstream reactor receives the **time average** of the upstream reactor outflow |
| `writeInterval` | One output row every `writeInterval` seconds instead of every step |
| Per-pulse file | `<reactor>_pulses.csv`, one row per pulse |
| Exact end of the discharge | The discharge stops exactly when the deposited energy reaches `pulseEnergy`. The time t* is located on the CVODE interpolant; from t* E/N goes back to the base value, a new EEDF is computed and CVODE is restarted |
| CVODE tolerances applied | `reltol` and `abstol` of the `odes` block are passed to CVODE before its initialization (they used to be ignored) |
| Stop criterion | `dTdt_min < 0` disables the stop on dT/dt |

---

## 1. Folder layout

```
CRN_v8/
├── masterCRN.cpp       main program: case input, network, time loop, freezing of the reactors, choice of the integrator
├── cvodeJac.h          CRN::CVodeJac: CVODES integrator with the analytic Jacobian (or CVODE's finite differences)
├── steadyState.h       steady-state monitor: graph groups (Tarjan), criterion d, counters, release
├── plasmaReactor.h     ChemPlasReactor class: gas-plasma reactor equations (right-hand side and Jacobian for CVODE)
├── crnDict.h           crnDict reader (reactors and controllers), shared by ChemPlasKin and EEDFTable
├── eedfTable.h         EEDF table: point, settings, lookup, interpolation, file input and output
├── eedfTableBuilder.cpp EEDFTable program: builds the table in advance
├── CMakeLists.txt      build (ChemPlasKin and EEDFTable)
├── controlDict         example case (time control)
├── chemPlasProperties  example case (physics and plasma)
├── crnDict             example case (reactors and mass flow controllers)
├── massflowsolver.py   example flow-rate solver
├── massflow_params.txt parameters read by massflowsolver.py
├── mechanisms.txt      notes on the available mechanisms
└── build/              ChemPlasKin executable (after the build)
```

Dependencies, with paths relative to the `CRN_v8` folder:
- `../../../cantera`: the ChemPlasKin fork of Cantera (headers in `include/`, libraries in `build/lib`, `libcantera_shared`), which contains CppBOLOS;
- `../../src/utilities.h`: `readParameter` and the reader of OpenFOAM-style files;
- Python 3 with development files (`Python.h`): the flow-rate solver is called from C++ through the embedded interpreter;
- for `massflowsolver.py`: `numpy` and `scipy`.

## 2. Build

```bash
cd CRN_v8
mkdir -p build && cd build
cmake ..
make            # builds build/ChemPlasKin and build/EEDFTable
```

## 3. Run

```bash
cd CRN_v8/build
./ChemPlasKin -case <case_folder> [-log]
```

- `-case`: folder with the case files. Without `-case`, `..` is used, i.e. `CRN_v8/`.
- `-log`: detailed output at every step and CppBOLOS info messages. Very slow: for debugging only.

During the run a progress line is printed (`runTime`, pulse, number of steps). At the end the program prints the reason for stopping, the CPU time and a summary of the table: `EEDF: N from table, M computed and added to <file>`.

## 4. Case files

### 4.1 Required files

| File | Content |
|---|---|
| `controlDict` | times, steps, CVODE tolerances, update and output intervals |
| `chemPlasProperties` | default conditions, thermal model, plasma and CppBOLOS parameters, data files |
| `crnDict` | reactors, reservoirs and mass flow controllers |
| `massflowsolver.py` | function `mainSolver(ox_Y, fuel_Y)` returning the flow rates |
| mechanism `.yaml` | gas + plasma kinetics, in the extended Cantera format of ChemPlasKin (`type: Boltzmann` reactions). Path in `mechFile` or `mechanismFile` |
| LXCat cross sections `.dat` | data for CppBOLOS. Path in `csDataFile` |

### 4.2 Optional files

| File | Content |
|---|---|
| `massflow_params.txt` | parameters for `massflowsolver.py` (in the test solver: air flow rate [kg/s] on the first line, global equivalence ratio on the second). Needed only if the solver reads it |
| `eedfTable.yaml` | EEDF table (section 5). Created if missing; if present it is read and new points are appended at the end |

All dictionary files use the OpenFOAM-style syntax: `key value;`, blocks `name { ... }`, comments `//` and `/* */`.

### 4.3 `controlDict`

| Key | Req. | Unit | Meaning |
|---|---|---|---|
| `startTime` | yes | s | start time (in `autoplay`, the end time of the Cantera phase) |
| `endTime` | yes | s | end time. The run ends at max(`endTime`, `startTime` + `nPulses`/`f_NRP`), or earlier when the network reaches steady state |
| `nPulses` | yes | – | number of pulses. With the monitor on it is a maximum (in `autoplay`: 1e5) |
| `endTimeForce` | no | s | if given, upper bound of the end time |
| `dt1` | yes | s | fixed step during the discharge of the first two pulses |
| `dt_max` | yes | s | largest step |
| `pulseExten` | yes | s | length of the fine-step phase from the start of each pulse. Clipped between `pulseWidth` and the period |
| `NstepInPulse` | yes | – | number of (logarithmic) steps during the discharge, from the third pulse on |
| `NstepAfPulse` | yes | – | number of (logarithmic) steps from the end of the discharge up to `pulseExten` |
| `dT_max` | yes | K | ignition detection threshold: if dT > `dT_max` and dT/dt > 1000 K/s the step is reduced. At least 1 K |
| `dTdt_min` | yes | K/s | stop when max(dT/dt) ≤ `dTdt_min`. A negative value disables the stop |
| `kOutPressure` | yes | – | read but not used (left over from the Cantera pressure controller) |
| `flowUpdateInterval` | **yes, > 0** | s | flow update interval; also the interval over which the reactors are integrated independently. Each downstream reactor receives the source state averaged over the previous interval. Recommended: `1/f_NRP`. With 0 the program stops with an error |
| `writeInterval` | no | s | interval between output rows. `0` (default) = every step |
| `reinitEveryStep` | no | 0/1 | `1` = restart CVODE at every step, as in v3. Default `0` |
| `integrateToFlowUpdate` | no | 0/1 | `1` (default) = reactors without discharge (`discharge no`) take one integration up to the flow update; `0` = they take the `findTimeStep` steps as the discharge reactors. Always 0 with `reinitEveryStep 1` |
| `steadyStateW` | no | pulses | W: window of the comparison (each sample is compared with the sample W pulses before). `0` (default) = monitor off |
| `steadyStateSamples` | no | – | k: samples per window, i.e. one sample (and one comparison) every W/k pulses. Default 10 |
| `steadyStateCluster` | no | pulses | interval between two samples, instead of W/k |
| `steadyStateConfirm` | no | – | consecutive comparisons below the threshold to declare steady state. Default 1 (a comparison already spans W pulses) |
| `steadyStateTol` | no | – | threshold on d. Default 1e-3 |
| `steadyStateYmin` | no | – | mass fractions below this value (in both samples) are not checked. Default 1e-12 |
| `steadyStateTau` | no | s | time scale τ of d for all reactors. `0` (default) = τ of each reactor (section 6.2) |
| `steadyStateStop` | no | 0/1 | `1` (default) = stop the run when the whole network is at steady state |
| `jacobian` | no | `cantera` / `fd` / `analytic` | integrator of the reactors (section 7). `cantera` (default) = Cantera's CVodesIntegrator, as in v7; `fd` = CVodeJac with the finite-difference Jacobian of CVODE (same steps as `cantera`); `analytic` = CVodeJac with `ChemPlasReactor::jacobian` |
| `jacobianCheck` | no | N | `N > 0`: every N-th analytic Jacobian is compared with a finite-difference one, and the largest scaled difference is printed (diagnostics, slow). Default 0 |
| `jacobianSkipFalloff`, `jacobianSkipThirdBodies` | no | 0/1 | leave the falloff or third-body terms out of the kinetics derivatives (approximate Jacobian). No gain on CRN_serialPlasma. Default 0 |

With the monitor on, `flowUpdateInterval` and the period 1/`f_NRP` must be multiples of each other, because the checks are taken on pulse boundaries.

`autoplay` values: `steadyStateW 1000; steadyStateSamples 10; steadyStateTol 1e-3; steadyStateYmin 1e-10; jacobian analytic;`.
| `odes { reltol; abstol; }` | no | – | CVODE tolerances. Defaults 1e-9 and 1e-15. `abstol` applies to mass fractions, and electrons have Y ≈ 1e-15: with `abstol` > 1e-13 the electron density is no longer controlled |

Example:
```
startTime   3.14891; endTime 3.64891; nPulses 7500;
dt1 1e-11;  dt_max 1e-5;
pulseExten 1e-5; NstepInPulse 400; NstepAfPulse 200;
dT_max 10.0; dTdt_min -1e-6; kOutPressure 0.0;
flowUpdateInterval 6.666666666666667e-05;   // 1/f_NRP
writeInterval 1e-6;
odes { reltol 1e-9; abstol 1e-15; }
```

### 4.4 `chemPlasProperties`

| Key | Req. | Unit / values | Meaning |
|---|---|---|---|
| `Temperature` | yes | K | default T of reactors without `initialTemperature` |
| `Pressure` | yes | Pa | default p of reactors without `pressure` |
| `CP_or_CV` | yes | `CP` / `CV` | constant-pressure or constant-volume reactor |
| `nonThermal` | yes | yes/no | `yes` = T, electron energy and N2 vibrational energy frozen: chemistry only |
| `plasmaHeatModel` | no | `Flitti` / other | gas heating by the plasma. `Flitti` = global model of Flitti & Pancheshnyi 2009; otherwise the detailed model, reaction by reaction |
| `heatLoss` | yes | yes/no | radial heat loss q = C0·λ·(T − T0). When on, the reactor starts at constant volume, and 100 ns after each pulse an isentropic expansion brings it back to the reactor pressure |
| `C0` | if `heatLoss` | 1/m² | heat loss coefficient (in `autoplay`: 1/r²) |
| `Mixture { species x; }` | yes | mole fractions | default composition; also used to initialize the Boltzmann species |
| `fuelName`, `oxName` | yes | element | elements (e.g. `H`, `O`) whose elemental mass fraction is computed in every reactor and reservoir and passed to `mainSolver` |
| `electronDens` | yes | #/cm³ | default initial electron density, if `Electron` is not in the composition |
| `inertSpecies` | yes | species | species that closes the sum of the mass fractions |
| `BoltzmannSpecies` | yes | list | target species of CppBOLOS (they must exist in the cross-section file) |
| `gridType` | yes | e.g. `QuadraticGrid` | energy grid of CppBOLOS |
| `gridSize { start; end; points; }` | yes | eV, eV, – | grid bounds and number of points |
| `E/N(NRP)` | yes | Td | reduced field during the discharge |
| `E/N(DC)` | yes | Td | base field between pulses (max(E/N(DC), 0.1 Td) is used) |
| `stayDC` | yes | yes/no | after the last pulse: `yes` keeps E/N(DC), `no` switches to 0.1 Td |
| `f_NRP` | yes | Hz | pulse frequency |
| `pulseEnergy` | yes | mJ/cm³ | energy per pulse: the discharge stops when it is reached |
| `pulseWidth` | yes | s | largest discharge duration |
| `TGAS_TOLERANCE` | yes | K | change of T beyond which the EEDF is recomputed |
| `EN_TOLERANCE` | yes | Td | change of E/N beyond which the EEDF is recomputed |
| `csDataFile` | yes | path | LXCat cross sections |
| `mechFile` | yes | path | default mechanism |

The EEDF is also recomputed when the mole fraction of a Boltzmann species changes by more than max(0.01, 10%) (internal criterion, `plasmaReactor.h`).

Keys added by v5, all optional:

| Key | Default | Meaning |
|---|---|---|
| `eedfTable` | `<case>/eedfTable.yaml` | table file. `<case>` is replaced by the case folder; `none` = table in memory only, not saved |
| `eedfTableMaxGap` | 5 | [K] largest T distance between two points for interpolation. `0` = exact T only; `< 0` = no reuse, every EEDF is computed as in v4 |
| `eedfTableXrel` | 0 | relative tolerance on the mixture: a point is valid if \|x − x_point\| ≤ max(Xabs, Xrel·x_point) for every Boltzmann species |
| `eedfTableXabs` | 0 | absolute tolerance on the mixture (0 and 0 = identical mixture only) |
| `eedfRtol` | 1e-5 | convergence tolerance of CppBOLOS (`converge()`) |
| `eedfMaxIterations` | 200 | largest number of CppBOLOS iterations |

`autoplay` uses `eedfTableMaxGap 5`, `eedfTableXrel 0.01`, `eedfTableXabs 0.001`.

### 4.5 `crnDict`

```
reactors
{
    R1
    {
        type                PSR;          // PSR (integrated) or reservoir (imposed state)
        volume              0.00025;      // m^3, required for PSRs
        discharge           yes;          // optional (default no): NRP pulses in this reactor
        initialTemperature  800;          // K   (default: Temperature)
        pressure            101325;       // Pa  (default: Pressure)
        initialElectronDensity 1e10;      // #/cm^3 (default: electronDens)
        composition { N2 0.712; O2 0.200; H2 0.088; }   // mole fractions (default: Mixture)
        mechanismFile       "mech.yaml";  // optional (default: mechFile)
        outputFile          "out/R1.csv"; // optional (default: <case>/output/R1.csv)
    }
    res_air { type reservoir; initialTemperature 800; pressure 101325; composition { N2 0.781; O2 0.219; } }
    outlet  { type reservoir; ... }
}
massFlowControllers
{
    airI { start res_air; end R1; value 0.0; }   // value: ignored, overwritten by mainSolver
    mA   { start R1;      end R2; value 0.0; }
}
```

- **Reservoirs.** Their state is fixed. A reservoir that receives flows (e.g. `outlet`) takes the mixture of its inflows at every update, at its own pressure.
- **Electrons.** If the composition does not contain `Electron`, the initial electron density is `initialElectronDensity` (or `electronDens`).

### 4.6 `massflowsolver.py`

The script is imported from the case folder. It must define:

```python
def mainSolver(ox_Y, fuel_Y):
    # ox_Y, fuel_Y: {node name: elemental mass fraction of oxName / fuelName} for every reactor and reservoir
    return {"airI": 2.744e-3, "fuelI": 1.85e-5, "mA": 2.76e-3, "out": 2.76e-3}   # kg/s, one key per controller of crnDict
```

The function is called at every flow update. If a controller is missing, the run stops with an error. The example file contains several variants (`mainSolver2`, `mainSolverFOLC`, `mainSolverMAD`); the one called is always `mainSolver`.

## 5. EEDF table

### 5.1 Content

A point of the table is a converged CppBOLOS solution for one set (E/N, T, mixture of the Boltzmann species). It holds what the reactors need:
- the electron temperature Te;
- the power absorbed per particle;
- the rate constant of every Boltzmann process.

The file is YAML: a header with the settings, then one point per line. New points are appended without rewriting the file.

```yaml
format: CRN-EEDF-table-1
settings:
  csDataFile: "bolsigdb_air_NH3_H2.dat"
  csDataFileSize: 753880
  grid: {type: QuadraticGrid, start: 0, end: 60, points: 160}
  convergence: {maxIterations: 200, rtol: 1e-05}
points:
  - {EN: 948.8, T: 800.0, Te: 123860.8, power: 5.93e-13, X: {"H2": 0.088, "N2": 0.712, ...}, k: {"N2 -> N2(B3)": 1.22e12, ...}}
```

- **Units:** EN in Td, T and Te in K, power in eV·m³/s, X mole fractions, k in m³/kmol/s. Numbers are written with 17 digits, so a mixture read back is identical to the one stored.
- **Different settings:** if the grid, the cross sections (file name and size) or the convergence of the case differ from those of the table, the table is rejected with an error. Do not edit the settings by hand.

### 5.2 Use in ChemPlasKin

When a reactor needs a new EEDF (see the criteria in section 6), it is looked up in the table through an index on E/N and mixture, with the points sorted by T:
1. **Interpolation:** if two points with the same E/N (within `EN_TOLERANCE`) and the same mixture (within Xrel/Xabs) bracket T and are at most `eedfTableMaxGap` apart, the EEDF is interpolated linearly in T.
2. **Nearest point:** otherwise the closest point within `TGAS_TOLERANCE` is used.
3. **Computation:** otherwise the Boltzmann equation is solved, starting from the previous EEDF of the reactor, and the point is added to the table and to the file.

Reusing the table in a second run changes the results by less than ~3·10⁻⁴.

### 5.3 Building the table in advance: `EEDFTable`

```bash
cd CRN_v8/build
./EEDFTable -case <case> -T <list> [-EN <list>] [-comp case|table|<file>]... [-table <file>] [-tolT <K>] [-tolEN <Td>]
```

| Option | Meaning |
|---|---|
| `-case` | case folder: settings from `chemPlasProperties`, mixtures from `crnDict` |
| `-T` | temperatures [K]: values and ranges `start:end:step` separated by commas, e.g. `800:2500:1` |
| `-EN` | E/N [Td], same syntax. Default: the values used by the case, i.e. E/N(NRP), max(E/N(DC), 0.1) and 0.1 |
| `-comp` | mixtures: `case` = compositions of the reactors in `crnDict` (default); `table` = mixtures already in the table; `<file>` = one mixture per line, e.g. `H2:0.07, O2:0.20, N2:0.73`. Can be repeated |
| `-table` | table file. Default: `eedfTable` of `chemPlasProperties`, or `<case>/eedfTable.yaml` |
| `-tolT`, `-tolEN` | distance below which a point counts as already present. Defaults: `TGAS_TOLERANCE`, `EN_TOLERANCE` |

The program computes every combination (E/N, T, mixture) not yet present, sweeping T in increasing order so that each EEDF starts from the previous one. The functions are `addPoint` (one point) and `buildTable` (the grid), in `eedfTableBuilder.cpp`.

## 6. Steady state of the network

### 6.1 Groups

At start the flow graph is built (one node per reactor and reservoir, one edge per controller), and its strongly connected components are computed with Tarjan's algorithm. A **group** is a single reactor or the set of reactors of a recirculation loop. Groups are ordered from upstream to downstream: for CRN_serialPlasma `[res_air] [res_H2] [R1] [R2] [outlet]`.

Reservoirs are at steady state by definition and are not checked.

### 6.2 Criterion

Every w = W/k pulses (k = `steadyStateSamples`, at a flow update on a pulse boundary) each reactor not frozen takes a **sample**:
- T, p and the mass fractions Y_k, **averaged over the w pulses** since the previous sample: this is the state passed downstream. The instantaneous state of the short-lived species at the end of a period jumps between pulses without any drift;
- ṁ_in and ṁ_out, constant between flow updates, at the time of the sample.

Each sample is compared with the sample of **W pulses before** (k samples back). The last k + 1 samples are kept, and the oldest one is dropped once it has been used. So a comparison over W pulses is made every W/k pulses, not every W. With W = 100 and k = 4: samples at 25, 50, 75, 100, 125, …, comparisons 25–125, 50–150, 75–175, … For each quantity x:

    d = |x(n) − x(n−k)| / max(|x(n)|, |x(n−k)|) · τ / (t(n) − t(n−k))

A real drift gives the same d whether two adjacent samples or two samples W pulses apart are compared. The noise between samples, instead, is divided by the whole window W, so it weighs k times less than between adjacent samples. Until W pulses of samples are available (at the start, or after a release) no comparison is made.

and for the mass balance:

    d_m = |ṁ_in − ṁ_out| · τ / m,   m = ρ·V

d is the relative change the quantity would have over τ if the trend of the last block went on: an estimate of the distance from steady state. One criterion covers both needs:
- the change is relative, so a small change of a small quantity stays large;
- the derivative is taken over a whole block of pulses, never over a small time step.

Mass fractions below `steadyStateYmin` are not checked, because at that level the check would measure noise.

**Time scale τ** (when `steadyStateTau` is not given):

| Case | τ | d_m |
|---|---|---|
| ṁ_out > 0 | m/ṁ_out (residence time) | \|ṁ_in − ṁ_out\|·τ/m |
| ṁ_out = 0, ṁ_in > 0 | m/ṁ_in (filling time) | 1: never at steady state |
| ṁ_out > 0, ṁ_in = 0 | m/ṁ_out | 1: never at steady state |
| closed reactor | W·period (the window) | 0 |

### 6.3 Counter and freezing

At every comparison, for every group not yet at steady state:
- if d < `steadyStateTol` for every quantity of **every** reactor of the group, and every node upstream of the group is at steady state, the counter n grows by 1;
- otherwise n goes back to 0: a single reactor of the group above the threshold resets it.

The group is at steady state when n ≥ `steadyStateConfirm` (default 1). A single comparison already spans W pulses, so steady state is detected at most W/k pulses after the whole window has settled. Its reactors are no longer integrated, and their outflows keep the average of the last sample. The downstream groups go on; with fixed inflows CVODE is not restarted, so they advance very fast. When all groups are at steady state (and `steadyStateStop 1`) the run stops.

### 6.4 Release

At every flow update a frozen reactor is integrated again, together with its whole group, if:
- any flow rate changes, inflow or outflow;
- or the composition, T or h of an inflow from a reservoir or from a non-frozen reactor changes.

Changes of state that come from a frozen reactor do not count: its last average reaches the downstream reactors once, right after it freezes. This also holds for the reactors of a recirculation loop that freeze together.

When a group is released, the counter and the samples start again. The pulse numbering is moved to the current pulse: if the release falls on a pulse boundary, the discharge starts at once. Ep is reset. If an upstream reactor is released, its outflow changes and the downstream reactors are released at the next update.

### 6.5 Choice of the parameters

The noise between samples is amplified by τ/Δt, where Δt is the interval between the two samples compared:
- **Species below ~10⁻¹⁰.** The pulse averages of ions and excited species (O(1D), NO2⁺, OH⁺, …) change by 10⁻³–5·10⁻² between pulses without any drift, because of the discrete EEDF updates. `steadyStateYmin 1e-10` leaves them out.
- **Excited species above 10⁻¹⁰.** Species such as O2(b1) and O2(a1) in R1 (Y ~10⁻⁹) still change by ~4·10⁻⁴ between 10-pulse averages, because of small differences between the discharges. Over Δt = 10 pulses this gives d ≈ 2·10⁻² for R1 at steady state, far above tol = 10⁻³. Over W = 100 it is ~2·10⁻³, still above the threshold. Over W = 1000 it is ~2·10⁻⁴, below.
- **Delay.** With the comparison over W, d measures the average drift over the last W pulses, so it falls below tol somewhat after the instantaneous drift does (~550 pulses for R1 with W = 1000). Steady state is then detected within W/k pulses (100 with k = 10). A counter requiring W more pulses below the threshold would add W pulses: R1 would freeze at pulse ~9070 instead of ~8200.
- **Accuracy.** With tol = 1e-3, at freezing the main species of R2 still move by ≤ 10⁻³ over one residence time.

Caveat: τ is the time scale of the mass renewal. When the chemistry of a reactor is slower (in R2 NO approaches steady state with a time constant of about 3τ), d underestimates the remaining distance. Use a tight threshold, or a longer `steadyStateTau`, in that case.

## 7. Analytic Jacobian

### 7.1 Why

CVODE (BDF, Newton) needs the Jacobian J = ∂f/∂y of the right-hand side. Cantera's integrator supports only a finite-difference Jacobian with the dense solver: N = K + 4 evaluations of the right-hand side per Jacobian (74 for the mechanism of CRN_serialPlasma). These evaluations were 58% of all the right-hand-side evaluations of R1.

### 7.2 Integrator (`cvodeJac.h`)

`CRN::CVodeJac` derives from `Cantera::Integrator` and drives CVODES directly, with the same set-up as Cantera's `CVodesIntegrator` with the `DENSE` linear solver:
- BDF, Newton, dense LAPACK solver, scalar tolerances;
- at most 20000 steps per `integrate()`;
- linear solver rebuilt at every restart.

The `step()`, `integrate()` and `derivative()` methods have the same semantics as Cantera's, so the rest of the code does not change. With `jacobian fd` the results are bit-identical to `jacobian cantera`. With `jacobian analytic`, CVODE calls `ChemPlasReactor::jacobian` through `CVodeSetJacFn`.

### 7.3 Jacobian (`ChemPlasReactor::jacobian`)

State y = [T, h_e, Ep, e_vib(N2), Y_1 … Y_K]; J is built column by column:

| Block | How |
|---|---|
| Species rows, species columns | analytic: dω/dC from Cantera (`Kinetics::netProductionRates_ddCi`: law of mass action, third bodies, falloff, any rate type, since at fixed T the EEDF and Te-dependent rate constants are constant). The chain rule through the mass fractions normalized by `setMassFractions` gives dω/dY = G·D·(I − Ŷ1ᵀ)/S, with D = diag(ρ/W). Then dY_k/dt = ω_k W_k/ρ + Σ ṁ/V (Y_src,k − Y_k)/ρ |
| T, h_e, Ep rows, species columns | analytic: derivatives of −Σu_kω_k, of P_ext = K_p c_tot C_e, of P_elec = 1.5 R T_e ω_e, of the outflow work and of c_v |
| T column | one forward difference (1 evaluation of the right-hand side): it covers the T dependence of rates, thermodynamics and flow terms |
| h_e, Ep, e_vib(N2) columns | zero: these variables do not enter the right-hand side with the Flitti heating model |

- **Negative mass fractions.** `setMassFractions` sets them to 0, so their columns have no composition derivative, only the flow term. The Jacobian check found this.
- **Neglected term.** The dependence of the heat-loss conductivity on the composition is neglected; `heatLoss` is off in practice.
- **Coverage.** The analytic Jacobian covers constant volume, the Flitti heating model, `nonThermal no` and no imposed n_e (`analyticJacobianSupported()`). In the other cases CVODE keeps its finite-difference Jacobian, and the program says so at start.

### 7.4 Validation (CRN_serialPlasma, 300 pulses, reltol 1e-9)

| | `cantera` | `analytic` |
|---|---|---|
| Right-hand-side evaluations, R1 | 601 000 | 260 000 (−57%) |
| Right-hand-side evaluations, R2 | 62 000 | 11 000 (−82%) |
| Steps, Jacobians, Newton iterations (R1) | 212 000, 4600, 260 000 | 210 000, 4500, 255 000 |
| CPU time (3 repetitions) | 11.6–11.9 s | 8.86–8.89 s (−24%) |
| Difference from reltol 1e-10, main species | ≤ 7.6·10⁻⁵ | ≤ 4·10⁻⁵ |

- **Check against finite differences.** With `jacobianCheck`, J agrees with the finite-difference Jacobian within 7.5·10⁻⁵ (scaled with the CVODE error weights). Larger differences appear only on species with Y ≈ −10⁻²⁷, where the finite-difference step crosses the clipping of `setMassFractions`.
- **Cost of the Jacobians.** Building them costs ~0.23 ms each (~12 right-hand-side evaluations, against 74), and 80% of that is `netProductionRates_ddCi`.
- **Where the time goes now.** Most of it is inside CVODE (factorizations, steps) and in the interpolation of the state at every internal step for the outflow averages.

## 8. How it works

```
read case → reactor defaults → CppBOLOS (grid, cross sections) → for every PSR: Cantera Solution,
ChemPlasReactor, CVODE, initial EEDF → time loop
```

**Time loop.** The outer loop advances one flow update interval at a time, from t_k to t_k+1 = startTime + (k+1)·`flowUpdateInterval`.

1. **Flows (at t_k).**
   - `mainSolver` is called, and each controller gets the source state (Y, h, T) averaged over the previous interval. Reservoirs with inflows are re-mixed.
   - For each reactor the new flow terms are compared with the previous ones; CVODE is restarted only if they changed.
2. **Each reactor, on its own, from t_k to t_k+1.**
   - *Discharge reactor* (or `integrateToFlowUpdate 0`): a sequence of steps as in v4 or v5.
     - `findTimeStep`: `dt1` in the first two pulses, then logarithmic steps during the discharge (`NstepInPulse`) and after it (`NstepAfPulse`, up to `pulseExten`), `dt_max` between pulses, reduced on ignition. The last step is cut at t_k+1.
     - EEDF check (table or computation) and CVODE restart only when needed.
     - During the discharge: CVODE internal steps and exact end of the discharge at `pulseEnergy`, then a new EEDF at the base E/N and a restart. Outside the discharge: CVODE internal steps.
     - Density (constant volume), output, per-pulse file, start of the new pulse (Ep reset and restart).
   - *Reactor without discharge* (`integrateToFlowUpdate 1`): a single iteration:
     - EEDF check at the start of the interval;
     - CVODE internal steps up to t_k+1; after a restart, a second restart at t_k + `dt_max`;
     - density, output and pulse counter at t_k+1.
   - In both cases the outflow of the reactor is integrated in time with the trapezoidal rule on the CVODE internal steps, and its average becomes the inflow of the downstream reactors in the next interval.
3. **End of the interval:** reservoir rows, stop check, steady-state check (section 6) with possible freezing, progress line.

In addition, in v7:
- frozen reactors are skipped in step 2;
- the release of frozen reactors is decided in step 1, by comparing the inflows;
- in reactors without discharge Ep is not reset at every pulse: the energy of a single pulse is Ep minus an offset, and CVODE is not restarted.

The one-interval lag of the coupling between reactors has no effect at periodic steady state, when the average outflow of a reactor is the same every period. In the transient it moves NO and H₂O of R2 by 0.3–0.6% compared with a coupling 4–10 times finer.

**Stop:** the end time is reached, max(dT/dt) ≤ `dTdt_min`, or (monitor on and `steadyStateStop 1`) the whole network is at steady state.

## 9. Output (`<case>/output/`)

| File | Content |
|---|---|
| `<reactor>.csv` | `Time(s), T_gas(K), p(Pa), N_gas(#/cm^3), MW(kg/kmol)`, number density [#/cm³] of each species, `Power_input(W/cm^3), Energy_deposited(mJ/cm^3)`. One row every `writeInterval` and at the last time. Reactors without discharge with `integrateToFlowUpdate 1` write only at the flow updates (with `flowUpdateInterval = 1/f_NRP`, one row per pulse) |
| `<reactor>_pulses.csv` | discharge reactors only, one row per pulse: `iPulse, Time(s)` (end of the period), `T_gas, p, n_e_end, n_e_max, Energy_deposited` (at the end of the period), `tau_NRP(s)` (discharge duration), `Energy_discharge(mJ/cm^3)` (energy at the end of the discharge) |
| `massflow.csv` | time and flow rate [kg/s] of every controller at every update |
| `<reservoir>.csv` | state of the reservoirs (same columns, without power and energy) |
| `<case>/eedfTable.yaml` | EEDF table, updated during the run (unless `eedfTable` is `none`) |
| `steadyState.csv` | with the monitor on, one row per reactor checked at every check: `Time(s), node, d_max, worst` (quantity with the largest d, or `mass`), `counter, converged` |

At the end of the run the work of each integrator is printed: right-hand-side evaluations and, for `fd`/`analytic`, steps, Jacobians, Newton iterations and failures, restarts and the time spent in the analytic Jacobians.

At the end of the run every reactor has a last row at the final time, even if it has been frozen for a long time. Freezing and release events are printed on the terminal (`Steady state at runTime = … s: R1 (no longer integrated)`, `New inflows at runTime = …`). If the run reaches the largest time without steady state, the last d and the worst quantity of every reactor are printed.

## 10. Practical notes

- **Run time without the monitor.** The full CRN_serialPlasma case (7500 pulses) takes about 310–320 s of CPU time with reltol 1e-9, about 296 s with reltol 1e-8. Over 80% of the time is the CVODE integration of the discharge reactor.
- **Finer coupling.** With `flowUpdateInterval` shorter than the period (e.g. 4 updates per pulse) the coupling between reactors becomes finer. The conditional restart avoids restarting reactors fed only by reservoirs, but the cost still grows (about +16% with 4 updates per period).
- **Run time with the monitor.** On CRN_serialPlasma (W 1000, k 10, tol 1e-3, `jacobian analytic`) the run reaches steady state in about 230 s of CPU time (299 s with `jacobian cantera`): ~8200 pulses with R1, then R2 alone up to pulse ~35 100. The phase with R2 alone costs a few seconds every 10⁴ pulses.
- **Checking the parameters.** To see how d decreases without freezing anything, run with an unreachable threshold (`steadyStateTol 1e-30`) and analyse `steadyState.csv`.
- **Timing under WSL.** Times measured with `/usr/bin/time` jitter by ±1–2 s. Compare runs with the CPU time printed at the end.
- **Species outside CppBOLOS.** If a species not listed in `BoltzmannSpecies` exceeds a mole fraction of 1e-2, a warning is printed.
- **Unsupported mechanisms.** Mechanisms with `branch_ratio` are rejected, because the Boltzmann rate constants stored per reactor do not support it. Use CRN_v3 for them.
- **Sharing the table.** The table depends only on the cross sections, the grid and the convergence settings, not on the case: it can be copied between cases with the same settings.
