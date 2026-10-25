# CRN_v4 – ChemPlasKin for reactor networks

ChemPlasKin solves coupled gas–plasma kinetics in 0D reactors:
- the chemistry comes from Cantera;
- the electron energy distribution function (EEDF) comes from the CppBOLOS Boltzmann solver;
- the integrator is CVODE.

`CRN_v4` extends the code to a **chemical reactor network** (CRN). Perfectly stirred reactors (PSR) are connected by mass flow controllers. The flow rates are computed at every update by a Python script of the case. One or more reactors can have a nanosecond repetitively pulsed (NRP) discharge.

Compared with the original version (CRN_v3), v4 adds:

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
CRN_v4/
├── masterCRN.cpp       main program: case input, network, time loop
├── plasmaReactor.h     ChemPlasReactor class: gas-plasma reactor equations (right-hand side for CVODE)
├── CMakeLists.txt      build
├── controlDict         example case (time control)
├── chemPlasProperties  example case (physics and plasma)
├── crnDict             example case (reactors and mass flow controllers)
├── massflowsolver.py   example flow-rate solver
├── massflow_params.txt parameters read by massflowsolver.py
├── mechanisms.txt      notes on the available mechanisms
└── build/              ChemPlasKin executable (after the build)
```

Dependencies, with paths relative to the `CRN_v4` folder:
- `../../../cantera`: the ChemPlasKin fork of Cantera (headers in `include/`, libraries in `build/lib`, `libcantera_shared`), which contains CppBOLOS;
- `../../src/utilities.h`: `readParameter` and the reader of OpenFOAM-style files;
- Python 3 with development files (`Python.h`): the flow-rate solver is called from C++ through the embedded interpreter;
- for `massflowsolver.py`: `numpy` and `scipy`.

## 2. Build

```bash
cd CRN_v4
mkdir -p build && cd build
cmake ..
make            # builds build/ChemPlasKin
```

## 3. Run

```bash
cd CRN_v4/build
./ChemPlasKin -case <case_folder> [-log]
```

- `-case`: folder with the case files. Without `-case`, `..` is used, i.e. `CRN_v4/`.
- `-log`: detailed output at every step and CppBOLOS info messages. Very slow: for debugging only.

During the run a progress line is printed (`runTime`, pulse, number of steps). At the end the program prints the reason for stopping and the CPU time.

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

All dictionary files use the OpenFOAM-style syntax: `key value;`, blocks `name { ... }`, comments `//` and `/* */`.

### 4.3 `controlDict`

| Key | Req. | Unit | Meaning |
|---|---|---|---|
| `startTime` | yes | s | start time (in `autoplay`, the end time of the Cantera phase) |
| `endTime` | yes | s | end time. The run ends at max(`endTime`, `startTime` + `nPulses`/`f_NRP`) |
| `nPulses` | yes | – | number of pulses |
| `endTimeForce` | no | s | if given, upper bound of the end time |
| `dt1` | yes | s | fixed step during the discharge of the first two pulses |
| `dt_max` | yes | s | largest step |
| `pulseExten` | yes | s | length of the fine-step phase from the start of each pulse. Clipped between `pulseWidth` and the period |
| `NstepInPulse` | yes | – | number of (logarithmic) steps during the discharge, from the third pulse on |
| `NstepAfPulse` | yes | – | number of (logarithmic) steps from the end of the discharge up to `pulseExten` |
| `dT_max` | yes | K | ignition detection threshold: if dT > `dT_max` and dT/dt > 1000 K/s the step is reduced. At least 1 K |
| `dTdt_min` | yes | K/s | stop when max(dT/dt) ≤ `dTdt_min`. A negative value disables the stop |
| `kOutPressure` | yes | – | read but not used (left over from the Cantera pressure controller) |
| `flowUpdateInterval` | no | s | flow update interval. `0` (default) = every step, with the instantaneous source state. `> 0` = every interval, with the time-averaged source state. Recommended: `1/f_NRP` |
| `writeInterval` | no | s | interval between output rows. `0` (default) = every step |
| `reinitEveryStep` | no | 0/1 | `1` = restart CVODE at every step, as in v3. Default `0` |
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

## 5. How it works

```
read case → reactor defaults → CppBOLOS (grid, cross sections) → for every PSR: Cantera Solution,
ChemPlasReactor, CVODE, initial EEDF → time loop
```

**Time loop.** All reactors advance together, with the same step. At every step:
1. **Flows.** If an update is due, `mainSolver` is called. Each controller gets the source state (Y, h, T) averaged over the past interval. Reservoirs with inflows are re-mixed, and the inflow and outflow terms enter the equations of every reactor.
2. **Step.** `findTimeStep` gives a step for each reactor, and the common step is the minimum:
   - during the discharge: `dt1` in the first two pulses, then logarithmic steps (`NstepInPulse`) scaled on the growth of n_e;
   - after the discharge, up to `pulseExten`: logarithmic steps (`NstepAfPulse`) scaled on the decay of n_e;
   - between pulses: `dt_max`, without crossing the start of the next pulse;
   - reduced on ignition (`dT_max`).
3. **EEDF.** For each reactor, if E/N, T or the mixture of the Boltzmann species changed beyond the tolerances, the Boltzmann equation is solved again and the new rate constants are installed; then CVODE is restarted.
4. **Integration.**
   - Discharging reactor: CVODE internal steps (`step()`), watching the deposited energy Ep. When Ep reaches `pulseEnergy`, the time t* is located on the interpolant. Then E/N switches to the base value, a new EEDF is computed, CVODE is restarted and the integration continues to the end of the step.
   - Other reactors: `integrate()`.
5. **After the step.**
   - in constant-volume reactors, the density is updated with the net flow;
   - power and deposited energy, output;
   - at the end of a period: row in `_pulses.csv`, Ep reset and start of the new pulse;
   - the time integrals of the outflow are accumulated (state at the end of the step × dt). At the next update they become the average.

**Stop:** the end time is reached, or max(dT/dt) ≤ `dTdt_min`.

## 6. Output (`<case>/output/`)

| File | Content |
|---|---|
| `<reactor>.csv` | `Time(s), T_gas(K), p(Pa), N_gas(#/cm^3), MW(kg/kmol)`, number density [#/cm³] of each species, `Power_input(W/cm^3), Energy_deposited(mJ/cm^3)`. One row every `writeInterval` and at the last time |
| `<reactor>_pulses.csv` | discharge reactors only, one row per pulse: `iPulse, Time(s)` (end of the period), `T_gas, p, n_e_end, n_e_max, Energy_deposited` (at the end of the period), `tau_NRP(s)` (discharge duration), `Energy_discharge(mJ/cm^3)` (energy at the end of the discharge) |
| `massflow.csv` | time and flow rate [kg/s] of every controller at every update |
| `<reservoir>.csv` | state of the reservoirs (same columns, without power and energy) |

## 7. Practical notes

- **Run time.** The full CRN_serialPlasma case (7500 pulses) takes about 500 s of CPU time. Over 80% of it is the CVODE integration of the discharge reactor.
- **Timing under WSL.** Times measured with `/usr/bin/time` jitter by ±1–2 s. Compare runs with the CPU time printed at the end.
- **Species outside CppBOLOS.** If a species not listed in `BoltzmannSpecies` exceeds a mole fraction of 1e-2, a warning is printed.
- **Unsupported mechanisms.** Mechanisms with `branch_ratio` are rejected, because the Boltzmann rate constants stored per reactor do not support it. Use CRN_v3 for them.
