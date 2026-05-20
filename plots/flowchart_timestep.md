# ChemPlasKin CRN — findTimeStep Flowchart

`findTimeStep` is called **once per reactor** per main time-loop iteration.
It modifies `dt` (passed by reference) and `EN_temp`.
The outer loop then takes `dt = min(dt, dt_element)` across all reactors.

```mermaid
flowchart TD
    START(["findTimeStep(state, config, runTime, ..."]) --> INIT["dt_ig = dt_max\nEN_temp = max(EN_DC, 0.1)"]

    INIT --> IGN{"dT > dT_max\nAND dTdt > 1000?"}
    IGN -->|"Yes — ignition detected"| IGDT["dt_ig = min(dT_max / dTdt, dt_max)\nPrint IGNITION DETECTED"]
    IGN -->|"No"| PCHECK
    IGDT --> PCHECK

    PCHECK{"iPulse < nPulses?\n(still within NRP pulse sequence)"}

    %% ---- NO MORE PULSES ----
    PCHECK -->|"No — all pulses done"| NOPULSE["if not stayEN_DC: EN_temp = 0.1\ndt = dt_max"]
    NOPULSE --> CLAMP

    %% ---- STILL PULSING ----
    PCHECK -->|"Yes"| AVOID["dt3 = (iPulse + 1) × pulsePeriod − runTime\n[avoid overstepping next pulse start]\ndt_ig = min(dt3, dt_ig)"]

    AVOID --> PEXT{"runTime + SMALL < iPulse × pulsePeriod\n+ pulseExten?\n(within pulse or fast relaxation window)"}

    PEXT -->|"No — between pulses"| BTWN["dt = min(dt3, dt_max)"]
    BTWN --> CLAMP

    PEXT -->|"Yes"| PDISCH{"runTime + SMALL < iPulse × pulsePeriod\n+ pulseWidth\nAND dischargeOn\nAND config.discharge?"}

    %% ---- DURING DISCHARGE ----
    PDISCH -->|"Yes — active discharge"| EARLY{"iPulse <= 1\nOR nE_1 < 1e10?"}
    EARLY -->|"Yes — early pulses"| FIXDT["dt = dt1  (fixed discharge step)"]
    EARLY -->|"No"| LOGDT["K = max(round(log10(nE_2/nE_1)×100)/100, 0.2)\ndt = logDynamicTimestep(K, tau_NRP,\n     NstepInPulse, iStepInPulse)"]
    FIXDT --> ENSTEP["EN_temp = EN\niStepInPulse++"]
    LOGDT --> ENSTEP
    ENSTEP --> CLAMP

    %% ---- AFTER DISCHARGE (pulse extension) ----
    PDISCH -->|"No — discharge off, in pulseExten"| AFIRST{"iStepAfPulse == 0?\n(first step after discharge ends)"}

    AFIRST -->|"Yes — record end-of-pulse state"| REINIT["integrator->reinitialize(runTime, odes)\ntau_NRP = runTime − iPulse × pulsePeriod\nnE_2 = n_e_current\nfastExpansion = true"]
    REINIT --> KDECAY

    AFIRST -->|"No"| KDECAY

    KDECAY["K_decay:\n  if iPulse<=1 or nE_2<1e10: K = -1  (fixed, slow decay)\n  else: K = min(round(log10(nE_3/nE_2)×100)/100, -0.1)"]

    KDECAY --> DTDECAY{"iStepAfPulse >= NstepAfPulse?"}
    DTDECAY -->|"Yes — plateau: dt unchanged"| HSCHECK
    DTDECAY -->|"No"| LOGAFT["dt = logDynamicTimestep(K, pulseExten,\n     NstepAfPulse, iStepAfPulse)"]
    LOGAFT --> HSCHECK

    HSCHECK{"heatLoss AND fastExpansion\nAND runTime > iPulse×pulsePeriod + 1e-7?"}
    HSCHECK -->|"Yes — isentropic expansion"| ISENTR["T_new = T × (P0/P)^(1 − cv/cp)\ngas->setState_TPY(T_new, P0, Y)\nodes->setConstPD(P0, rho_new)\nif CP: constPressure = true\nfastExpansion = false  (only once per pulse)"]
    ISENTR --> STEPACC
    HSCHECK -->|"No"| STEPACC

    STEPACC["iStepAfPulse++\nif iStepAfPulse == NstepAfPulse:\n  nE_3 = n_e_current"]
    STEPACC --> CLAMP

    CLAMP["dt = max(dt, 1e-14)\ndt = min(dt, dt_ig)"]
    CLAMP --> TINY{"dt < 1e-16?"}
    TINY -->|"Yes"| ABORT(["ABORT — dt cascade"])
    TINY -->|"No"| RET(["Return (dt updated by reference)"])
```

## Pulse phase timeline

```
 |<-- pulseWidth -->|<--- pulseExten (total) ------>|<--- inter-pulse gap --->|
 [  DISCHARGE ON   ][  DISCHARGE OFF + relaxation  ][     normal dt_max      ]
  iStepInPulse++     iStepAfPulse++
  EN_temp = EN       EN_temp -> DC
  dt = logDynamic    dt = logDynamic (decay K)
  (growth K)
```

## logDynamicTimestep behaviour

| Parameter | Meaning |
|---|---|
| `K > 0` | Logarithmically growing dt during discharge (electron density rises) |
| `K < 0` | Logarithmically decaying dt during relaxation (electron density falls) |
| `K = -1` | Fixed early-pulse fallback (first two pulses or very low nE) |
| `iStepAfPulse >= NstepAfPulse` | Plateau: dt frozen at last value to prevent dt runaway |
