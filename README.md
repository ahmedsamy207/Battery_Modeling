# SHP320-35 Hybrid Solid–Liquid Cell — 2-D (SOC × T) 2-RC ECM & LUT pipeline

Extension of the original 25 °C 1-D pipeline to **temperature-dependent 2-D
look-up tables** `param = f(SOC, T)`. The pipeline keeps the exact extraction
methodology of the 25 °C scripts (GITT-rested OCV, 10 s 2-RC `fminsearch`
pulse fits, true-capacity Coulomb counting) and generalises every stage to a
list of temperatures.

```
SHP320_2RC_2D/
├── data/                                   ← raw cycler CSVs (all temps)
├── matlab/
│   ├── testConfig.m                        ← the ONLY file you edit (test matrix)
│   ├── step1_Capacity_and_OCV_2D.m         ← Q_Ah(T), OCV chg/dis/avg (SOC×T)
│   ├── step2_Extract_Perfect_2RC_2D.m      ← 2-D averaged LUT  (chg+dis mean)
│   ├── step2_Extract_Perfect_2RC_DualLUT_2D.m ← 2-D charge & discharge LUTs
│   ├── step3_Validation_2D.m               ← validates averaged LUT per T
│   ├── step3_Validation_2D_DualLUT.m       ← validates dual LUTs per T
│   ├── step4_Build_SimscapeModel.m         ← Simscape cell block (dual LUT + 1-state hysteresis)
│   ├── step4_Validate.m                    ← validates dual LUTs for the cell block
│   ├── step5_BuildPackParams.m             ← Scale to 192S4P pack
│   ├── step5_Validate.m                    ← validates the pack
│   ├── Capacity_OCV_Validate.m             ← Capacity/OCV replay on cell model
│   └── Aircraft_Mission.m                  ← NASA X-57 / eVTOL mission on pack
│                                           (also exports mission_<type>.csv)
├── dashboard/
│   ├── index.html                          ← cell passport / digital twin dashboard
│   ├── mission.html                        ← aircraft mission dashboard (192S4P pack, playback)
│   └── README.md
├── tools/verify_pipeline_python.py         ← numerical mirror used for testing
└── results/                                ← all outputs are written here
```

## Browser dashboard

Static, GitHub Pages-compatible dashboards live in `dashboard/` (plain HTML +
inline CSS/JS, Plotly from CDN — no build system or backend):

- **Cell dashboard** (`dashboard/index.html`) — the Cell Passport / Digital Twin
  dashboard. Loads the generated dual-LUT model from `results/`, accepts an
  uploaded cycler CSV, replays it through the time-domain 2-RC ECM, and plots
  measured vs simulated voltage, residuals, SOC, parameters along the uploaded
  profile, and dedicated SOC-domain Cell Passport LUT curves.
- **Mission dashboard** (`dashboard/mission.html`) — aircraft mission /
  pack-load visualization for the **192S4P pack** with **dynamic time playback**.
  Supports NASA X-57 style and eVTOL demo mission profiles (phase data from
  `matlab/Aircraft_Mission.m`), plotted against pack voltage, current, power,
  SOC, altitude, range, C-rate and energy, with aircraft visuals, KPI/telemetry
  readouts, a phase table, and a mission timeline. Playback includes a
  Play/Pause/Reset control, a draggable timeline slider, and a **playback speed
  factor** (`1×`, `5×`, `10×`, `25×`, `50×`, `100×`) that only affects dynamic
  playback. Without a MATLAB export the profiles are clearly labeled **demo
  mission profiles** with **browser-side pack estimates**; `Aircraft_Mission.m`
  can export exact `results/mission_<type>.csv` profiles that the dashboard then
  replays as-is.

Run them from the repository root with:

```bash
python3 -m http.server 8080 --bind 0.0.0.0
```

Then open `http://localhost:8080/dashboard/` (cell) or
`http://localhost:8080/dashboard/mission.html` (mission).

### GitHub Pages hosting

The dashboard is ready to host from the repository root with GitHub Pages:
`index.html` redirects to `dashboard/`, and the dashboard reads the required
CSV artifacts from `results/` plus the bundled 25 °C HPPC demo from `data/`.

After this branch is merged, enable Pages in the repository settings with:

- Source: **Deploy from a branch**
- Branch: `main`
- Folder: `/ (root)`

The dashboard will then be available at:

```text
https://ahmedsamy207.github.io/Battery_Modeling/
```

## How to run

1. Put the four characterization CSVs of each temperature into `data/`,
   keeping the naming convention
   `Solid State - Cell 2 - {Capacity|OCV Charge|OCV Discharge|HPPC} <T> Celsius.csv`.
2. Register the temperature in `matlab/testConfig.m` (`Tlist = [15 25 45];`
   is pre-configured). **Missing files are skipped automatically** — the
   pipeline runs with any subset (e.g. 15 °C + 25 °C today, re-run when the
   45 °C files arrive).
3. Run in MATLAB, in order:
   ```matlab
   step1_Capacity_and_OCV_2D
   step2_Extract_Perfect_2RC_2D
   step2_Extract_Perfect_2RC_DualLUT_2D
   step3_Validation_2D
   step3_Validation_2D_DualLUT
   ```

> Run order matters: step 2 needs `results/True_OCV_Capacity_2D.mat`,
> step 3 needs `results/LUT2D_*.mat`.

## LUT format

All 2-D LUTs are MATLAB structs saved with `-struct` (so `load` gives you the
fields directly), oriented exactly as Simulink / Simscape 2-D lookup tables
expect — **rows = SOC breakpoints, columns = temperature breakpoints**:

| field       | size        | meaning                                   |
|-------------|-------------|-------------------------------------------|
| `soc_grid`  | 201 × 1     | SOC breakpoints (0:0.005:1)               |
| `Tbp`       | 1 × nT      | temperature breakpoints [°C]              |
| `Q_Ah`      | 1 × nT      | true capacity vs temperature              |
| `OCV_V`     | 201 × nT    | OCV (avg / dis-branch / chg-branch)       |
| `R0_ohm`, `R1_ohm`, `R2_ohm` | 201 × nT | resistances                    |
| `tau1_s`, `tau2_s`           | 201 × nT | time constants                   |
| `C1_F`, `C2_F`               | 201 × nT | capacitances                     |

Human-readable wide CSVs (`LUT2D_*_2RC.csv`, one column-block per
temperature) and raw per-pulse extraction CSVs (`PerPulse_RAW_<T>C.csv`) are
written alongside for inspection in Excel.

In Simulink use a *2-D Lookup Table* with breakpoint vectors `{soc_grid, Tbp}`
and table data `LUT2D.R0_ohm` (etc.); set extrapolation to **Clip**.

## What changed vs. the 25 °C scripts (and why)

1. **Test matrix instead of hard-coded files** — every script loops over
   `testConfig.m`; SOC tracking uses each temperature's own `Q_Ah`.
2. **2-D lookup in validation** — `griddedInterpolant` over `(SOC, T)` with
   *clamped* extrapolation in T (so you can simulate temperatures between and
   slightly outside the measured breakpoints).
3. **Safe grid extrapolation (bug fix).** The 25 °C scripts used
   `interp1(...,'linear','extrap')`, i.e. *slope* extrapolation of noisy pulse
   fits. At 15 °C the first +70 A charge pulse is aborted (10 ms), leaving no
   charge pulse between 85 → 100 % SOC; the slope extrapolation of R2/C2 there
   blew the RC states up to **~3 V error**. The safe rule is:
   * linear interpolation inside the measured SOC coverage,
   * **constant nearest-edge extension** outside it,
   * OCV: where one hysteresis branch has no data, use the other branch
     (hysteresis collapses at the SOC edges) — no artificial slope overshoot.
4. Pulse-selection thresholds (`maxDur 20 s`, `|I| > 40 A`, `min 5 samples`,
   10 s fit window) are now parameters at the top of each step-2 script — the
   same values that worked at 25 °C, validated to catch exactly the
   18 dis + 16 chg valid pulses at 15 °C (the 10 ms aborted pulse is rejected
   by the ≥5-sample rule).

## Verified numbers on the uploaded 15 °C + 45 °C data

(tools/verify_pipeline_python.py — an exact numerical mirror of the MATLAB code;
25 °C is skipped here because those raw CSVs are only on your machine — on your
PC all three temperatures will be picked up automatically)

| stage | 15 °C | 45 °C |
|---|---|---|
| True capacity Q | **33.8961 Ah** | **35.2611 Ah** |
| OCV rest points | 21 chg + 21 dis | 21 chg + 21 dis |
| HPPC pulses extracted | 34 (18 dis + 16 chg) | 37 (20 dis + 17 chg) |
| R0 range | 2.69–3.17 mΩ | **1.57–1.94 mΩ** |
| Validation, averaged 2-D LUT | RMSE 64.0 mV | RMSE 38.2 mV |
| Validation, **dual 2-D LUT** | **RMSE 34.5 mV** (chg 33.3 / dis 51.6) | **RMSE 11.2 mV** (chg 22.3 / dis 12.7) |

Physics captured by the 2-D tables:
* **Capacity grows +1.37 Ah (≈ 4 %) from 15 → 45 °C** — each temperature tracks
  SOC with its own Q.
* **R0 drops ≈ 45 %** from 15 → 45 °C (thermally activated kinetics of the
  hybrid electrolyte); R1 shows the classic mid-SOC bowl with a much deeper
  cold-temperature rise; R2 (diffusion) rises sharply below ~30 % SOC,
  especially at 15 °C.
* OCV shows a clear hysteresis loop (≈ 30–60 mV mid-SOC) and a mild
  temperature offset, both absorbed into `OCV_V(SOC, T)`.

Residual structure: ~10–25 mV in the 40–100 % SOC band, growing toward low
SOC where R(SOC) curvature between sparse pulse points dominates
(worst ≈ 220 mV on the 70 A pulse at 9 % SOC, 15 °C). This is the expected
behaviour of an HPPC-trained 2-RC ECM.

Result artifacts: `results/verify_LUT2D_*.csv` (SOC × {15 °C, 45 °C} tables),
`verify_PerPulse_RAW_<T>C.csv`, and plots `verify_step1_OCV.png`,
`verify_step2_params.png`, `verify_step3_validation_ALL_T.png`.

## When the 25 °C files are added on your machine

Keep the 4 CSVs in `data/` (names already match `testConfig.m`) and re-run
steps 1→3 — all LUTs get their third column and you get a 3-row RMSE report.
