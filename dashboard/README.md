# SHP320-35 Dashboards (Cell Passport / Digital Twin + Aircraft Mission)

Static, GitHub Pages-compatible browser dashboards for the repository's SHP320-35
cell model and its 192S4P battery pack. Both pages are plain HTML + inline CSS/JS
with Plotly from CDN — no build system, no backend.

| Page | URL | What it does |
|---|---|---|
| Cell dashboard | `dashboard/index.html` | Replays an uploaded cycler CSV through the dual-LUT 2-RC ECM, plots measured vs simulated voltage, residuals, SOC and LUT parameters |
| Mission dashboard | `dashboard/mission.html` | Aircraft mission / pack-load visualization for the 192S4P pack (X-57 & eVTOL demo missions) with dynamic time playback |

Top navigation on both pages links between them.

## Run locally

From the repository root:

```bash
python3 -m http.server 8080 --bind 0.0.0.0
```

Then open:

```text
http://localhost:8080/dashboard/            # cell dashboard
http://localhost:8080/dashboard/mission.html  # aircraft mission dashboard
```

## Cell dashboard (`index.html`)

The page loads model artifacts directly from `results/`:

- `LUT2D_Charge_2RC.csv`
- `LUT2D_Discharge_2RC.csv`
- `Capacity_vs_T.csv`

### CSV input

The uploaded file should contain the same cycler columns used by the MATLAB workflow:

- `DWell Time(ms)` or another time column
- `Voltage(V)`
- `Current(A)`
- optional `Temperature`

Current sign convention follows the MATLAB project: `I > 0` is charge and `I < 0` is discharge.

### What it shows

- measured voltage vs simulated voltage from the dual-LUT 2-RC ECM
- voltage residual and error metrics
- estimated SOC trajectory
- charge/discharge values of OCV, R0, R1, R2, Tau1, Tau2, C1, and C2 evaluated along the uploaded profile
- dedicated Cell Passport LUT Curves showing OCV/ECM parameters versus SOC at a selected temperature
- sample-level readout controlled by a slider

The default SOC mode is the same HPPC/full-charge anchor used in the MATLAB validation scripts. For arbitrary drive cycles, switch to user initial SOC mode.

## Aircraft mission dashboard (`mission.html`)

New page for mission / pack-load visualization of the **192S4P pack** (192S × 4P,
≈ 140 Ah, ≈ 700.8 V, ≈ 98.1 kWh — scaling per `matlab/step5_BuildPackParams.m`).

### Supported mission profiles

- **X-57 Maxwell (cruise mission)** — phase powers/durations from the NASA X-57
  profile (Battery_Evaluation_EATS Table 1) as scripted in `matlab/Aircraft_Mission.m`:
  taxi → motor check → takeoff → climb → cruise → descent (≈ 25 min).
- **eVTOL / UAM (stress mission)** — representative eVTOL profile from
  `matlab/Aircraft_Mission.m`: idle → taxi → hover takeoff (2C peak) → climb →
  cruise → descent → hover landing (≈ 30 min).

Profiles are **demo mission profiles**: without a MATLAB export they are generated
client-side at 1 s resolution from the phase definitions, with demo altitude /
range / speed assumptions (X-57: 8,000 ft cruise at ≈ 172 kt; eVTOL: ≈ 300 m
cruise at ≈ 60 m/s). Pack quantities (voltage, current, SOC, energy) are a
**browser-side pack estimate**: current from power / V_conv (750 V, MATLAB
convention), SOC by Coulomb counting with the per-temperature pack capacity from
`results/Capacity_vs_T.csv`, OCV from `results/OCV_2D_SOCxT.csv` and an R0-only
voltage estimate (no RC dynamics, no Simscape). The page is clearly labeled as a
static browser estimate and lists all assumptions in the "Assumptions & method"
panel.

### Exact MATLAB export (optional)

`matlab/Aircraft_Mission.m` now exports `results/mission_x57.csv` and
`results/mission_evtol.csv` (columns: `t_s, phase, P_demand_kW, V_pack, I_pack,
SOC_pct, E_consumed_kWh, alt_m, range_km, v_mps, C_rate, E_remaining_kWh`).
When one of these files is present, the dashboard uses the exact export instead
of the built-in demo profile and disables the assumption controls (exported
profiles are fixed). Delete the CSV to go back to the interactive demo profile.

### Dynamic playback

- **Play / Pause / Reset** controls plus a timeline slider and current mission
  time readout.
- **Playback speed selector**: `1× real time`, `5×`, `10×`, `25×`, `50×`, `100×`.
  The factor applies to dynamic playback only; the slider can always be dragged
  manually to any mission time.
- Playback uses a `requestAnimationFrame` loop; plots and readouts follow the
  current mission time with a moving vertical cursor (Plotly shapes updated via
  `relayout`, throttled) and live KPI / telemetry updates.
- Playback stops automatically at mission end.
- Assumption controls: ambient/pack temperature, initial SOC, reserve SOC floor,
  load factor — changing them rebuilds the demo profile and resets to t = 0.

### What it shows

- aircraft cards with inline SVG silhouettes (X-57, eVTOL) that double as the
  mission selector
- whole-mission KPIs (duration, range, start/final SOC, energy, peak power,
  peak current, min voltage, reserve energy)
- live telemetry readout at the current mission time
- Plotly charts: pack voltage, pack current, pack power, SOC, altitude, range,
  C-rate, energy (consumed + remaining), aircraft speed, and a phase timeline —
  all with a current-time cursor and phase-colored backgrounds
- mission phase table (start/end/duration, nominal & peak power, altitude
  change, range, energy, notes) with the active phase highlighted
- progress visuals: aircraft marker moving along the route, altitude mini-profile
  with position dot, phase bar and phase chips

## GitHub Pages

Both dashboards are hosted from the repository root (`index.html` redirects to
`dashboard/`) and only fetch local artifacts from `results/` and `data/`, so the
whole site works under GitHub Pages without any server.
