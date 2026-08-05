# Cell Passport / Digital Twin Dashboard

Static browser dashboard for replaying an uploaded cycler CSV through the repository's identified SHP320-35 cell model.

## Run locally

From the repository root:

```bash
python3 -m http.server 8080 --bind 0.0.0.0
```

Then open:

```text
http://localhost:8080/dashboard/
```

The page loads model artifacts directly from `results/`:

- `LUT2D_Charge_2RC.csv`
- `LUT2D_Discharge_2RC.csv`
- `Capacity_vs_T.csv`

## CSV input

The uploaded file should contain the same cycler columns used by the MATLAB workflow:

- `DWell Time(ms)` or another time column
- `Voltage(V)`
- `Current(A)`
- optional `Temperature`

Current sign convention follows the MATLAB project: `I > 0` is charge and `I < 0` is discharge.

## What it shows

- measured voltage vs simulated voltage from the dual-LUT 2-RC ECM
- voltage residual and error metrics
- estimated SOC trajectory
- SOC-domain charge/discharge passport curves for OCV, R0, R1, R2, Tau1, Tau2, C1, and C2, with the selected parameter controlled in the passport panel
- sample-level readout of the uploaded profile controlled by the slider

The default SOC mode is the same HPPC/full-charge anchor used in the MATLAB validation scripts. For arbitrary drive cycles, switch to user initial SOC mode.
