#!/usr/bin/env python3
"""
verify_pipeline_python.py
Exact numerical mirror of the MATLAB 2D pipeline (step1 -> step2 -> step3),
used to sanity-check the algorithms on the uploaded 15 C data.
It produces results/verify15C_* artifacts (plots + CSVs) that match what the
MATLAB scripts will generate.
"""
import os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.optimize import minimize

ROOT   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA   = os.path.join(ROOT, 'data')
RESULT = os.path.join(ROOT, 'results')
os.makedirs(RESULT, exist_ok=True)

CELL = 'Solid State - Cell 2 - '
TLIST = [15, 25, 45]                      # same as testConfig.m
TESTS = []
for Tk in TLIST:
    f = {k: DATA + f'/{CELL}{k2} {Tk} Celsius.csv'
         for k, k2 in [('capacity','Capacity'), ('ocv_chg','OCV Charge'),
                        ('ocv_dis','OCV Discharge'), ('hppc','HPPC')]}
    if all(os.path.exists(v) for v in f.values()):
        TESTS.append(dict(T=Tk, **f))
    else:
        print(f'  (skipping {Tk}C: files not found)')
print('Temperatures found:', [t['T'] for t in TESTS])

SOC_GRID = np.arange(0, 1.0001, 0.005)          # 201 points

# ---------- column helpers (0-based equivalents of MATLAB T{:,n}) ----------
COL_STEP, COL_MODE, COL_T, COL_V, COL_I = 0, 4, 6, 7, 8

def load(f):
    T = pd.read_csv(f)
    cols = T.columns
    return dict(step=T[cols[COL_STEP]].to_numpy(),
                mode=T[cols[COL_MODE]].to_numpy(),
                t=T[cols[COL_T]].to_numpy()/1000.0,
                V=T[cols[COL_V]].to_numpy(),
                I=T[cols[COL_I]].to_numpy())

def cumAh(d):
    dt = np.diff(d['t'], prepend=0.0); dt[dt < 0] = 0
    return np.cumsum(d['I']*dt)/3600.0

def segments(d):
    step = d['step']
    dS = np.where(np.diff(step) != 0)[0]
    sStart = np.concatenate([[0], dS+1]); sEnd = np.append(dS, len(step)-1)
    return sStart, sEnd

# ----------------------------- STEP 1 --------------------------------------
def get_rested_ocv(f, Q, is_charge):
    d = load(f)
    Ah = cumAh(d)
    soc = (Ah - Ah.min())/Q if is_charge else 1 + (Ah - Ah.max())/Q
    _, sEnd = segments(d)
    isRest = d['mode'][sEnd] == 'REST'
    restIdx = sEnd[isRest]
    soc_pts, v_pts = soc[restIdx], d['V'][restIdx]
    u, ii = np.unique(soc_pts, return_index=True)
    return u, v_pts[ii]

print('--- Step 1 (2D): True Capacity & OCV per temperature ---')
Tbp, Qs = [], []
Vchg_all, Vdis_all, Vavg_all = [], [], []
for ts in TESTS:
    d = load(ts['capacity']); Ah = cumAh(d); Q = Ah.max()-Ah.min()
    sc, vc = get_rested_ocv(ts['ocv_chg'], Q, True)
    sd, vd = get_rested_ocv(ts['ocv_dis'], Q, False)
    # interpolate with NaN outside, then cross-fill missing branch regions
    # and flat-extend to the grid edges (NO slope extrapolation: see README)
    vcg = np.interp(SOC_GRID, sc, vc, left=np.nan, right=np.nan)
    vdg = np.interp(SOC_GRID, sd, vd, left=np.nan, right=np.nan)
    def filln(x):
        x = x.copy(); m = np.isnan(x); i = np.arange(len(x))
        x[m] = np.interp(i[m], i[~m], x[~m]); return x  # nearest fill of edges
    vavg = np.nanmean(np.vstack([vcg, vdg]), axis=0)     # uses single branch where other is NaN
    vavg = filln(vavg)
    vcg2 = vcg.copy(); m = np.isnan(vcg2); vcg2[m] = vavg[m]; vcg2 = filln(vcg2)
    vdg2 = vdg.copy(); m = np.isnan(vdg2); vdg2[m] = vavg[m]; vdg2 = filln(vdg2)
    vcg, vdg = vcg2, vdg2
    Tbp.append(ts['T']); Qs.append(Q)
    Vchg_all.append(vcg); Vdis_all.append(vdg); Vavg_all.append(vavg)
    print(f"  T={ts['T']}C: Q_Ah={Q:.4f} Ah | {len(sc)} chg + {len(sd)} dis rest points "
          f"| SOC coverage chg {sc.min():.3f}-{sc.max():.3f}, dis {sd.min():.3f}-{sd.max():.3f}")

Tbp = np.array(Tbp); Qs = np.array(Qs)
Vchg_all = np.array(Vchg_all).T; Vdis_all = np.array(Vdis_all).T; Vavg_all = np.array(Vavg_all).T

plt.figure(figsize=(8,5.5))
for j,Tk in enumerate(Tbp):
    plt.plot(SOC_GRID*100, Vavg_all[:,j], lw=2, label=f'True OCV @ {Tk}°C')
    plt.plot(SOC_GRID*100, Vchg_all[:,j], '--', lw=0.8)
    plt.plot(SOC_GRID*100, Vdis_all[:,j], ':',  lw=0.8)
plt.xlabel('SOC [%]'); plt.ylabel('OCV [V]'); plt.grid(alpha=.3)
plt.title('Step 1 (2D): True thermodynamic OCV'); plt.legend(); plt.tight_layout()
plt.savefig(f'{RESULT}/verify_step1_OCV.png', dpi=140); plt.close()

# ----------------------------- STEP 2 --------------------------------------
MAXDUR, MINI, MINN, FITWIN = 20, 40, 5, 10.0

def fitRC(tt, Vv, OCV0, Ibar, R0):
    def model(p):
        return OCV0 + Ibar*(R0 + np.exp(p[0])*(1-np.exp(-tt/np.exp(p[1])))
                                + np.exp(p[2])*(1-np.exp(-tt/np.exp(p[3]))))
    p0 = np.log([2e-4, 0.6, 9e-4, 11])
    r = minimize(lambda p: np.sum((model(p)-Vv)**2), p0, method='Nelder-Mead',
                 options={'xatol':1e-12,'fatol':1e-16,'maxfev':400000,'maxiter':40000})
    p = r.x
    return dict(R1=np.exp(p[0]), tau1=np.exp(p[1]), R2=np.exp(p[2]), tau2=np.exp(p[3]))

print('\n--- Step 2 (2D): HPPC pulse extraction ---')
perPulse = []
for j, ts in enumerate(TESTS):
    d = load(ts['hppc']); Ah = cumAh(d)
    SOC = 1 - (Ah.max()-Ah)/Qs[j]
    sStart, sEnd = segments(d)
    segMode = d['mode'][sStart]; segDur = d['t'][sEnd]-d['t'][sStart]
    segN = sEnd-sStart+1
    segI = np.array([d['I'][a:b+1].mean() for a,b in zip(sStart,sEnd)])
    SOCstart = SOC[sStart]
    isPulse = (segDur < MAXDUR) & (np.abs(segI) > MINI) & (segN >= MINN) & \
              np.isin(segMode, ['CC Discharge','CC Charge'])
    rows = []
    for k in np.where(isPulse)[0]:
        a,b = sStart[k], sEnd[k]
        tp = d['t'][a:b+1]-d['t'][a]; Vp = d['V'][a:b+1]; Ip = d['I'][a:b+1]
        sel = tp <= FITWIN; tp, Vp, Ip = tp[sel], Vp[sel], Ip[sel]
        if len(tp) < MINN: continue
        if k > 0 and segMode[k-1] == 'REST':
            OCV0 = d['V'][sEnd[k-1]]
        else:
            continue
        Ibar = Ip.mean(); R0 = (Vp[0]-OCV0)/Ibar
        f = fitRC(tp, Vp, OCV0, Ibar, R0)
        rows.append(dict(Step=int(d['step'][a]), Dir='chg' if Ibar>0 else 'dis',
                         SOC=SOCstart[k], R0=R0, R1=f['R1'], tau1=f['tau1'],
                         R2=f['R2'], tau2=f['tau2']))
    Tall = pd.DataFrame(rows)
    Tall.to_csv(f"{RESULT}/verify_PerPulse_RAW_{ts['T']}C.csv", index=False)
    perPulse.append(Tall)
    nd = (Tall.Dir=='dis').sum(); nc = (Tall.Dir=='chg').sum()
    print(f"  T={ts['T']}C: {len(Tall)} pulses ({nd} dis, {nc} chg), "
          f"SOC {Tall.SOC.min():.3f}..{Tall.SOC.max():.3f}")
    print(f"    R0 dis: {Tall[Tall.Dir=='dis'].R0.min()*1e3:.3f}..{Tall[Tall.Dir=='dis'].R0.max()*1e3:.3f} mOhm | "
          f"R0 chg: {Tall[Tall.Dir=='chg'].R0.min()*1e3:.3f}..{Tall[Tall.Dir=='chg'].R0.max()*1e3:.3f} mOhm")

def interp_ex(x, xp, fp):
    """linear interpolation inside coverage, CONSTANT (nearest-edge) extension
    outside. Never slope-extrapolate noisy pulse fits: at 15C the first
    charge pulse was aborted, leaving a 15%-SOC hole where a linear slope
    extrapolation explodes the RC states."""
    o = np.argsort(xp); xp, fp = np.asarray(xp)[o], np.asarray(fp)[o]
    out = np.interp(x, xp, fp, left=np.nan, right=np.nan)
    out[x < xp[0]] = fp[0]; out[x > xp[-1]] = fp[-1]
    return out

PARAMS = ['R0','R1','tau1','R2','tau2']
# averaged 2D LUT
G_mean = {p: np.full((len(SOC_GRID), len(Tbp)), np.nan) for p in PARAMS}
G_dis  = {p: np.full((len(SOC_GRID), len(Tbp)), np.nan) for p in PARAMS}
G_chg  = {p: np.full((len(SOC_GRID), len(Tbp)), np.nan) for p in PARAMS}
for j, Tall in enumerate(perPulse):
    dis = Tall[Tall.Dir=='dis']; chg = Tall[Tall.Dir=='chg']
    for p in PARAMS:
        gd = interp_ex(SOC_GRID, dis.SOC.to_numpy(), dis[p].to_numpy())
        gc = interp_ex(SOC_GRID, chg.SOC.to_numpy(), chg[p].to_numpy())
        G_mean[p][:,j] = (gd+gc)/2; G_dis[p][:,j] = gd; G_chg[p][:,j] = gc

def addCaps(G):
    G['C1'] = G['tau1']/G['R1']; G['C2'] = G['tau2']/G['R2']
for G in (G_mean, G_dis, G_chg): addCaps(G)

def wideCSV(G, OCV, fname):
    T = pd.DataFrame({'SOCbp': SOC_GRID})
    for j, Tk in enumerate(Tbp):
        T[f'OCV_V_{Tk}C']  = OCV[:,j]
        for pname, mat in [('R0_ohm',G['R0']),('R1_ohm',G['R1']),('C1_F',G['C1']),
                           ('tau1_s',G['tau1']),('R2_ohm',G['R2']),('C2_F',G['C2']),('tau2_s',G['tau2'])]:
            T[f'{pname}_{Tk}C'] = mat[:,j]
    T.to_csv(f'{RESULT}/{fname}', index=False)
    return T

LUT_mean = wideCSV(G_mean, Vavg_all, 'verify_LUT2D_Averaged_batteryParams_2RC.csv')
LUT_dis  = wideCSV(G_dis , Vdis_all, 'verify_LUT2D_Discharge_2RC.csv')
LUT_chg  = wideCSV(G_chg , Vchg_all, 'verify_LUT2D_Charge_2RC.csv')

plt.figure(figsize=(13,4.6))
cmap = plt.cm.coolwarm
tcols = {Tk: cmap(i/max(len(perPulse)-1,1)) for i,Tk in enumerate(Tbp)}
for i, pn in enumerate(['R0','R1','R2']):
    plt.subplot(1,3,i+1)
    for j, pp in enumerate(perPulse):
        Tk = Tbp[j]; col = tcols[Tk]
        dd = pp[pp.Dir=='dis']; cc = pp[pp.Dir=='chg']
        plt.plot(dd.SOC*100, dd[pn]*1e3, 'o', color=col, ms=4, mfc='none',
                 label=f'{Tk}°C raw')
        plt.plot(cc.SOC*100, cc[pn]*1e3, 'x', color=col, ms=3, alpha=.6)
        plt.plot(SOC_GRID*100, G_mean[pn][:,j]*1e3, '-', color=col, lw=1.4,
                 label=f'{Tk}°C grid (avg)')
    plt.xlabel('SOC [%]'); plt.ylabel(f'{pn} [mOhm]'); plt.grid(alpha=.3)
    if i==0: plt.legend(fontsize=7)
plt.suptitle('Step 2 (2D): extracted pulse parameters vs SOC (all temperatures)')
plt.tight_layout(); plt.savefig(f'{RESULT}/verify_step2_params.png', dpi=140); plt.close()

# ----------------------------- STEP 3 --------------------------------------
from scipy.interpolate import RegularGridInterpolator

print('\n--- Step 3 (2D): validation against each temperature\'s HPPC ---')

def lut_factory(M):
    """griddedInterpolant({soc_grid,Tbp}, M, 'linear', clamped extrap) mirror"""
    if M.shape[1] == 1 or len(Tbp) == 1:
        return lambda s, Tq: np.interp(s, SOC_GRID, M[:, 0])
    F = RegularGridInterpolator((SOC_GRID, Tbp), M, method='linear',
                                bounds_error=False, fill_value=None)
    return lambda s, Tq: F([min(max(s,SOC_GRID[0]),SOC_GRID[-1]),
                            min(max(Tq,Tbp.min()),Tbp.max())]).item() if np.isscalar(s) else            F(np.c_[np.clip(s,SOC_GRID[0],SOC_GRID[-1]),
                   np.full_like(s, min(max(Tq,Tbp.min()),Tbp.max()), dtype=float)])

def simulate(d, Q, Tk, G, OCV, dual=None):
    t, V, I = d['t'], d['V'], d['I']
    dt = np.diff(t, prepend=0.0); dt[dt<0] = 0
    Ah = np.cumsum(I*dt)/3600.0
    SOC = np.clip(1-(Ah.max()-Ah)/Q, 0, 1)
    fOCV = lut_factory(OCV)
    fR0  = lut_factory(G['R0']);  fR1 = lut_factory(G['R1']); ft1 = lut_factory(G['tau1'])
    fR2  = lut_factory(G['R2']);  ft2 = lut_factory(G['tau2'])
    if dual is not None:
        Gg, OO = dual
        gOCV = lut_factory(OO);    gR0 = lut_factory(Gg['R0'])
        gR1  = lut_factory(Gg['R1']); gt1 = lut_factory(Gg['tau1'])
        gR2  = lut_factory(Gg['R2']); gt2 = lut_factory(Gg['tau2'])
    Vsim = np.zeros_like(t); Vc1 = np.zeros_like(t); Vc2 = np.zeros_like(t)
    Vsim[0] = fOCV(SOC[0],Tk) + I[0]*fR0(SOC[0],Tk)
    for k in range(1, len(t)):
        s = SOC[k]
        if dual is not None and I[k] > 0:
            Oc, R0k = gOCV(s,Tk), gR0(s,Tk)
            R1k, t1 = gR1(s,Tk), max(gt1(s,Tk),1e-3)
            R2k, t2 = gR2(s,Tk), max(gt2(s,Tk),1e-3)
        else:
            Oc, R0k = fOCV(s,Tk), fR0(s,Tk)
            R1k, t1 = fR1(s,Tk), max(ft1(s,Tk),1e-3)
            R2k, t2 = fR2(s,Tk), max(ft2(s,Tk),1e-3)
        dec1 = np.exp(-dt[k]/t1); Vc1[k] = Vc1[k-1]*dec1 + R1k*I[k]*(1-dec1)
        dec2 = np.exp(-dt[k]/t2); Vc2[k] = Vc2[k-1]*dec2 + R2k*I[k]*(1-dec2)
        Vsim[k] = Oc + I[k]*R0k + Vc1[k] + Vc2[k]
    return t, V, Vsim

nT = len(TESTS)
fig, ax = plt.subplots(2, nT, figsize=(7.2*nT, 6.2), squeeze=False)
summary = dict(T=Tbp.tolist(), Q_Ah=np.round(Qs,4).tolist(), temps=[], rmse_avg_mV=[],
               rmse_dual_mV=[], rmse_dual_chg_mV=[], rmse_dual_dis_mV=[])
for j, ts in enumerate(TESTS):
    Tk = ts['T']; d = load(ts['hppc'])
    t, V, Vs_m = simulate(d, Qs[j], Tk, G_mean, Vavg_all, dual=None)
    _, _, Vs_d = simulate(d, Qs[j], Tk, G_dis, Vdis_all, dual=(G_chg, Vchg_all))
    I = d['I']
    rm_m  = np.sqrt(np.mean((Vs_m-V)**2))*1e3
    rm_d  = np.sqrt(np.mean((Vs_d-V)**2))*1e3
    rm_c  = np.sqrt(np.mean((Vs_d[I>0.5]-V[I>0.5])**2))*1e3
    rm_dd = np.sqrt(np.mean((Vs_d[I<-0.5]-V[I<-0.5])**2))*1e3
    print(f'  T={Tk}C: averaged RMSE = {rm_m:.2f} mV | dual RMSE = {rm_d:.2f} mV '
          f'(chg {rm_c:.2f} / dis {rm_dd:.2f}), max|e| dual = {np.max(np.abs(Vs_d-V))*1e3:.1f} mV')
    summary['temps'].append(Tk); summary['rmse_avg_mV'].append(round(rm_m,2))
    summary['rmse_dual_mV'].append(round(rm_d,2))
    summary['rmse_dual_chg_mV'].append(round(rm_c,2)); summary['rmse_dual_dis_mV'].append(round(rm_dd,2))
    ax[0][j].plot(t/3600, V, 'k', lw=1.0, label='Measured')
    ax[0][j].plot(t/3600, Vs_m, 'r--', lw=.7, label=f'Averaged ({rm_m:.1f} mV)')
    ax[0][j].plot(t/3600, Vs_d, 'b:',  lw=.8, label=f'Dual ({rm_d:.1f} mV)')
    ax[0][j].set_title(f'HPPC @ {Tk}\u00b0C — 2-D LUT validation')
    ax[0][j].set_ylabel('V [V]'); ax[0][j].legend(loc='lower left', fontsize=8)
    ax[1][j].plot(t/3600, (Vs_m-V)*1e3, 'r', lw=.4, label='avg err')
    ax[1][j].plot(t/3600, (Vs_d-V)*1e3, 'b', lw=.4, label='dual err')
    ax[1][j].set_ylabel('Error [mV]'); ax[1][j].set_xlabel('Time [h]'); ax[1][j].legend(fontsize=8)
    for a in (ax[0][j], ax[1][j]): a.grid(alpha=.3)
plt.tight_layout()
plt.savefig(f'{RESULT}/verify_step3_validation_ALL_T.png', dpi=140); plt.close()

with open(f'{RESULT}/verify_summary.json','w') as f: json.dump(summary, f, indent=2)
tag = '-'.join(str(t['T']) for t in TESTS)
import shutil
shutil.copy(f'{RESULT}/verify_step3_validation_ALL_T.png', f'{RESULT}/verify{tag}C_step3_validation.png')
print('\nWrote results/verify_step3_validation_ALL_T.png and verify_summary.json')
print(json.dumps(summary, indent=2))
