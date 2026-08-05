%% step5_BuildPackParams.m
%  Scale the single-cell Simscape tables to a 192S4P pack represented by
%  ONE Battery Equivalent Circuit block.
%
%  The block mask uses the R + tau (resistance + time-constant) form, so
%  there are NO capacitance (C1/C2) tables to scale -- only R, tau, OCV,
%  capacity and hysteresis voltage. (tau is invariant under pack scaling
%  because tau = R*C: R scales x Ns/Np while C scales x Np/Ns.)
%
%  Scaling rules (identical, balanced cells, one equivalent block):
%    series  (Ns) : voltage-like quantities scale x Ns
%    parallel(Np) : capacity scales x Np
%    R            : scales x Ns/Np
%    tau = R*C    : invariant, scales x 1
%    hysteresis V : scales x Ns
%
%  Run AFTER step4_Build_SimscapeModel.m (needs Simscape_BatteryParams.mat
%  and LUT2D_Charge_2RC.mat).

clear; clc
[~, ~, resultsDir] = testConfig();

Ns = 192;   % series cells
Np = 4;     % parallel strings

scaleV   = Ns;        % OCV, hysteresis voltage  -> x 192
scaleR   = Ns/Np;     % resistances              -> x 48
scaleTau = 1;         % tau = R*C, invariant     -> x 1
scaleCap = Np;        % capacity                 -> x 4

%% --- single-cell Simscape tables (R + tau form) -----------------------
S  = load(fullfile(resultsDir, 'Simscape_BatteryParams.mat'));
Lc = load(fullfile(resultsDir, 'LUT2D_Charge_2RC.mat'));

SOC_bp = S.SOC_bp;
T_bp   = S.T_bp;                       % Kelvin (block uses K)

% OCV and resistance tables (SOC x T)
OCV_tab   = S.OCV_tab;
R0_chg = S.R0_chg;  R0_dis = S.R0_dis;
R1_chg = S.R1_chg;  R1_dis = S.R1_dis;
R2_chg = S.R2_chg;  R2_dis = S.R2_dis;

% RC time constants (SOC x T) -- used directly by the mask, no C tables
tau1_chg = S.tau1_chg; tau1_dis = S.tau1_dis;
tau2_chg = S.tau2_chg; tau2_dis = S.tau2_dis;

% hysteresis (voltage) tables (SOC x T) -- scale like OCV (x Ns)
MaxHyst_tab  = S.MaxHyst_tab;
InstHyst_tab = S.InstHyst_tab;
HystRate_tab = S.HystRate_tab;         % dimensionless rate, x 1

% per-temperature full-cell capacity vector -> pack capacity vector
Q_Ah   = Lc.Q_Ah;
Q_Ah_p = Q_Ah * scaleCap;

%% --- pack-scaled tables ----------------------------------------------
OCV_tab_p   = OCV_tab * scaleV;

R0_chg_p = R0_chg * scaleR;  R0_dis_p = R0_dis * scaleR;
R1_chg_p = R1_chg * scaleR;  R1_dis_p = R1_dis * scaleR;
R2_chg_p = R2_chg * scaleR;  R2_dis_p = R2_dis * scaleR;

tau1_chg_p = tau1_chg * scaleTau;  tau1_dis_p = tau1_dis * scaleTau;
tau2_chg_p = tau2_chg * scaleTau;  tau2_dis_p = tau2_dis * scaleTau;

MaxHyst_tab_p  = MaxHyst_tab  * scaleV;
InstHyst_tab_p = InstHyst_tab * scaleV;
HystRate_tab_p = HystRate_tab * 1;

% representative pack capacity (25 C entry) + full per-T vector
Cap_Ah_p = Q_Ah_p(2);                  % ~140 Ah
Q_pack   = Cap_Ah_p;

save(fullfile(resultsDir, 'Simscape_BatteryPackParams.mat'), ...
     'SOC_bp','T_bp', ...
     'OCV_tab_p','MaxHyst_tab_p','InstHyst_tab_p','HystRate_tab_p', ...
     'R0_chg_p','R0_dis_p','R1_chg_p','R1_dis_p','tau1_chg_p','tau1_dis_p', ...
     'R2_chg_p','R2_dis_p','tau2_chg_p','tau2_dis_p', ...
     'Q_Ah_p','Cap_Ah_p','Q_pack','Ns','Np');

fprintf('Saved Simscape_BatteryPackParams.mat\n');
fprintf('  Pack: %dS%dP   OCV~%.1f V   R x%.0f   tau x%.0f   Cap~%.1f Ah (%.1f kWh)\n', ...
        Ns, Np, OCV_tab_p(end,2), scaleR, scaleTau, Cap_Ah_p, ...
        OCV_tab_p(end,2)*Cap_Ah_p/1000);
