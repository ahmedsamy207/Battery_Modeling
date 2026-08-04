%% step4_Build_SimscapeModel.m
%  Stage 4 (data prep ONLY). The Simscape model itself is built by hand —
%  see step4_ManualBuild.md for the exact blocks, wiring and mask variables.
%  This script just loads the 2-D LUTs and creates the base-workspace
%  variables that the 'Battery Equivalent Circuit' block mask will reference,
%  and saves them to results/Simscape_BatteryParams.mat. Pure MATLAB — no
%  Simulink / add_block calls — so it always runs.
%
%  Run AFTER step2 (needs LUT2D_Charge/Discharge/Averaged_2RC.mat).
%  Then build the model manually and run step4_Validate.m.

clear; clc
[~, ~, resultsDir] = testConfig();

USE_HYSTERESIS = true;     % recover the measured OCV hysteresis loop
HYST_RATE      = 1.0;      % one-state hysteresis rate (tunable; see .md)

%% load DUAL LUTs (charge & discharge R/RC) + AVERAGED OCV
Lc = load(fullfile(resultsDir, 'LUT2D_Charge_2RC.mat'));
Ld = load(fullfile(resultsDir, 'LUT2D_Discharge_2RC.mat'));
La = load(fullfile(resultsDir, 'LUT2D_Averaged_2RC.mat'));

soc_grid = Lc.soc_grid;                  % 201 x 1
Tbp      = Lc.Tbp(:).' + 273.15;         % 1 x nT, in KELVIN (block uses K)
Cap_Ah   = mean(Lc.Q_Ah);               % representative capacity [Ah]
T_sim    = Tbp(1);                       % default operating temperature [K]

% charge / discharge resistance + RC time-constant tables (SOC x T)
R0_chg = Lc.R0_ohm;   R0_dis = Ld.R0_ohm;
R1_chg = Lc.R1_ohm;   R1_dis = Ld.R1_ohm;
tau1_chg = Lc.tau1_s; tau1_dis = Ld.tau1_s;
R2_chg = Lc.R2_ohm;   R2_dis = Ld.R2_ohm;
tau2_chg = Lc.tau2_s; tau2_dis = Ld.tau2_s;

% single OCV = thermodynamic center; hysteresis = half the branch gap
OCV_base = La.OCV_V;                                  % averaged (center) OCV
MaxHyst_tab  = abs(0.5*(Lc.OCV_V - Ld.OCV_V));       % half-loop amplitude [V] (must be >= 0)
InstHyst_tab = zeros(size(MaxHyst_tab));              % no instantaneous term
HystRate_tab = HYST_RATE * ones(size(MaxHyst_tab));    % (SOC,T), constant rate

% expose to BASE workspace (the block mask references these by name)
SOC_bp = soc_grid;  T_bp = Tbp;  Cap_Ah = Cap_Ah;  T_sim = T_sim;
OCV_tab = OCV_base;  MaxHyst_tab = MaxHyst_tab;  InstHyst_tab = InstHyst_tab;  HystRate_tab = HystRate_tab;
R0_chg = R0_chg; R0_dis = R0_dis;
R1_chg = R1_chg; R1_dis = R1_dis; tau1_chg = tau1_chg; tau1_dis = tau1_dis;
R2_chg = R2_chg; R2_dis = R2_dis; tau2_chg = tau2_chg; tau2_dis = tau2_dis;

% save a standalone copy (used by step4_Validate.m and for manual build)
save(fullfile(resultsDir, 'Simscape_BatteryParams.mat'), ...
     'SOC_bp','T_bp','Cap_Ah','T_sim','OCV_tab','MaxHyst_tab','InstHyst_tab','HystRate_tab', ...
     'R0_chg','R0_dis','R1_chg','R1_dis','tau1_chg','tau1_dis', ...
     'R2_chg','R2_dis','tau2_chg','tau2_dis');

fprintf('Created base-workspace variables and saved Simscape_BatteryParams.mat\n');
fprintf('Variables: SOC_bp, T_bp, Cap_Ah, T_sim, OCV_tab, MaxHyst_tab,\n');
fprintf('          InstHyst_tab, HystRate_tab, R0/R1/tau1/R2/tau2 _chg / _dis\n\n');
fprintf('Now build the model by hand — see step4_ManualBuild.md.\n');
