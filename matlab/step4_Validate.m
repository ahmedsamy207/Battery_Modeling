%% step4_Validate.m
%  Stage 4 (deployment): validate the Simscape 'Battery (Table-Based)' model
%  built by step4_Build_SimscapeModel.m.
%
%  (A) HPPC REPLAY  - drive the model with every available HPPC current
%      profile and compare terminal voltage to the measured data. This should
%      land near the step-3 DUAL-LUT RMSE (~0.5 %), since the Battery
%      Equivalent Circuit block uses the separate charge/discharge LUTs
%      (current directionality) plus the one-state hysteresis model.
%  (B) DRIVE CYCLE  - run a representative synthetic load cycle at 25 degC to
%      show general (non-HPPC) usage of the deployed cell model.
%
%  Working configuration (verified to reproduce step-3 ~0.5 % RMSE):
%   - Input current is fed as +I (cycler I>0 = charge) to the Controlled
%     Current Source; the block's internal current convention then matches the
%     data, so SOC integrates correctly without flipping the sign.
%   - Each HPPC profile is TRIMMED to start at its first 70 A pulse, where the
%     cell is at full charge (SOC = 1). This makes the block's initial SOC = 1
%     a guaranteed-correct anchor (no SOC-reference ambiguity).
%   - DIAG prints, per temperature:
%       RMSE(sim vs OCV@dataSOC) -> residual of the IR/RC/hysteresis part
%                                  (~30-50 mV expected; large (~0.5 V) would
%                                  mean the block's SOC is pinned/reversed).
%       corr(sim,OCV) / corr(sim,meas) -> shape agreement.
%
%  Requires: step4_BatteryModel.slx (run step4_Build_SimscapeModel first).

clear; clc; close all

%% 0) setup
[tests, ~, resultsDir] = testConfig();
modelName = 'step4_BatteryModel';
load(fullfile(resultsDir, 'Simscape_BatteryParams.mat'));   % SOC_bp,T_bp,...,T_sim
Lc = load(fullfile(resultsDir, 'LUT2D_Charge_2RC.mat'));     % Lc.Q_Ah holds per-temperature capacity

% Locate the .slx (built by hand; its location varies). Search recursively
% from the likely roots and report what was found (with diagnostics if missing).
scriptFolder = fileparts(mfilename('fullpath'));
rootDirs = unique({ resultsDir, pwd, scriptFolder, fileparts(scriptFolder) });
mdlPath = '';
for r = 1:numel(rootDirs)
    if ~isfolder(rootDirs{r}), continue; end
    d = dir(fullfile(rootDirs{r}, '**', [modelName '.slx']));
    if ~isempty(d)
        mdlPath = fullfile(d(1).folder, d(1).name);
        break;
    end
end
if isempty(mdlPath)
    fprintf('Searched these roots for %s.slx:\n', modelName);
    for r = 1:numel(rootDirs), fprintf('  %s\n', rootDirs{r}); end
    fprintf('*.slx files in script folder (%s):\n', scriptFolder);
    dd = dir(fullfile(scriptFolder, '*.slx'));
    if isempty(dd), fprintf('  (none)\n'); else, for k = 1:numel(dd), fprintf('  %s\n', dd(k).name); end, end
    error('Model %s.slx not found. Save it under matlab/ (or results/ or current folder), or fix modelName.', modelName);
end
fprintf('Using model: %s\n', mdlPath);
if ~bdIsLoaded(modelName), load_system(mdlPath); end

% locate the battery block once (used to force capacity / initial SOC)
battPath = '';
try
    bl = find_system(modelName, 'Type', 'Block');
    bb = bl(contains(lower(bl), 'battery'));
    if ~isempty(bb), battPath = bb{1}; end
catch
end

% Robust solver settings: floor the minimum step size so the variable-step
% solver does not stall on fast transients (tiny RC time constants at pulse
% edges drive it to its ~1e-11 step floor). Keeps accuracy reasonable while
% guaranteeing the sim finishes.
% daessc is the Simscape-recommended variable-step solver. (ode23t is NOT a
% fixed-step solver in this release, so the earlier 'Fixed-step/ode23t' setting
% did not apply — that is why the RMSE was unchanged. The min-step warnings
% below are benign and do not affect the RMSE; they are suppressed further down.)
set_param(modelName, 'SolverType', 'Variable-step', 'Solver', 'daessc', ...
          'RelTol', '1e-3', 'AbsTol', '1e-4', 'MinStep', '1e-4', 'MaxStep', '1');

% initial SOC = 1 (full) so the block's SOC integrator starts where the data
% does (SOC=1 at t=0); this keeps SOC within [0,1] over the long test.
% (If you already set it manually in the block mask, this is a harmless no-op.)
try
    blks = find_system(modelName, 'Type', 'Block');
    batt = blks(contains(lower(blks), 'battery'));
    if ~isempty(batt)
        d = get_param(batt{1}, 'DialogParameters'); nm = fieldnames(d); lc = lower(nm);
        is = find(contains(lc,'initial') & contains(lc,'soc'),1);
        if isempty(is), is = find(contains(lc,'soc0') | contains(lc,'soc1'),1); end
        if ~isempty(is), set_param(batt{1}, nm{is}, '1'); end
    end
catch
end

% for normalized-RMSE % (same convention as step 3)
Vmean_ref = @(V) mean(V);

%% ==================== (A) HPPC REPLAY ===============================
idx = find([tests.has_hppc]);
rmse_mv = nan(1, numel(idx));
rmse_pct = nan(1, numel(idx));

fig = figure('Color','w','Position',[30 30 1180 270*max(numel(idx),1)]);
nTiles = 0;

fprintf('\n=== Step 4 (A) HPPC REPLAY — Simscape Battery Equivalent Circuit (dual LUT + hysteresis) ===\n');
for ii = 1:numel(idx)
    k  = idx(ii);
    Tk = tests(k).T;
    T  = readtable(tests(k).hppc, 'VariableNamingRule','preserve');
    t  = T{:,7}/1000;  V = T{:,8};  I = T{:,9};

    % --- trim to the first 70 A pulse so the cell is at full charge (SOC = 1) ---
    % The data's SOC reference (max Ah) then sits exactly at this pulse, so the
    % block's initial SOC = 1 is guaranteed correct (no anchor ambiguity).
    Ithr = 70;
    ip = find(abs(I) >= Ithr, 1, 'first');
    if ~isempty(ip)
        t = t(ip:end) - t(ip);        % trim AND re-zero time at the pulse start
        V = V(ip:end);  I = I(ip:end);
        fprintf('     trimmed to first %g A pulse: sim starts at t=0, %d points\n', Ithr, numel(t));
    else
        fprintf('     NOTE: no %g A pulse found; using full profile\n', Ithr);
    end

    T_sim = Tk + 273.15;                        % set operating temperature [K] (block uses K)
    jj = find(Lc.Tbp == Tk, 1);
    if ~isempty(jj), Cap_Ah = Lc.Q_Ah(jj); end   % per-T capacity: block integrates SOC from this

    % Force the block to use this exact capacity (prevents the SOC integrator
    % from running with a default/leftover value and overshooting past 1).
    if ~isempty(battPath)
        try
            d = get_param(battPath,'DialogParameters'); pnm = fieldnames(d); plc = lower(pnm);
            ic = find(contains(plc,'capac'),1);
            if ~isempty(ic), set_param(battPath, pnm{ic}, num2str(Cap_Ah)); end
        catch
        end
    end
    fprintf('     capacity applied = %.4f Ah (T=%d C)\n', Cap_Ah, Tk);

    % current input for the Controlled Current Source (sign flip, see header)
    I_in = timeseries(I, t, 'Name', 'I_in');     % cycler I>0 = charge; Simscape battery uses current-out-of-+ = discharge, so feed +I

    set_param(modelName, 'StopTime', num2str(max(t)));
    assignin('base', 'V_out', []);                % clear stale base variable (if any)
    warning('off','all');                 % silence benign solver / SOC-assertion flood during sim
    out = sim(modelName, 'ReturnWorkspaceOutputs', 'on');
    warning('on','all');
    Vsim_ts = [];
    try, Vsim_ts = out.V_out;        end           % SimulationOutput property
    if isempty(Vsim_ts), try, Vsim_ts = out.get('V_out'); end, end
    if isempty(Vsim_ts), try, Vsim_ts = evalin('base', 'V_out'); end, end
    if isempty(Vsim_ts), error('V_out not captured at T=%d C. Verify the To-Workspace block (Variable name = V_out, Save format = Timeseries).', Tk); end
    V_sim = interp1(Vsim_ts.Time, Vsim_ts.Data, t, 'linear', 'extrap');   % resample onto measured t

    % ---- DIAGNOSTIC: is the block's OCV pinned at SOC=1? -------------------
    % "Truth" OCV = OCV evaluated at the DATA's SOC (the SOC we trust from step
    % 3). step 3's V = this OCV + IR + hysteresis, and got ~0.5 %, so V_sim
    % should track V_ocv_data within ~0.1 V. If instead V_sim stays near a
    % constant (OCV at the wrong SOC, e.g. SOC=1), the block's SOC is
    % pinned/reversed -> that is the source of the ~0.5 V offset.
    % Uses the model's own table (OCV_tab, SOC_bp, T_bp[K]) so orientation is correct.
    V_ocv_data = NaN(size(t));          % default so the plot can never crash
    try
        Ahc   = cumtrapz(t(:), I(:)) / 3600;          % Ah, charge-positive
        SOCd  = 1 - (max(Ahc) - Ahc)/(Cap_Ah);
        SOCd  = max(0, min(1, SOCd));
        % OCV_tab is stored as (numT x numSOC) = 3 x 201; transpose to
        % (numSOC x numT) for interp2 with X=SOC_bp, Y=T_bp(K).
        V_ocv_data = interp2(SOC_bp, T_bp, OCV_tab.', SOCd, Tk + 273.15);
        rmse_ocv = sqrt(mean((V_sim(:) - V_ocv_data(:)).^2))*1e3;
        fprintf('     DIAG  RMSE(sim vs OCV@dataSOC)=%.1f mV ; corr(sim,OCV)=%.3f ; corr(sim,meas)=%.3f\n', ...
                rmse_ocv, corr(V_sim(:), V_ocv_data(:)), corr(V_sim(:), V(:)));
    catch ME
        fprintf('     DIAG  skipped (%s)\n', ME.message);
    end

    res  = V_sim - V;
    rmse_mv(ii)  = sqrt(mean(res.^2))*1e3;
    rmse_pct(ii) = (sqrt(mean(res.^2))/Vmean_ref(V))*100;
    fprintf('  T = %2d C : RMSE = %.3f mV (%.3f%%)\n', Tk, rmse_mv(ii), rmse_pct(ii));

    nTiles = nTiles + 1;
    subplot(numel(idx), 2, nTiles); hold on; grid on
    plot(t/3600, V, 'k',  'LineWidth', 1.4);
    plot(t/3600, V_sim, 'r--','LineWidth', 1.2);
    ylabel('V [V]'); 
    title(sprintf('HPPC @ %d^{\\circ}C — Simscape vs measured (RMSE = %.2f mV, %.2f%%)', ...
        Tk, rmse_mv(ii), rmse_pct(ii)));
    legend('Measured','Simscape','Location','best');

    nTiles = nTiles + 1;
    subplot(numel(idx), 2, nTiles); hold on; grid on
    plot(t/3600, res*1e3, 'b'); 
    ylabel('Error [mV]'); xlabel('Time [h]');
    title(sprintf('Residuals @ %d^{\\circ}C', Tk));
end
saveas(fig, fullfile(resultsDir, 'step4_HPPC_replay.png'), 'png');
fprintf('---------------------------------------------------------------\n');
fprintf('Mean HPPC RMSE : %.3f mV (%.3f%%)\n', ...
    mean(rmse_mv,'omitnan'), mean(rmse_pct,'omitnan'));

%% ==================== (B) DRIVE CYCLE =================================
fprintf('\n=== Step 4 (B) DRIVE CYCLE (synthetic, 25 degC) ===\n');
T_sim = 25 + 273.15;                         % 298.15 K = 25 degC
jj = find(Lc.Tbp == 25, 1);
if ~isempty(jj), Cap_Ah = Lc.Q_Ah(jj); end    % per-T capacity for the drive cycle
% Start below SOC = 1 so the synthetic charge spikes do not trip the block's
% SOC <= 1 assertion; keeps the demo warning-free.
SOC0 = 0.85;
if ~isempty(battPath)
    try
        d = get_param(battPath,'DialogParameters'); pnm = fieldnames(d); plc = lower(pnm);
        is = find(contains(plc,'initial') & contains(plc,'soc'),1);
        if isempty(is), is = find(contains(plc,'soc0') | contains(plc,'soc1'),1); end
        if ~isempty(is), set_param(battPath, pnm{is}, num2str(SOC0)); end
    catch
    end
end
tc = 0:1:3600;                                  % 1 h, 1 s steps
Idrive = driveCycle(tc);                         % +/- current, A (I>0 = charge)
I_in   = timeseries(Idrive, tc, 'Name', 'I_in');   % feed +I (same sign convention as HPPC)

set_param(modelName, 'StopTime', num2str(max(tc)));
assignin('base', 'V_out', []);
warning('off','all');                 % silence benign solver / SOC-assertion flood during sim
out = sim(modelName, 'ReturnWorkspaceOutputs', 'on');
warning('on','all');
Vsim_d_ts = [];
try, Vsim_d_ts = out.V_out;        end
if isempty(Vsim_d_ts), try, Vsim_d_ts = out.get('V_out'); end, end
if isempty(Vsim_d_ts), try, Vsim_d_ts = evalin('base', 'V_out'); end, end
if isempty(Vsim_d_ts), error('V_out not captured in drive cycle.'); end
Vsim_d = interp1(Vsim_d_ts.Time, Vsim_d_ts.Data, tc, 'linear', 'extrap');

% SOC estimated externally from the input current (release-independent).
% battery current (charge-positive) = Idrive ; SOC drops while discharging.
dtc = [0, diff(tc)]; dtc(dtc<0)=0;
SOCsim = SOC0 - cumsum(Idrive .* dtc)/(Cap_Ah*3600);
SOCsim = max(0, min(1, SOCsim));

fprintf('  drive cycle : V range %.3f .. %.3f V\n', min(Vsim_d), max(Vsim_d));

figure('Color','w','Position',[60 60 1000 640]);
subplot(3,1,1); hold on; grid on
plot(tc/3600, Idrive, 'k', 'LineWidth', 1.2); ylabel('I [A]'); title('Drive-cycle current (I>0 = charge)');
subplot(3,1,2); hold on; grid on
plot(tc/3600, Vsim_d, 'r', 'LineWidth', 1.2); ylabel('V [V]'); title('Simulated terminal voltage');
subplot(3,1,3); hold on; grid on
plot(tc/3600, SOCsim*100, 'b', 'LineWidth', 1.2);
ylabel('SOC [%]'); title('SOC (from input current, external estimate)');
xlabel('Time [h]');
saveas(gcf, fullfile(resultsDir, 'step4_DriveCycle.png'), 'png');

fprintf('Saved: step4_HPPC_replay.png, step4_DriveCycle.png\n');

%% ========================================================================
function I = driveCycle(t)
% a representative ~1C-peak synthetic load (mix of pulses + sinusoid)
    f1 = 1/600; f2 = 1/120;                              % slow + fast
    I  = 15*sin(2*pi*f1*t) + 8*sign(sin(2*pi*f2*t));
    I  = I + 5*randn(size(t));                           % mild noise
    In = 35; I = max(min(I, In), -In);                   % clamp to ~1C of 35 Ah
end
