%% step5_Validate.m
%  Stage 5 (deployment): validate the ONE-block 192S4P pack model
%  (batteryPack.slx) by replaying the SAME trimmed HPPC profiles used in
%  step4, but at the pack level.
%
%  Scaling (identical, balanced cells -> one equivalent block):
%    - Feed the PACK terminal current  I_in = I_cell * Np.  Each of the Np
%      parallel strings then carries the original cell HPPC current (same
%      C-rate per cell). With Cap_pack = Cap_cell * Np the SOC integrator
%      evolves identically, so the SOC=1 anchor at the first 70 A pulse holds.
%    - Compare the pack terminal voltage to the measured cell voltage scaled
%      by Ns:  V_meas_pack = V_cell * Ns.  (The block's OCV is already x Ns.)
%
%  Expect the SAME ~0.9 % RMSE as step4 (absolute error scales x Ns, % does
%  not). The DIAG residual (sim vs OCV@dataSOC) scales x Ns too (~7 V).
%
%  Requires: batteryPack.slx (built by hand, referencing the *_p vars from
%            Simscape_BatteryPackParams.mat) and step5_BuildPackParams.m run.

clear; clc; close all

%% 0) setup
[tests, ~, resultsDir] = testConfig();
modelName = 'step4_BatteryModel';
load(fullfile(resultsDir, 'Simscape_BatteryPackParams.mat'));   % *_p tables, Q_Ah_p, Ns, Np

% --- pack scaling factors (from the saved .mat) -------------------------
Ns = Ns; Np = Np;                       % 192 , 4  (kept explicit for clarity)
scaleV   = Ns;                          % measured voltage & OCV reference
scaleI   = Np;                          % pack terminal current vs cell current
%   To apply the SAME ABSOLUTE current at the pack terminals instead (each
%   cell then sees I_cell/Np, a gentler test), set scaleI = 1.

% Locate the .slx (built by hand; its location varies). Search recursively.
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

% Same robust solver as step4 (daessc, variable-step; min-step warnings benign)
set_param(modelName, 'SolverType', 'Variable-step', 'Solver', 'daessc', ...
          'RelTol', '1e-3', 'AbsTol', '1e-4', 'MinStep', '1e-4', 'MaxStep', '1');

% initial SOC = 1 (full) so the block's SOC integrator starts where the data
% does (SOC=1 at the first 70 A pulse); keeps SOC within [0,1] over the test.
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

Vmean_ref = @(V) mean(V);

%% ==================== HPPC REPLAY (pack level) ========================
idx = find([tests.has_hppc]);
rmse_V    = nan(1, numel(idx));
rmse_pct  = nan(1, numel(idx));

fig = figure('Color','w','Position',[30 30 1180 270*max(numel(idx),1)]);
nTiles = 0;

fprintf('\n=== Step 5 HPPC REPLAY — 192S4P pack (one block), dual LUT + hysteresis ===\n');
fprintf('   scaleI = %g (pack terminal current = cell current x Np) ; scaleV = %g (meas V x Ns)\n', scaleI, scaleV);
for ii = 1:numel(idx)
    k  = idx(ii);
    Tk = tests(k).T;
    T  = readtable(tests(k).hppc, 'VariableNamingRule','preserve');
    t  = T{:,7}/1000;  V = T{:,8};  I = T{:,9};     % cell-level: t[s], V[V], I[A] (I>0 = charge)

    % --- trim to the first 70 A pulse so the cell is at full charge (SOC = 1) ---
    Ithr = 70;
    ip = find(abs(I) >= Ithr, 1, 'first');
    if ~isempty(ip)
        t = t(ip:end) - t(ip);        % trim AND re-zero time at the pulse start
        V = V(ip:end);  I = I(ip:end);
        fprintf('     trimmed to first %g A pulse: sim starts at t=0, %d points\n', Ithr, numel(t));
    else
        fprintf('     NOTE: no %g A pulse found; using full profile\n', Ithr);
    end

    T_sim = Tk + 273.15;                        % operating temperature [K] (block uses K)
    % per-T pack capacity: Cap_pack = Cap_cell * Np  -> pick the right entry
    jj = find(abs(T_bp - T_sim) < 0.5, 1);
    if ~isempty(jj), Cap_Ah = Q_Ah_p(jj); end   % pack capacity at this T

    % Force the block to use this exact (per-T) pack capacity.
    if ~isempty(battPath)
        try
            d = get_param(battPath,'DialogParameters'); pnm = fieldnames(d); plc = lower(pnm);
            ic = find(contains(plc,'capac'),1);
            if ~isempty(ic), set_param(battPath, pnm{ic}, num2str(Cap_Ah)); end
        catch
        end
    end
    fprintf('     pack capacity applied = %.4f Ah (T=%d C)\n', Cap_Ah, Tk);

    % ---- pack terminal current = cell current x Np (same C-rate per cell) ----
    I_pack = I * scaleI;
    I_in   = timeseries(I_pack, t, 'Name', 'I_in');   % +I, cycler I>0 = charge

    set_param(modelName, 'StopTime', num2str(max(t)));
    assignin('base', 'V_out', []);
    warning('off','all');                 % silence benign solver / SOC-assertion flood
    out = sim(modelName, 'ReturnWorkspaceOutputs', 'on');
    warning('on','all');
    Vsim_ts = [];
    try, Vsim_ts = out.V_out;        end
    if isempty(Vsim_ts), try, Vsim_ts = out.get('V_out'); end, end
    if isempty(Vsim_ts), try, Vsim_ts = evalin('base', 'V_out'); end, end
    if isempty(Vsim_ts), error('V_out not captured at T=%d C. Verify the To-Workspace block (Variable name = V_out, Save format = Timeseries).', Tk); end
    V_sim = interp1(Vsim_ts.Time, Vsim_ts.Data, t, 'linear', 'extrap');   % resample onto measured t

    % ---- DIAGNOSTIC: is the block's OCV pinned at SOC=1? -------------------
    % Truth OCV at the DATA's SOC, using the PACK table (OCV_tab_p, x Ns).
    % Residual ~ Ns * (step4 cell residual ~38 mV) ~ 7 V is expected; a much
    % larger / flat residual would mean SOC pinned/reversed.
    V_ocv_data = NaN(size(t));
    try
        Ahc   = cumtrapz(t(:), I_pack(:)) / 3600;          % Ah, charge-positive (pack)
        SOCd  = 1 - (max(Ahc) - Ahc)/(Cap_Ah);
        SOCd  = max(0, min(1, SOCd));
        % OCV_tab_p is (SOC x T) = 201 x 3; transpose to (numT x numSOC).
        V_ocv_data = interp2(SOC_bp, T_bp, OCV_tab_p.', SOCd, Tk + 273.15);
        rmse_ocv = sqrt(mean((V_sim(:) - V_ocv_data(:)).^2));
        fprintf('     DIAG  RMSE(sim vs OCV@dataSOC)=%.3f V ; corr(sim,OCV)=%.3f ; corr(sim,meas)=%.3f\n', ...
                rmse_ocv, corr(V_sim(:), V_ocv_data(:)), corr(V_sim(:), V(:)*scaleV));
    catch ME
        fprintf('     DIAG  skipped (%s)\n', ME.message);
    end

    % ---- compare pack sim to measured cell voltage scaled by Ns ----
    V_meas_pack = V * scaleV;
    res  = V_sim - V_meas_pack;
    rmse_V(ii)   = sqrt(mean(res.^2));
    rmse_pct(ii) = (sqrt(mean(res.^2))/Vmean_ref(V_meas_pack))*100;
    fprintf('  T = %2d C : RMSE = %.3f V (%.3f%%)   [= %.1f x cell RMSE]\n', ...
            Tk, rmse_V(ii), rmse_pct(ii), rmse_V(ii)/0.032);

    nTiles = nTiles + 1;
    subplot(numel(idx), 2, nTiles); hold on; grid on
    plot(t/3600, V_meas_pack, 'k',  'LineWidth', 1.4);   % measured cell V x Ns
    plot(t/3600, V_sim,       'r--','LineWidth', 1.2);   % Simscape block output
    ylabel('V [V]');
    title(sprintf('HPPC @ %d^{\\circ}C — Pack (RMSE = %.2f V, %.2f%%)', ...
        Tk, rmse_V(ii), rmse_pct(ii)));
    legend('Measured xNs','Simscape','Location','best');

    nTiles = nTiles + 1;
    subplot(numel(idx), 2, nTiles); hold on; grid on
    plot(t/3600, res, 'b');
    ylabel('Error [V]'); xlabel('Time [h]');
    title(sprintf('Residuals @ %d^{\\circ}C', Tk));
end
saveas(fig, fullfile(resultsDir, 'step5_HPPC_replay.png'), 'png');
fprintf('---------------------------------------------------------------\n');
fprintf('Mean HPPC RMSE : %.3f V (%.3f%%)\n', ...
    mean(rmse_V,'omitnan'), mean(rmse_pct,'omitnan'));
fprintf('  (expect ~0.9%% — identical to the cell result, since only the scale changed)\n');
