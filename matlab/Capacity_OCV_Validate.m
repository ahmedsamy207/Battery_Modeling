%% Capacity_OCV_Validate.m
%  Cell-level validation (like step4_Validate) of the single-cell Battery
%  Equivalent Circuit block (step4_BatteryModel.slx) against a CAPACITY or
%  OCV (GITT) characterization test, instead of HPPC.
%
%  NO pack scaling here (this is the CELL): feed I_cell directly, compare to
%  measured cell V directly. Capacity is checked against Lc.Q_Ah (cell, ~35 Ah).
%
%  Set testType below to 'capacity' or 'ocv'. (For HPPC replay use step4_Validate.)

clear; clc; close all
testType = 'ocv';   % 'capacity' | 'ocv'   <-- switch this

[tests, ~, resultsDir] = testConfig();
modelName = 'step4_BatteryModel';
load(fullfile(resultsDir, 'Simscape_BatteryParams.mat'));   % SOC_bp,T_bp,OCV_tab,R0/1/2_*,tau*,MaxHyst_*
Lc = load(fullfile(resultsDir, 'LUT2D_Charge_2RC.mat'));    % Lc.Q_Ah (per-T cell capacity)

% ---- choose which characterization CSV to use ----
switch lower(testType)
    case 'capacity', fld = 'capacity';
    case 'ocv',      fld = 'ocv_dis';   % GITT discharge; change to 'ocv_chg' for charge side
    otherwise, error('testType must be ''capacity'' or ''ocv''.');
end

% locate the .slx
scriptFolder = fileparts(mfilename('fullpath'));
rootDirs = unique({ resultsDir, pwd, scriptFolder, fileparts(scriptFolder) });
mdlPath = '';
for r = 1:numel(rootDirs)
    if ~isfolder(rootDirs{r}), continue; end
    d = dir(fullfile(rootDirs{r}, '**', [modelName '.slx']));
    if ~isempty(d), mdlPath = fullfile(d(1).folder, d(1).name); break; end
end
if isempty(mdlPath), error('Model %s.slx not found.', modelName); end
fprintf('Using model: %s  (testType = %s)\n', mdlPath, testType);
if ~bdIsLoaded(modelName), load_system(mdlPath); end

battPath = '';
try
    bl = find_system(modelName, 'Type', 'Block');
    bb = bl(contains(lower(bl), 'battery'));
    if ~isempty(bb), battPath = bb{1}; end
catch
end

% variable-step daessc, tuned for the long profile (no fixed-step local solver)
set_param(modelName, 'SolverType', 'Variable-step', 'Solver', 'daessc', ...
          'RelTol', '1e-3', 'AbsTol', '1e-4', 'MinStep', '1e-2', 'MaxStep', '10');
try, set_param(modelName, 'ZeroCrossControl', 'Adaptive');  end
try, set_param(modelName, 'ConsecutiveZCs', 'none');      end
try, set_param(modelName, 'MaxConsecutiveZCs', '10000'); end

% initial SOC = 1 (full charge)
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

%% ==================== validate chosen test ============================
hasTest = false(1, numel(tests));
for k = 1:numel(tests)
    if isfield(tests(k), fld) && ~isempty(tests(k).(fld)) && isfile(tests(k).(fld))
        hasTest(k) = true;
    end
end
idx = find(hasTest);
if isempty(idx), error('No %s test files found (check testConfig / data folder).', testType); end

rmse_mv   = nan(1, numel(idx));
rmse_pct  = nan(1, numel(idx));

fig = figure('Color','w','Position',[30 30 1180 270*max(numel(idx),1)]);
nTiles = 0;

fprintf('\n=== Step 8 %s REPLAY — CELL model (step4_BatteryModel) ===\n', upper(testType));
for ii = 1:numel(idx)
    k  = idx(ii);
    Tk = tests(k).T;
    T  = readtable(tests(k).(fld), 'VariableNamingRule','preserve');
    t  = T{:,7}/1000;  V = T{:,8};  I = T{:,9};     % cell-level: t[s], V[V], I[A] (I>0 = charge)

    fprintf('  T = %d C : %d points, I range [%.2f, %.2f] A\n', Tk, numel(t), min(I), max(I));

    % ---- capacity / OCV tests: isolate the relevant sweep ----
    % The raw CSV usually includes a charge-to-full precondition + rests, so the
    % whole-profile replay from init-SOC=1 desyncs the sim. Trim to the part that
    % actually exercises the cell from a known full-charge state.
    if strcmpi(testType, 'capacity') || strcmpi(testType, 'ocv')
        try
            vnames = T.Properties.VariableNames;
            im = find(contains(lower(string(vnames)), 'stepmode'), 1);
            if isempty(im), im = 5; end                  % fallback: col 5
            StepMode = T{:, im};
            sm = lower(cellstr(StepMode));               % 'REST','CC Charge','CV Charge','CC Discharge'
            isDis = contains(sm, 'discharge') | contains(sm, 'dchg');
            if strcmpi(testType, 'capacity')
                seg = longestRun(isDis);                 % single dominant discharge step
            else
                f = find(isDis, 1, 'first');             % GITT: first pulse to end of sweep
                if isempty(f), seg = []; else seg = [f, numel(isDis)]; end
            end
            if ~isempty(seg)
                t = t(seg(1):seg(2)); V = V(seg(1):seg(2)); I = I(seg(1):seg(2));
                t = t - t(1);
                fprintf('     trimmed to discharge step: %d points (%.1f .. %.1f s)\n', numel(t), t(1), t(end));
            else
                fprintf('     WARNING: no discharge step found; using full profile.\n');
            end
        catch
            fprintf('     WARNING: step-column parse failed; using full profile.\n');
        end
    end

    T_sim = Tk + 273.15;
    jj = find(Lc.Tbp == Tk, 1);
    if ~isempty(jj), Cap_Ah = Lc.Q_Ah(jj); end          % CELL capacity at this T
    if ~isempty(battPath)
        try
            d = get_param(battPath,'DialogParameters'); pnm = fieldnames(d); plc = lower(pnm);
            ic = find(contains(plc,'capac'),1);
            if ~isempty(ic), set_param(battPath, pnm{ic}, num2str(Cap_Ah)); end
        catch
        end
    end

    I_in = timeseries(I, t, 'Name', 'I_in');             % cell current, no pack scaling

    set_param(modelName, 'StopTime', num2str(max(t)));
    assignin('base', 'V_out', []);
    warning('off','all');
    out = sim(modelName, 'ReturnWorkspaceOutputs', 'on');
    warning('on','all');
    Vsim_ts = [];
    try, Vsim_ts = out.V_out;        end
    if isempty(Vsim_ts), try, Vsim_ts = out.get('V_out'); end, end
    if isempty(Vsim_ts), try, Vsim_ts = evalin('base', 'V_out'); end, end
    if isempty(Vsim_ts), error('V_out not captured at T=%d C.', Tk); end
    V_sim = interp1(Vsim_ts.Time, Vsim_ts.Data, t, 'linear', 'extrap');

    % ---- DIAG: is the block's OCV pinned at SOC=1? (cell-level) ----
    V_ocv_data = NaN(size(t));
    try
        Ahc = cumtrapz(t(:), I(:)) / 3600;                % Ah, charge-positive (cell)
        SOCd = 1 - (max(Ahc) - Ahc)/(Cap_Ah);
        SOCd = max(0, min(1, SOCd));
        V_ocv_data = interp2(SOC_bp, T_bp, OCV_tab.', SOCd, T_sim);
        rmse_ocv = sqrt(mean((V_sim(:) - V_ocv_data(:)).^2))*1e3;
        fprintf('     DIAG  RMSE(sim vs OCV@dataSOC)=%.1f mV ; corr(sim,OCV)=%.3f ; corr(sim,meas)=%.3f\n', ...
                rmse_ocv, corr(V_sim(:), V_ocv_data(:)), corr(V_sim(:), V(:)));
    catch ME
        fprintf('     DIAG  skipped (%s)\n', ME.message);
    end

    % ---- compare cell sim to measured cell voltage (no pack scaling) ----
    res  = V_sim - V;
    rmse_mv(ii)  = sqrt(mean(res.^2))*1e3;
    rmse_pct(ii) = (sqrt(mean(res.^2))/Vmean_ref(V))*100;
    fprintf('     RMSE = %.3f mV (%.3f%%)\n', rmse_mv(ii), rmse_pct(ii));

    % ---- test-specific metric ----
    if strcmpi(testType, 'capacity')
        Cap_meas = -trapz(t, I)/3600;                    % discharged Ah (cell)
        fprintf('     measured cell capacity = %.3f Ah ; model cell cap = %.3f Ah\n', Cap_meas, Cap_Ah);
    elseif strcmpi(testType, 'ocv')
        iRest = abs(I) < 0.05*max(abs(I(:)));            % rest segments -> V ~= OCV
        if any(iRest)
            rmse_rest = sqrt(mean((V_sim(iRest) - V(iRest)).^2))*1e3;
            fprintf('     OCV (rest-segment) RMSE = %.2f mV\n', rmse_rest);
        end
    end

    nTiles = nTiles + 1;
    subplot(numel(idx), 2, nTiles); hold on; grid on
    plot(t/3600, V, 'k',  'LineWidth', 1.4);
    plot(t/3600, V_sim,'r--','LineWidth', 1.2);
    ylabel('V [V]');
    title(sprintf('%s @ %d^{\\circ}C — Cell (RMSE = %.2f mV, %.2f%%)', ...
        upper(testType), Tk, rmse_mv(ii), rmse_pct(ii)));
    legend('Measured','Simscape','Location','best');

    nTiles = nTiles + 1;
    subplot(numel(idx), 2, nTiles); hold on; grid on
    plot(t/3600, res*1e3, 'b');
    ylabel('Error [mV]'); xlabel('Time [h]');
    title(sprintf('Residuals @ %d^{\\circ}C', Tk));
end
saveas(fig, fullfile(resultsDir, sprintf('step8_%s_replay_cell.png', testType)), 'png');
fprintf('---------------------------------------------------------------\n');
fprintf('Mean %s RMSE : %.3f mV (%.3f%%)\n', upper(testType), ...
    mean(rmse_mv,'omitnan'), mean(rmse_pct,'omitnan'));

%% local helper: [start,end] indices of the longest contiguous true run in a logical mask
function seg = longestRun(mask)
    if ~any(mask), seg = []; return; end
    d = diff([0; mask(:); 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    [~, k] = max(ends - starts + 1);
    seg = [starts(k), ends(k)];
end
