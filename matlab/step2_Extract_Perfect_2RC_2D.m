%% step2_Extract_Perfect_2RC_2D.m
%  Stage 2 (2D version): HPPC parameter extraction for every temperature,
%  then assembly of 2-D LUTs: param(SOC, T)  for the averaged 2-RC ECM.
%
%  Output matrices are oriented (SOC x T): rows = soc_grid breakpoints,
%  columns = temperature breakpoints. This is exactly what Simscape /
%  Simulink 2-D Lookup Table blocks expect.

clear; clc; close all

%% 1) LOAD STAGE-1 DNA & TEST MATRIX
[tests, ~, resultsDir] = testConfig();
S = load(fullfile(resultsDir, 'True_OCV_Capacity_2D.mat'));   % soc_grid, Tbp, Q_Ah, True_all...
soc_grid = S.soc_grid;
Tbp      = S.Tbp(:).';                % row vector of temperatures
Q_Ah     = S.Q_Ah(:).';               % capacity at each temperature
nT       = numel(Tbp);

%% 2) EXTRACTION PARAMETERS (identical logic to your 25C script)
pulse.maxDur_s  = 20;     % a pulse is a short segment
pulse.minAbsI_A = 40;     % HPPC pulses are large (|70 A| here, threshold 40 A)
pulse.minN      = 5;      % must contain >= 5 samples (kills the aborted 10 ms pulse)
pulse.fitWin_s  = 10;     % fit the first 10 s of each pulse

%% 3) LOOP OVER TEMPERATURES: EXTRACT PER-PULSE PARAMETERS
perPulse = cell(1, nT);   % each cell: table of extracted pulse parameters

fprintf('--- Step 2 (2D): HPPC pulse extraction per temperature ---\n');
for j = 1:nT
    Tk = Tbp(j);
    k  = find([tests.T] == Tk, 1);
    if isempty(k) || ~tests(k).has_hppc
        warning('No HPPC file for T = %d C -- column %d of the LUTs will be NaN-safe extrapolated later.', Tk, j);
        continue;
    end

    T    = readtable(tests(k).hppc, 'VariableNamingRule', 'preserve');
    Step = T{:,1};
    t    = T{:,7}/1000;
    V    = T{:,8};
    I    = T{:,9};

    % SOC tracking with THIS temperature's true capacity
    dt = [0; diff(t)]; dt(dt < 0) = 0;
    Ah = cumsum(I .* dt) / 3600;
    SOC_continuous = 1 - (max(Ah) - Ah) / Q_Ah(j);

    % segment mapping
    dS      = find(diff(Step)~=0);
    sStart  = [1; dS+1];
    sEnd    = [dS; numel(Step)];
    segMode = T{sStart,5};
    segDur  = t(sEnd)-t(sStart);
    segN    = sEnd-sStart+1;
    segI    = arrayfun(@(a,b) mean(I(a:b)), sStart, sEnd);
    SOCstart= SOC_continuous(sStart);

    % pulse selection (same rule as the 25C script, now parameterized)
    isPulse = segDur < pulse.maxDur_s & abs(segI) > pulse.minAbsI_A & segN >= pulse.minN & ...
              (strcmp(segMode,'CC Discharge') | strcmp(segMode,'CC Charge'));
    plist = find(isPulse);

    res = [];
    for kk = plist.'
        a = sStart(kk); b = sEnd(kk);
        tp = t(a:b)-t(a);  Vp = V(a:b);  Ip = I(a:b);
        sel = tp <= pulse.fitWin_s;  tp = tp(sel);  Vp = Vp(sel);  Ip = Ip(sel);
        if numel(tp) < pulse.minN, continue; end

        if kk > 1 && strcmp(segMode{kk-1},'REST')
            OCV0 = V(sEnd(kk-1));
        else
            continue;
        end

        Ibar = mean(Ip);
        R0   = (Vp(1)-OCV0)/Ibar;
        fit  = fitRC(tp, Vp, OCV0, Ibar, R0);

        dir = "dis"; if Ibar > 0, dir = "chg"; end
        res = [res; {Step(a), dir, SOCstart(kk), R0, fit.R1, fit.tau1, fit.R2, fit.tau2}]; %#ok<AGROW>
    end

    Tall = cell2table(res, 'VariableNames', {'Step','Dir','SOC','R0','R1','tau1','R2','tau2'});
    perPulse{j} = Tall;
    fprintf('  T = %2d C : %3d pulses (%d dis, %d chg), SOC %.3f .. %.3f\n', ...
        Tk, height(Tall), sum(Tall.Dir=="dis"), sum(Tall.Dir=="chg"), ...
        min(Tall.SOC), max(Tall.SOC));
end

%% 4) MAP PARAMETERS ONTO THE 201-POINT SOC GRID, ONE COLUMN PER TEMPERATURE
paramNames = {'R0','R1','tau1','R2','tau2'};
grids = struct();
for p = 1:numel(paramNames)
    grids.(paramNames{p}) = nan(numel(soc_grid), nT);
end

for j = 1:nT
    Tall = perPulse{j};
    if isempty(Tall), continue; end
    dis = Tall(Tall.Dir=="dis",:);
    chg = Tall(Tall.Dir=="chg",:);
    for p = 1:numel(paramNames)
        pn = paramNames{p};
        % average of the two directions. interpGrid = linear inside the
        % measured SOC coverage, CONSTANT edge extension outside (never
        % slope-extrapolate noisy fits -- blows up the RC states, e.g. the
        % aborted first charge pulse at 15 C leaves a 15%-SOC hole).
        gi = (interpGrid(dis.SOC, dis.(pn), soc_grid) + ...
              interpGrid(chg.SOC, chg.(pn), soc_grid)) / 2;
        grids.(pn)(:, j) = gi;
    end
end

% capacitances from the averaged quantities
grids.C1 = grids.tau1 ./ grids.R1;
grids.C2 = grids.tau2 ./ grids.R2;

%% 5) ASSEMBLE THE 2D LUT (SOC x T) AND SAVE
LUT2D = struct();
LUT2D.soc_grid  = soc_grid;                 % 201 x 1
LUT2D.Tbp       = Tbp;                      % 1 x nT
LUT2D.Q_Ah      = Q_Ah;                     % capacity vs T
LUT2D.OCV_V     = S.True_all;               % averaged OCV branch  (SOC x T)
LUT2D.R0_ohm    = grids.R0;
LUT2D.R1_ohm    = grids.R1;
LUT2D.C1_F      = grids.C1;
LUT2D.tau1_s    = grids.tau1;
LUT2D.R2_ohm    = grids.R2;
LUT2D.C2_F      = grids.C2;
LUT2D.tau2_s    = grids.tau2;
LUT2D.perPulse  = perPulse;                 % raw extraction, for audit

fMAT = fullfile(resultsDir, 'LUT2D_Averaged_2RC.mat');
save(fMAT, '-struct', 'LUT2D');
fprintf('Saved: %s\n', fMAT);

% wide human-readable CSV (SOC rows, temperature blocks of columns)
LUT = table(soc_grid, 'VariableNames', {'SOCbp'});
for j = 1:nT
    tag = sprintf('_%dC', Tbp(j));
    LUT.(['OCV_V'  tag]) = LUT2D.OCV_V(:,j);
    LUT.(['R0_ohm' tag]) = LUT2D.R0_ohm(:,j);
    LUT.(['R1_ohm' tag]) = LUT2D.R1_ohm(:,j);
    LUT.(['C1_F'   tag]) = LUT2D.C1_F(:,j);
    LUT.(['tau1_s' tag]) = LUT2D.tau1_s(:,j);
    LUT.(['R2_ohm' tag]) = LUT2D.R2_ohm(:,j);
    LUT.(['C2_F'   tag]) = LUT2D.C2_F(:,j);
    LUT.(['tau2_s' tag]) = LUT2D.tau2_s(:,j);
end
writetable(LUT, fullfile(resultsDir, 'LUT2D_Averaged_batteryParams_2RC.csv'));
fprintf('Saved: %s\n', fullfile(resultsDir, 'LUT2D_Averaged_batteryParams_2RC.csv'));

% raw per-pulse extraction CSV (audit trail)
for j = 1:nT
    if ~isempty(perPulse{j})
        writetable(perPulse{j}, fullfile(resultsDir, ...
            sprintf('PerPulse_RAW_%dC.csv', Tbp(j))));
    end
end

%% 6) QC PLOTS: parameter surfaces vs SOC and T
plotLUT2D(LUT2D, 'Averaged', resultsDir);

%% LOCAL FUNCTIONS --------------------------------------------------------
function g = interpGrid(xp, fp, xq)
% Linear interpolation inside the measured SOC coverage, CONSTANT
% (nearest-edge) extension outside. Safe for building LUT grids.
    [xp, ii] = unique(xp(:)); fp = fp(ii);
    g = interp1(xp, fp, xq, 'linear', nan);
    lo = isnan(g) & (xq < xp(1));   g(lo) = fp(1);
    hi = isnan(g) & (xq > xp(end)); g(hi) = fp(end);
    g  = fillmissing(g, 'nearest');     % safety net (no NaN may remain)
end

function out = fitRC(tt, Vv, OCV0, Ibar, R0)
model = @(p) OCV0 + Ibar*( R0 + exp(p(1))*(1-exp(-tt/exp(p(2)))) ...
    + exp(p(3))*(1-exp(-tt/exp(p(4)))) );
p0 = log([2e-4 0.6 9e-4 11]);
p  = fminsearch(@(p) sum((model(p)-Vv).^2), p0, ...
    optimset('TolX',1e-12,'TolFun',1e-16,'MaxFunEvals',400000,'MaxIter',40000));
out.R1   = exp(p(1));  out.tau1 = exp(p(2));
out.R2   = exp(p(3));  out.tau2 = exp(p(4));
end

function plotLUT2D(L, styleTag, resultsDir)
fields2 = {'R0_ohm','R1_ohm','R2_ohm','tau1_s','tau2_s','C1_F','C2_F'};
ttl = {'R_0 [\Omega]','R_1 [\Omega]','R_2 [\Omega]', ...
       '\tau_1 [s]','\tau_2 [s]','C_1 [F]','C_2 [F]'};
nF = numel(fields2);

% --- surface (imagesc) overview: 3 x 4 grid, OCV spans the first row ---
fig = figure('Color','w','Position',[50 50 1250 780]);
subplot(3,4,1:4);
imagesc(L.Tbp, L.soc_grid*100, L.OCV_V); set(gca,'YDir','normal');
cb = colorbar; cb.Label.String = 'V';
xlabel('T [^{\circ}C]'); ylabel('SOC [%]');
title(sprintf('2-RC %s 2-D LUTs — OCV(SOC,T)', styleTag));
for f = 1:nF
    subplot(3,4,4+f);
    imagesc(L.Tbp, L.soc_grid*100, L.(fields2{f})); set(gca,'YDir','normal');
    colorbar; xlabel('T [^{\circ}C]'); ylabel('SOC [%]');
    title(ttl{f});
end
savePlots(fig, fullfile(resultsDir, sprintf('step2_LUT2D_%s', styleTag)));

% --- 1-D line overview per parameter (one curve per temperature) --------
fig = figure('Color','w','Position',[50 50 1250 780]);
subplot(3,4,1:4); hold on; grid on
for j = 1:numel(L.Tbp)
    plot(L.soc_grid*100, L.OCV_V(:,j), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%d^{\\circ}C', L.Tbp(j)));
end
ylabel('OCV [V]'); xlabel('SOC [%]');
title(sprintf('2-RC %s LUTs — curves per temperature', styleTag));
legend('Location','best');
for f = 1:nF
    subplot(3,4,4+f); hold on; grid on
    for j = 1:numel(L.Tbp)
        plot(L.soc_grid*100, L.(fields2{f})(:,j), 'LineWidth', 1.2, ...
            'DisplayName', sprintf('%d^{\\circ}C', L.Tbp(j)));
    end
    xlabel('SOC [%]'); title(ttl{f});
end
savePlots(fig, fullfile(resultsDir, sprintf('step2_LUT2D_lines_%s', styleTag)));
end

function savePlots(fig, fpathNoExt)
% works on old and new MATLAB releases
    savefig(fig, [fpathNoExt '.fig']);
    print(fig, [fpathNoExt '.png'], '-dpng', '-r200');
end
