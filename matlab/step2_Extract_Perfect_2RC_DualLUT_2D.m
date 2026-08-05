%% step2_Extract_Perfect_2RC_DualLUT_2D.m
%  Stage 2 (2D, dual): separate Charge and Discharge 2-D LUTs.
%  Each parameter is a matrix param(SOC, T): rows = 201-point SOC grid,
%  columns = temperature breakpoints from step1.
%
%  Charge LUT pairs charge-pulse parameters with the CHARGE OCV branch;
%  Discharge LUT pairs discharge-pulse parameters with the DISCHARGE branch
%  (exactly like your 1-D dual script, now at every temperature).

clear; clc; close all

%% 1) LOAD STAGE-1 DNA & TEST MATRIX
[tests, ~, resultsDir] = testConfig();
S = load(fullfile(resultsDir, 'True_OCV_Capacity_2D.mat'));
soc_grid = S.soc_grid;
Tbp      = S.Tbp(:).';
Q_Ah     = S.Q_Ah(:).';
nT       = numel(Tbp);

%% 2) EXTRACTION PARAMETERS
pulse.maxDur_s  = 20;
pulse.minAbsI_A = 40;
pulse.minN      = 5;
pulse.fitWin_s  = 10;

%% 3) LOOP OVER TEMPERATURES: EXTRACT PER-PULSE PARAMETERS
perPulse = cell(1, nT);
fprintf('--- Step 2 (2D, dual LUT): HPPC pulse extraction ---\n');
for j = 1:nT
    Tk = Tbp(j);
    k  = find([tests.T] == Tk, 1);
    if isempty(k) || ~tests(k).has_hppc
        warning('No HPPC file for T = %d C -- LUT column stays NaN (clamped at lookup).', Tk);
        continue;
    end

    T    = readtable(tests(k).hppc, 'VariableNamingRule', 'preserve');
    Step = T{:,1}; t = T{:,7}/1000; V = T{:,8}; I = T{:,9};
    dt = [0; diff(t)]; dt(dt < 0) = 0;
    Ah = cumsum(I .* dt) / 3600;
    SOC_continuous = 1 - (max(Ah) - Ah) / Q_Ah(j);

    dS      = find(diff(Step)~=0);
    sStart  = [1; dS+1];  sEnd = [dS; numel(Step)];
    segMode = T{sStart,5};
    segDur  = t(sEnd)-t(sStart);
    segN    = sEnd-sStart+1;
    segI    = arrayfun(@(a,b) mean(I(a:b)), sStart, sEnd);
    SOCstart= SOC_continuous(sStart);

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

    perPulse{j} = cell2table(res, 'VariableNames', ...
        {'Step','Dir','SOC','R0','R1','tau1','R2','tau2'});
    fprintf('  T = %2d C : %3d pulses (%d dis, %d chg)\n', Tk, ...
        height(perPulse{j}), sum(perPulse{j}.Dir=="dis"), sum(perPulse{j}.Dir=="chg"));
end

%% 4) BUILD THE TWO 2D LUTs  (SOC x T)
paramNames = {'R0','R1','tau1','R2','tau2'};
gDis = struct(); gChg = struct();
for p = 1:numel(paramNames)
    gDis.(paramNames{p}) = nan(numel(soc_grid), nT);
    gChg.(paramNames{p}) = nan(numel(soc_grid), nT);
end

for j = 1:nT
    Tall = perPulse{j};
    if isempty(Tall), continue; end
    dis = Tall(Tall.Dir=="dis",:);
    chg = Tall(Tall.Dir=="chg",:);
    for p = 1:numel(paramNames)
        pn = paramNames{p};
        % linear inside measured coverage, constant edge extension outside
        gDis.(pn)(:, j) = interpGrid(dis.SOC, dis.(pn), soc_grid);
        gChg.(pn)(:, j) = interpGrid(chg.SOC, chg.(pn), soc_grid);
    end
end
gDis.C1 = gDis.tau1 ./ gDis.R1;   gDis.C2 = gDis.tau2 ./ gDis.R2;
gChg.C1 = gChg.tau1 ./ gChg.R1;   gChg.C2 = gChg.tau2 ./ gChg.R2;

LUT2D_Dis = assembleLUT(soc_grid, Tbp, Q_Ah, S.V_dis_all, gDis);
LUT2D_Chg = assembleLUT(soc_grid, Tbp, Q_Ah, S.V_chg_all, gChg);

save(fullfile(resultsDir, 'LUT2D_Discharge_2RC.mat'), '-struct', 'LUT2D_Dis');
save(fullfile(resultsDir, 'LUT2D_Charge_2RC.mat'),    '-struct', 'LUT2D_Chg');
fprintf('Saved: LUT2D_Discharge_2RC.mat & LUT2D_Charge_2RC.mat\n');

writeWideCSV(LUT2D_Dis, fullfile(resultsDir, 'LUT2D_Discharge_2RC.csv'));
writeWideCSV(LUT2D_Chg, fullfile(resultsDir, 'LUT2D_Charge_2RC.csv'));
fprintf('Saved: LUT2D_Discharge_2RC.csv & LUT2D_Charge_2RC.csv\n');

%% 5) QC PLOTS
fig = figure('Color','w','Position',[50 50 1250 420]);

subplot(1,3,1); hold on; grid on
for j = 1:numel(Tbp)
    if ~isempty(perPulse{j})
        d = perPulse{j}(perPulse{j}.Dir=="dis",:);
        plot(d.SOC*100, d.R0*1e3, 'o', 'DisplayName', sprintf('%d^{\\circ}C dis', Tbp(j)));
        plot(soc_grid*100, LUT2D_Dis.R0_ohm(:,j)*1e3, '-', 'HandleVisibility','off');
    end
end
xlabel('SOC [%]'); ylabel('R_0 [m\Omega]');
title('Discharge R_0: raw pulses + grid'); legend('Location','best');

subplot(1,3,2); hold on; grid on
for j = 1:numel(Tbp)
    if ~isempty(perPulse{j})
        d = perPulse{j}(perPulse{j}.Dir=="chg",:);
        plot(d.SOC*100, d.R0*1e3, 'o', 'DisplayName', sprintf('%d^{\\circ}C chg', Tbp(j)));
        plot(soc_grid*100, LUT2D_Chg.R0_ohm(:,j)*1e3, '-', 'HandleVisibility','off');
    end
end
xlabel('SOC [%]'); ylabel('R_0 [m\Omega]');
title('Charge R_0: raw pulses + grid'); legend('Location','best');

subplot(1,3,3); hold on; grid on
for j = 1:numel(Tbp)
    plot(soc_grid*100, LUT2D_Dis.OCV_V(:,j), '-', 'LineWidth', 1.5, ...
        'DisplayName', sprintf('%d^{\\circ}C dis-branch', Tbp(j)));
    plot(soc_grid*100, LUT2D_Chg.OCV_V(:,j), '--', 'HandleVisibility','off');
end
xlabel('SOC [%]'); ylabel('OCV [V]');
title('OCV branches (solid = dis, dashed = chg)'); legend('Location','best');

savefig(fig, fullfile(resultsDir, 'step2_LUT2D_DualLUT_QC.fig'));
print(fig, fullfile(resultsDir, 'step2_LUT2D_DualLUT_QC.png'), '-dpng', '-r200');

%% LOCAL FUNCTIONS --------------------------------------------------------
function g = interpGrid(xp, fp, xq)
% Linear interpolation inside the measured SOC coverage, CONSTANT
% (nearest-edge) extension outside. Safe for building LUT grids.
    [xp, ii] = unique(xp(:)); fp = fp(ii);
    g = interp1(xp, fp, xq, 'linear', nan);
    lo = isnan(g) & (xq < xp(1));   g(lo) = fp(1);
    hi = isnan(g) & (xq > xp(end)); g(hi) = fp(end);
    g  = fillmissing(g, 'nearest');
end

function L = assembleLUT(soc_grid, Tbp, Q_Ah, OCVmat, g)
    L.soc_grid = soc_grid;
    L.Tbp      = Tbp;
    L.Q_Ah     = Q_Ah;
    L.OCV_V    = OCVmat;
    L.R0_ohm   = g.R0;
    L.R1_ohm   = g.R1;
    L.C1_F     = g.C1;
    L.tau1_s   = g.tau1;
    L.R2_ohm   = g.R2;
    L.C2_F     = g.C2;
    L.tau2_s   = g.tau2;
end

function writeWideCSV(L, fpath)
    LUT = table(L.soc_grid, 'VariableNames', {'SOCbp'});
    for j = 1:numel(L.Tbp)
        tag = sprintf('_%dC', L.Tbp(j));
        LUT.(['OCV_V'  tag]) = L.OCV_V(:,j);
        LUT.(['R0_ohm' tag]) = L.R0_ohm(:,j);
        LUT.(['R1_ohm' tag]) = L.R1_ohm(:,j);
        LUT.(['C1_F'   tag]) = L.C1_F(:,j);
        LUT.(['tau1_s' tag]) = L.tau1_s(:,j);
        LUT.(['R2_ohm' tag]) = L.R2_ohm(:,j);
        LUT.(['C2_F'   tag]) = L.C2_F(:,j);
        LUT.(['tau2_s' tag]) = L.tau2_s(:,j);
    end
    writetable(LUT, fpath);
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
