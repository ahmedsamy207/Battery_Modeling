%% step1_Capacity_and_OCV_2D.m
%  Stage 1 (2D version): True Capacity Q(T) and True OCV(SOC,T)
%  Loops over every temperature defined in testConfig.m. Files that do
%  not exist yet (e.g. the 45C data) are skipped automatically, so you
%  can simply drop the CSVs in the data folder and re-run.

clear; clc; close all

%% 1) TEST MATRIX
[tests, ~, resultsDir] = testConfig();
soc_grid = (0:0.005:1)';          % 201-point common SOC grid (same as 25C pipeline)

%% 2) LOOP OVER TEMPERATURES
idx       = find([tests.has_step1]);
Tbp       = [];
Q_Ah      = [];
V_chg_all = [];                   % SOC x nT matrices
V_dis_all = [];
True_all  = [];

fprintf('--- Step 1 (2D): True Capacity & OCV per temperature ---\n');
for k = idx
    Tk = tests(k).T;

    % ---- 2a) capacity at this temperature ------------------------------
    T_cap  = readtable(tests(k).capacity, 'VariableNamingRule', 'preserve');
    t_cap  = T_cap{:,7}/1000;
    I_cap  = T_cap{:,9};
    dt_cap = [0; diff(t_cap)]; dt_cap(dt_cap < 0) = 0;
    Ah_cap = cumsum(I_cap .* dt_cap) / 3600;
    Qk     = max(Ah_cap) - min(Ah_cap);

    % ---- 2b) rested OCV points (charge & discharge hysteresis branches) -
    [SOC_chg, V_chg] = getRestedOCV(tests(k).ocv_chg, Qk, true);
    [SOC_dis, V_dis] = getRestedOCV(tests(k).ocv_dis, Qk, false);

    [SOC_chg, i1] = unique(SOC_chg); V_chg = V_chg(i1);
    [SOC_dis, i2] = unique(SOC_dis); V_dis = V_dis(i2);

    % ---- 2c) map both branches onto the common 201-point SOC grid ------
    % SAFE GRIDDING (never slope-extrapolate):
    %  - interpolate inside each branch's SOC coverage (NaN outside)
    %  - the thermodynamic average uses the single available branch where
    %    the other branch has no data (hysteresis collapses at the edges)
    %  - finally, flat-extend the nearest edge value to cover 0% and 100%
    V_chg_grid = interp1(SOC_chg, V_chg, soc_grid, 'linear', nan);
    V_dis_grid = interp1(SOC_dis, V_dis, soc_grid, 'linear', nan);
    True_OCV   = mean([V_chg_grid, V_dis_grid], 2, 'omitnan');
    True_OCV   = fillmissing(True_OCV, 'nearest');
    m = isnan(V_chg_grid); V_chg_grid(m) = True_OCV(m);
    V_chg_grid = fillmissing(V_chg_grid, 'nearest');
    m = isnan(V_dis_grid); V_dis_grid(m) = True_OCV(m);
    V_dis_grid = fillmissing(V_dis_grid, 'nearest');

    % ---- 2d) accumulate -------------------------------------------------
    Tbp       = [Tbp,  Tk];           %#ok<AGROW>
    Q_Ah      = [Q_Ah, Qk];           %#ok<AGROW>
    V_chg_all = [V_chg_all, V_chg_grid]; %#ok<AGROW>
    V_dis_all = [V_dis_all, V_dis_grid]; %#ok<AGROW>
    True_all  = [True_all,  True_OCV];   %#ok<AGROW>

    fprintf('  T = %2d C : Q_Ah = %.4f Ah,  %d chg + %d dis rest points\n', ...
            Tk, Qk, numel(SOC_chg), numel(SOC_dis));

    % legacy per-temperature file (same variable names as your 1-D scripts,
    % so the old validation code keeps working unchanged)
    Q_Ah_1T        = Qk;            %#ok<NASGU>
    soc_grid_1T    = soc_grid;      %#ok<NASGU>
    True_OCV_1T    = True_OCV;      %#ok<NASGU>
    V_chg_grid_1T  = V_chg_grid;    %#ok<NASGU>
    V_dis_grid_1T  = V_dis_grid;    %#ok<NASGU>
    fLeg = fullfile(resultsDir, sprintf('True_OCV_Capacity_%dC.mat', Tk));
    save(fLeg, 'Q_Ah_1T', 'soc_grid_1T', 'True_OCV_1T', 'V_chg_grid_1T', 'V_dis_grid_1T');
end

%% 3) SORT BY TEMPERATURE & SAVE THE 2D RESULT
[Tbp, ord] = sort(Tbp);
Q_Ah      = Q_Ah(ord);
V_chg_all = V_chg_all(:, ord);
V_dis_all = V_dis_all(:, ord);
True_all  = True_all(:, ord);

save(fullfile(resultsDir, 'True_OCV_Capacity_2D.mat'), ...
     'soc_grid', 'Tbp', 'Q_Ah', 'V_chg_all', 'V_dis_all', 'True_all');
fprintf('Saved: %s  (%d temperatures)\n', ...
        fullfile(resultsDir, 'True_OCV_Capacity_2D.mat'), numel(Tbp));

% human-readable 2D table (SOC x T)
Tbl = table(soc_grid, 'VariableNames', {'SOC'});
for j = 1:numel(Tbp)
    Tbl.(sprintf('OCVchg_V_%dC', Tbp(j)))  = V_chg_all(:, j);
    Tbl.(sprintf('OCVdis_V_%dC', Tbp(j)))  = V_dis_all(:, j);
    Tbl.(sprintf('OCVavg_V_%dC', Tbp(j)))  = True_all(:, j);
end
writetable(Tbl, fullfile(resultsDir, 'OCV_2D_SOCxT.csv'));

% capacity-vs-temperature table
writetable(table(Tbp(:), Q_Ah(:), 'VariableNames', {'T_degC','Q_Ah'}), ...
           fullfile(resultsDir, 'Capacity_vs_T.csv'));

%% 4) PLOT: OCV(SOC) FOR EVERY TEMPERATURE
c = lines(max(numel(Tbp), 1));
figure('Color','w','Position',[100 100 900 620]); hold on; grid on
for j = 1:numel(Tbp)
    plot(soc_grid*100, V_chg_all(:,j), '--', 'Color', c(j,:), 'LineWidth', 1, ...
         'HandleVisibility','off');
    plot(soc_grid*100, V_dis_all(:,j), ':',  'Color', c(j,:), 'LineWidth', 1, ...
         'HandleVisibility','off');
    plot(soc_grid*100, True_all(:,j),  '-',  'Color', c(j,:), 'LineWidth', 2, ...
         'DisplayName', sprintf('True OCV @ %d^{\\circ}C', Tbp(j)));
end
xlabel('State of Charge [%]'); ylabel('Voltage [V]');
title('True Thermodynamic OCV vs SOC (all temperatures)');
legend('Location','best');
savefig(gcf, fullfile(resultsDir, 'step1_OCV_all_temperatures.fig'));
print(gcf, fullfile(resultsDir, 'step1_OCV_all_temperatures.png'), '-dpng', '-r200');

%% LOCAL FUNCTION: extract only the rested points (unchanged from 25C code)
function [soc_pts, v_pts] = getRestedOCV(file, Q_Ah, is_charge)
    T = readtable(file, 'VariableNamingRule', 'preserve');
    t = T{:,7}/1000; I = T{:,9}; V = T{:,8};
    Step = T{:,1}; Mode = T{:,5};
    dt = [0; diff(t)]; dt(dt < 0) = 0;
    Ah = cumsum(I .* dt) / 3600;

    if is_charge
        soc = (Ah - min(Ah)) / Q_Ah;
    else
        soc = 1 + (Ah - max(Ah)) / Q_Ah;
    end

    dS = find(diff(Step)~=0);
    sEnd = [dS; numel(Step)];
    isRest = strcmp(Mode(sEnd), 'REST');
    restIdx = sEnd(isRest);

    soc_pts = soc(restIdx);
    v_pts   = V(restIdx);
end
