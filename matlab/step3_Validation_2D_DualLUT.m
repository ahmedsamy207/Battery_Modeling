%% step3_Validation_2D_DualLUT.m
%  Stage 3 (2D, dual): validate the separate CHARGE & DISCHARGE 2-D LUTs
%  against every available HPPC temperature, with dynamic LUT switching
%  on the sign of the current (I > 0 -> charge LUT), like your 1-D script.
%
%  RMSE is reported in [mV] AND as a percentage of the mean measured cell
%  voltage, the normalization most commonly used in cell-level ECM
%  validation reports:
%        RMSE[%] = RMSE[V] / mean(V_measured) * 100
%  (If you prefer textbook NRMSE-by-range, replace Vmean below with
%   (max(V) - min(V)).)

clear; clc; close all

%% 1) LOAD LUTs & TEST MATRIX
[tests, ~, resultsDir] = testConfig();
Lc = load(fullfile(resultsDir, 'LUT2D_Charge_2RC.mat'));
Ld = load(fullfile(resultsDir, 'LUT2D_Discharge_2RC.mat'));

%% 2) VALIDATE EACH TEMPERATURE THAT HAS HPPC DATA
idx      = find([tests.has_hppc]);
rmse_all = nan(1, numel(idx));
rmse_chg = nan(1, numel(idx));
rmse_dis = nan(1, numel(idx));
rmse_pct_all = nan(1, numel(idx));
rmse_pct_chg = nan(1, numel(idx));
rmse_pct_dis = nan(1, numel(idx));

fig = figure('Color','w','Position',[30 30 1180 270*max(numel(idx),1)]);
nTiles = 0;

fprintf('\n--- 2-D DUAL-LUT VALIDATION (HPPC at each temperature) ---\n');
for ii = 1:numel(idx)
    k  = idx(ii);
    Tk = tests(k).T;
    j  = find(Ld.Tbp == Tk, 1);
    assert(~isempty(j), 'Temperature %d C missing in LUT breakpoints!', Tk);

    T = readtable(tests(k).hppc, 'VariableNamingRule', 'preserve');
    t = T{:,7}/1000;  V = T{:,8};  I = T{:,9};
    Vmean = mean(V);                       % reference voltage for % RMSE

    dt  = [0; diff(t)]; dt(dt < 0) = 0;
    Ah  = cumsum(I .* dt) / 3600;
    SOC = 1 - (max(Ah) - Ah) / Ld.Q_Ah(j);
    SOC = max(0, min(1, SOC));

    V_sim = zeros(size(t));
    V_c1  = zeros(size(t));
    V_c2  = zeros(size(t));
    V_sim(1) = lut2(Ld, Ld.OCV_V, SOC(1), Tk) + I(1) * lut2(Ld, Ld.R0_ohm, SOC(1), Tk);

    for n = 2:numel(t)
        % --- dynamic LUT switching ------------------------------------
        if I(n) > 0, L = Lc; else, L = Ld; end

        OCV_k  = lut2(L, L.OCV_V,  SOC(n), Tk);
        R0_k   = lut2(L, L.R0_ohm, SOC(n), Tk);
        R1_k   = lut2(L, L.R1_ohm, SOC(n), Tk);
        tau1_k = max(lut2(L, L.tau1_s, SOC(n), Tk), 0.001);
        R2_k   = lut2(L, L.R2_ohm, SOC(n), Tk);
        tau2_k = max(lut2(L, L.tau2_s, SOC(n), Tk), 0.001);

        decay1 = exp(-dt(n)/tau1_k);
        V_c1(n) = V_c1(n-1)*decay1 + R1_k*I(n)*(1-decay1);
        decay2 = exp(-dt(n)/tau2_k);
        V_c2(n) = V_c2(n-1)*decay2 + R2_k*I(n)*(1-decay2);

        V_sim(n) = OCV_k + I(n)*R0_k + V_c1(n) + V_c2(n);
    end

    res = V_sim - V;
    iCat = I > 0.5;   iDat = I < -0.5;
    rmse_all(ii) = sqrt(mean(res.^2))*1e3;
    rmse_chg(ii) = sqrt(mean(res(iCat).^2))*1e3;
    rmse_dis(ii) = sqrt(mean(res(iDat).^2))*1e3;

    % --- RMSE as a percentage of the mean cell voltage ---------------
    rmse_pct_all(ii) = (rmse_all(ii)/1e3) / Vmean * 100;
    rmse_pct_chg(ii) = (rmse_chg(ii)/1e3) / Vmean * 100;
    rmse_pct_dis(ii) = (rmse_dis(ii)/1e3) / Vmean * 100;

    fprintf('  T = %2d C : overall RMSE = %.3f mV (%.3f%%) | chg = %.3f mV (%.3f%%) | dis = %.3f mV (%.3f%%)\n', ...
            Tk, rmse_all(ii), rmse_pct_all(ii), ...
                rmse_chg(ii), rmse_pct_chg(ii), ...
                rmse_dis(ii), rmse_pct_dis(ii));

    nTiles = nTiles + 1;
    subplot(numel(idx), 2, nTiles); hold on; grid on
    plot(t/3600, V,     'k',   'LineWidth', 1.2);
    plot(t/3600, V_sim, 'r--', 'LineWidth', 1.0);
    ylabel('V [V]');
    title(sprintf('Dual-LUT HPPC @ %d^{\\circ}C (RMSE = %.2f mV, %.2f%%)', ...
            Tk, rmse_all(ii), rmse_pct_all(ii)));
    legend('Measured','Simulated','Location','best');

    nTiles = nTiles + 1;
    subplot(numel(idx), 2, nTiles); hold on; grid on
    plot(t/3600, res*1e3, 'b');
    ylabel('Error [mV]'); xlabel('Time [h]');
    title(sprintf('Residuals @ %d^{\\circ}C', Tk));
end
savefig(fig, fullfile(resultsDir, 'step3_Validation_2D_DualLUT.fig'));
print(fig, fullfile(resultsDir, 'step3_Validation_2D_DualLUT.png'), '-dpng', '-r200');

fprintf('---------------------------------------------------------------\n');
fprintf('Mean overall RMSE : %.3f mV (%.3f%%) | mean chg : %.3f mV (%.3f%%) | mean dis : %.3f mV (%.3f%%)\n', ...
        mean(rmse_all,'omitnan'),     mean(rmse_pct_all,'omitnan'), ...
        mean(rmse_chg,'omitnan'),     mean(rmse_pct_chg,'omitnan'), ...
        mean(rmse_dis,'omitnan'),     mean(rmse_pct_dis,'omitnan'));

%% LOCAL FUNCTION: 2-D lookup with clamping --------------------------------
function v = lut2(L, M, soc, Tq)
    soc = min(max(soc, L.soc_grid(1)),   L.soc_grid(end));
    Tq  = min(max(Tq,  L.Tbp(1)),        L.Tbp(end));
    if numel(L.Tbp) == 1
        v = interp1(L.soc_grid, M(:,1), soc, 'linear');
    else
        F = griddedInterpolant({L.soc_grid, L.Tbp}, M, 'linear', 'nearest');
        v = F(soc, Tq);
    end
end
