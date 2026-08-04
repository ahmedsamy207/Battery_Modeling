%% step3_Validation_2D.m
%  Stage 3 (2D): validate the AVERAGED 2-D LUT against the HPPC test of
%  EVERY available temperature. Uses griddedInterpolant over (SOC, T) with
%  clamped (nearest) extrapolation outside the measured temperature range.
%
%  RMSE is reported in [mV] AND as a percentage of the mean measured cell
%  voltage, the normalization most commonly used in cell-level ECM
%  validation reports:
%        RMSE[%] = RMSE[V] / mean(V_measured) * 100
%  (If you prefer textbook NRMSE-by-range, replace Vmean below with
%   (max(V) - min(V)).)

clear; clc; close all

%% 1) LOAD LUT & TEST MATRIX
[tests, ~, resultsDir] = testConfig();
L = load(fullfile(resultsDir, 'LUT2D_Averaged_2RC.mat'));   % soc_grid,Tbp,Q_Ah,OCV_V,R0_ohm,...

%% 2) VALIDATE EACH TEMPERATURE THAT HAS HPPC DATA
idx     = find([tests.has_hppc]);
rmse_mv = nan(1, numel(idx));
rmse_pct = nan(1, numel(idx));

fig = figure('Color','w','Position',[30 30 1180 270*max(numel(idx),1)]);
nTiles = 0;

fprintf('\n--- 2-D AVERAGED LUT VALIDATION (HPPC at each temperature) ---\n');
for ii = 1:numel(idx)
    k  = idx(ii);
    Tk = tests(k).T;
    j  = find(L.Tbp == Tk, 1);                       % LUT column for Q_Ah
    assert(~isempty(j), 'Temperature %d C missing in LUT breakpoints!', Tk);

    T = readtable(tests(k).hppc, 'VariableNamingRule', 'preserve');
    t = T{:,7}/1000;  V = T{:,8};  I = T{:,9};
    Vmean = mean(V);                      % reference voltage for % RMSE

    % SOC tracking with the capacity of THIS temperature
    dt  = [0; diff(t)]; dt(dt < 0) = 0;
    Ah  = cumsum(I .* dt) / 3600;
    SOC = 1 - (max(Ah) - Ah) / L.Q_Ah(j);
    SOC = max(0, min(1, SOC));

    % --- simulate the 2-RC model with 2D parameter lookup --------------
    V_sim = zeros(size(t));
    V_c1  = zeros(size(t));
    V_c2  = zeros(size(t));

    V_sim(1) = lut2(L, L.OCV_V,  SOC(1), Tk) + I(1) * lut2(L, L.R0_ohm, SOC(1), Tk);

    for n = 2:numel(t)
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

    res  = V_sim - V;
    rmse = sqrt(mean(res.^2));
    rmse_mv(ii) = rmse*1e3;
    rmse_pct(ii) = rmse / Vmean * 100;              % RMSE (V) as % of mean cell V
    fprintf('  T = %2d C : RMSE = %.3f mV (%.3f%%) , max |e| = %.3f mV\n', ...
            Tk, rmse*1e3, rmse_pct(ii), max(abs(res))*1e3);

    % --- plots per temperature -----------------------------------------
    nTiles = nTiles + 1;
    subplot(numel(idx), 2, nTiles); hold on; grid on
    plot(t/3600, V,     'k',   'LineWidth', 1.2);
    plot(t/3600, V_sim, 'r--', 'LineWidth', 1.0);
    ylabel('V [V]');
    title(sprintf('HPPC @ %d^{\\circ}C — measured vs model (RMSE = %.2f mV, %.2f%%)', ...
            Tk, rmse*1e3, rmse_pct(ii)));
    legend('Measured','Simulated','Location','best');

    nTiles = nTiles + 1;
    subplot(numel(idx), 2, nTiles); hold on; grid on
    plot(t/3600, res*1e3, 'b');
    ylabel('Error [mV]'); xlabel('Time [h]');
    title(sprintf('Residuals @ %d^{\\circ}C', Tk));
end
savefig(fig, fullfile(resultsDir, 'step3_Validation_2D.fig'));
print(fig, fullfile(resultsDir, 'step3_Validation_2D.png'), '-dpng', '-r200');

fprintf('---------------------------------------------------------------\n');
fprintf('Overall RMSE across temperatures : %.3f mV (%.3f%%) (mean of per-T RMSE)\n', ...
        mean(rmse_mv,'omitnan'), mean(rmse_pct,'omitnan'));

%% LOCAL FUNCTION: 2-D lookup with clamping --------------------------------
function v = lut2(L, M, soc, Tq)
%M is (SOC x T). Clamp queries into the grid, then interpolate.
    soc = min(max(soc, L.soc_grid(1)),   L.soc_grid(end));
    Tq  = min(max(Tq,  L.Tbp(1)),        L.Tbp(end));
    if numel(L.Tbp) == 1
        v = interp1(L.soc_grid, M(:,1), soc, 'linear');
    else
        F = griddedInterpolant({L.soc_grid, L.Tbp}, M, 'linear', 'nearest');
        v = F(soc, Tq);
    end
end
