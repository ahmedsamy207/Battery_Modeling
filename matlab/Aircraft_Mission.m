%% step9_Mission.m
%  Mission-profile simulation of the SCRIPTED 192S4P pack (batteryPack.slx +
%  Simscape_BatteryPackParams.mat). missionType switches between:
%    'x57'   - NASA X-57 Maxwell profile (Table 1 of Battery_Evaluation_EATS):
%              Batt. Power [kW] + Duration [s] -> pack current = P / Vnom_pack.
%    'evtol' - representative eVTOL / UAM stress profile (synthetic): 2C takeoff
%              peak, steady cruise. Currents are pack currents [A].
%  Both report SOC, pack voltage, current and DC battery efficiency, with design
%  limits overlaid (2.75*Ns .. 4.2*Ns V, 2C = 2*Np*cellCap A, SOC floor 20%).
%
%  This targets the SCRIPTED pack. For the app pack (pack.slx) use pack_param and
%  the ModuleType1 struct instead (see the earlier step9 variant).

clear; clc; close all
missionType = 'x57';        % 'x57' (airplane baseline) | 'evtol' (stress)
modelName  = 'step4_BatteryModel';  % SCRIPTED pack
T_C        = 25;             % mission operating temperature [C]

% ---- mission profile selection ----
switch lower(missionType)
    case 'x57'
        Vnom_pack = 750;                          % V, converts X-57 powers -> pack current
        P_kW = [34.0, 187.8, 189.3, 156.8, 141.6, 31.4];
        Dur  = [300,    30,    15,    446,   417,   293];   % seconds
        names = {'Taxi','MotorCheck','TakeOff','Climb','Cruise','Descent'};
        tb = [0; cumsum(Dur(:))];
        t  = (0:tb(end))';
        I  = zeros(size(t));
        for p = 1:numel(P_kW)
            m = (t >= tb(p)) & (t < tb(p+1));
            I(m) = -P_kW(p)*1000 / Vnom_pack;      % discharge -> negative
        end
        % demo altitude [m] and ground-speed [m/s] segment end-values
        % (the browser mission dashboard uses the same demo assumptions)
        altSegEnd = [0, 0, 152, 2438, 2438, 15];   % 0 / 0 / 500 ft / 8000 ft / 8000 ft / ~50 ft
        vSegEnd   = [0, 0, 25,  60,   88,   55];   % m/s (cruise ~172 kt)
    case 'evtol'
        Vnom_pack = 750;
        I_A  = [-14,  -56,  -280,  -200,  -112,  -28,  -14];
        Dur  = [300,  180,   60,   180,   720,  180,  180];   % seconds
        names = {'Idle','Taxi','TakeOff','Climb','Cruise','Descent','Landing'};
        tb = [0; cumsum(Dur(:))];
        t  = (0:tb(end))';
        I  = zeros(size(t));
        for p = 1:numel(I_A)
            m = (t >= tb(p)) & (t < tb(p+1));
            I(m) = I_A(p);
        end
        % demo altitude [m] and ground-speed [m/s] segment end-values
        altSegEnd = [0, 0, 15, 300, 300, 30, 0];   % hover takeoff -> 300 m cruise -> hover landing
        vSegEnd   = [0, 2, 0.5, 50, 60, 15, 0];    % m/s
    otherwise
        error('missionType must be ''x57'' or ''evtol''.');
end
T_sim = T_C + 273.15;

% ---- results dir + SCRIPTED pack setup ----
try, [~,~,resultsDir] = testConfig(); catch, resultsDir = pwd; end
if ~isfolder(resultsDir), resultsDir = pwd; end
load(fullfile(resultsDir, 'Simscape_BatteryPackParams.mat'));   % Ns, Np, Q_Ah_p, T_bp, *_p tables

% design limits from pack scaling
Vmin_lim = 2.75 * Ns;       % discharge cutoff [V]
Vmax_lim = 4.20 * Ns;       % charge cutoff  [V]
jj = find(abs(T_bp - T_sim) < 0.5, 1);
if ~isempty(jj), Cap_Ah = Q_Ah_p(jj); end     % PACK capacity [Ah] at this T
Q_pack = Cap_Ah;                            % already pack-scaled
I_lim  = 2 * Cap_Ah;                        % 2C pack current [A] (= 2*Np*cellCap)

% ---- locate model + battery block ----
scriptFolder = fileparts(mfilename('fullpath'));
rootDirs = unique({ pwd, scriptFolder, resultsDir, fileparts(scriptFolder) });
mdlPath = '';
for r = 1:numel(rootDirs)
    if ~isfolder(rootDirs{r}), continue; end
    d = dir(fullfile(rootDirs{r}, '**', [modelName '.slx']));
    if ~isempty(d), mdlPath = fullfile(d(1).folder, d(1).name); break; end
end
if isempty(mdlPath), error('Model %s.slx not found.', modelName); end
fprintf('Using model: %s  (mission = %s)\n', mdlPath, missionType);
if ~bdIsLoaded(modelName), load_system(mdlPath); end
battPath = '';
try
    bl = find_system(modelName, 'Type', 'Block');
    bb = bl(contains(lower(bl), 'battery'));
    if ~isempty(bb), battPath = bb{1}; end
catch
end

% ---- solver (same tuned daessc as step5/step7) ----
set_param(modelName, 'SolverType','Variable-step','Solver','daessc', ...
          'RelTol','1e-3','AbsTol','1e-4','MinStep','1e-2','MaxStep','10');
try, set_param(modelName,'ZeroCrossControl','Adaptive'); end
try, set_param(modelName,'ConsecutiveZCs','none');        end
try, set_param(modelName,'MaxConsecutiveZCs','10000');    end

% ---- initial SOC = 1 + force per-T pack capacity ----
try
    if ~isempty(battPath)
        d = get_param(battPath,'DialogParameters'); nm = fieldnames(d); lc = lower(nm);
        is = find(contains(lc,'initial') & contains(lc,'soc'),1);
        if isempty(is), is = find(contains(lc,'soc0') | contains(lc,'soc1'),1); end
        if ~isempty(is), set_param(battPath, nm{is}, '1'); end
        ic = find(contains(lc,'capac'),1);
        if ~isempty(ic), set_param(battPath, nm{ic}, num2str(Cap_Ah)); end
    end
catch
end

% ---- simulate ----
I_in = timeseries(I, t, 'Name','I_in');
set_param(modelName, 'StopTime', num2str(t(end)));
assignin('base','V_out',[]);
warning('off','all');
out = sim(modelName, 'ReturnWorkspaceOutputs','on');
warning('on','all');
Vsim_ts = [];
try, Vsim_ts = out.V_out;        end
if isempty(Vsim_ts), try, Vsim_ts = out.get('V_out'); end, end
if isempty(Vsim_ts), try, Vsim_ts = evalin('base','V_out'); end, end
if isempty(Vsim_ts), error('V_out not captured.'); end
V = interp1(Vsim_ts.Time, Vsim_ts.Data, t, 'linear','extrap');

% ---- analysis ----
SOC   = 1 + cumtrapz(t,I)/3600 / Q_pack;     % I<0 (discharge) -> SOC decreases
P_del = -V .* I;                             % W, positive during discharge
E_del = trapz(t, P_del)/3.6e6;               % kWh delivered to load
dSOC  = SOC(1) - SOC(end);
E_from_cells = dSOC * Q_pack * Vnom_pack;    % Wh (stored-energy reference)
eta   = E_del*1000 / E_from_cells * 100;     % DC battery efficiency [%]
Imax  = max(abs(I));
Vmin  = min(V); Vmax = max(V);
SOCmin= min(SOC)*100;

% ---- demo altitude / range / speed + mission CSV export ----------------
%  The browser mission dashboard (dashboard/mission.html) replays
%  results/mission_<type>.csv when it exists; otherwise it builds the same
%  demo assumptions client-side. Altitude and ground speed are DEMO values.
alt_m  = zeros(size(t));   v_mps = zeros(size(t));
for p = 1:numel(names)
    if p == 1
        a0 = 0;  v0 = 0;
    else
        a0 = altSegEnd(p-1);  v0 = vSegEnd(p-1);
    end
    a1 = altSegEnd(p);  v1 = vSegEnd(p);
    m  = (t >= tb(p)) & (t <= tb(p+1));
    tt = (t(m) - tb(p)) / max(Dur(p), eps);
    alt_m(m) = a0 + (a1 - a0) .* tt;
    v_mps(m) = v0 + (v1 - v0) .* tt;
end
range_km = cumtrapz(t, v_mps) / 1000;               % integrate ground speed
P_dem    = -I .* Vnom_pack / 1000;                  % kW, positive = discharge demand
E_cum    = cumtrapz(t, P_del) / 3.6e6;              % kWh delivered to load
C_rate   = abs(I) / Q_pack;
E_rem    = SOC .* Q_pack .* Vnom_pack / 1000;       % kWh remaining (nominal-voltage reference)
phaseName = cell(size(t));
for p = 1:numel(names)
    m = (t >= tb(p)) & (t <= tb(p+1));
    [phaseName{m}] = deal(names{p});
end
MissionTbl = table(t, phaseName, P_dem, V, I, SOC*100, E_cum, alt_m, range_km, v_mps, C_rate, E_rem, ...
    'VariableNames', {'t_s','phase','P_demand_kW','V_pack','I_pack','SOC_pct', ...
                      'E_consumed_kWh','alt_m','range_km','v_mps','C_rate','E_remaining_kWh'});
writetable(MissionTbl, fullfile(resultsDir, sprintf('mission_%s.csv', lower(missionType))));
fprintf('  Exported mission CSV: results/mission_%s.csv (%d rows)\n', lower(missionType), height(MissionTbl));

% ---- plots ----
figure('Color','w','Position',[30 30 1000 760]);
subplot(4,1,1); plot(t/60, I, 'b','LineWidth',1.3); grid on; ylabel('I [A]'); hold on;
plot([0 t(end)/60],[-I_lim -I_lim],'r--'); text(0.2, -I_lim*0.98, '2C limit (-280 A)');
title(sprintf('%s mission - pack current (scripted)', upper(missionType)));
subplot(4,1,2); plot(t/60, V, 'b','LineWidth',1.3); grid on; ylabel('V [V]'); hold on;
plot([0 t(end)/60],[Vmin_lim Vmin_lim],'r--'); plot([0 t(end)/60],[Vmax_lim Vmax_lim],'r--');
title(sprintf('Pack voltage (window %.0f-%.0f V)', Vmin_lim, Vmax_lim));
subplot(4,1,3); plot(t/60, SOC*100, 'b','LineWidth',1.3); grid on; ylabel('SOC [%]'); hold on;
plot([0 t(end)/60],[20 20],'r--'); title('SOC (20 % floor)');
subplot(4,1,4); plot(t/60, P_del/1000, 'g','LineWidth',1.3); grid on;
ylabel('P [kW]'); xlabel('Time [min]'); title('Delivered power');

fprintf('\n=== %s mission on SCRIPTED 192S4P pack (T=%d C) ===\n', upper(missionType), T_C);
fprintf('  Duration: %.0f s (%.1f min)\n', t(end), t(end)/60);
fprintf('  Energy delivered: %.1f kWh ; SOC used: %.1f %% ; DC batt efficiency ~ %.1f %%\n', E_del, dSOC*100, eta);
fprintf('  Peak |I|: %.0f A (limit %.0f = 2C) ; V: %.0f..%.0f V (window %.0f..%.0f) ; min SOC: %.1f %%\n', Imax, I_lim, Vmin, Vmax, Vmin_lim, Vmax_lim, SOCmin);
saveas(gcf, fullfile(resultsDir, sprintf('step9_%s_mission_scriptedPack.png', missionType)), 'png');
