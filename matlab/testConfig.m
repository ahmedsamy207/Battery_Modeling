function [tests, dataDir, resultsDir] = testConfig()
% testConfig  Central definition of all characterization datasets.
%             To add a temperature, extend the struct below (or just drop
%             the four CSV files into the data folder using the naming
%             convention "... <T> Celsius.csv" and add the T value).
%
%   tests(k) fields:
%       T        - nominal test temperature [degC]
%       capacity - capacity test CSV
%       ocv_chg  - stepwise (GITT-style) OCV charge CSV
%       ocv_dis  - stepwise (GITT-style) OCV discharge CSV
%       hppc     - hybrid pulse power characterization CSV
%
%   Any dataset whose files are missing is skipped automatically, so the
%   pipeline works with 1, 2 or N temperatures (2D LUTs need >= 2 for real
%   interpolation in the temperature dimension).

    % ---- locate folders -------------------------------------------------
    matlabDir  = fileparts(mfilename('fullpath'));
    rootDir    = fileparts(matlabDir);
    dataDir    = fullfile(rootDir, 'data');
    resultsDir = fullfile(rootDir, 'results');
    if ~exist(resultsDir, 'dir'), mkdir(resultsDir); end

    cellName = 'Solid State - Cell 2 - ';
    mk = @(kind, T) fullfile(dataDir, sprintf('%s%s %d Celsius.csv', cellName, kind, T));

    Tlist = [15 25 45];                 % nominal temperatures present/planned

    tests = struct('T', {}, 'capacity', {}, 'ocv_chg', {}, 'ocv_dis', {}, 'hppc', {});
    for k = 1:numel(Tlist)
        tests(k).T        = Tlist(k);
        tests(k).capacity = mk('Capacity',     Tlist(k));
        tests(k).ocv_chg  = mk('OCV Charge',   Tlist(k));
        tests(k).ocv_dis  = mk('OCV Discharge',Tlist(k));
        tests(k).hppc     = mk('HPPC',         Tlist(k));
    end

    % ---- flag which datasets actually exist on disk ---------------------
    for k = 1:numel(tests)
        tests(k).has_step1 = isfile(tests(k).capacity) & ...
                             isfile(tests(k).ocv_chg)  & isfile(tests(k).ocv_dis);
        tests(k).has_hppc  = isfile(tests(k).hppc);
    end
end
