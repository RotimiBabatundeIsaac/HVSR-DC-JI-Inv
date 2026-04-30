%% a1_0_collate_results.m
%  Check which stations have completed inversion results
%  Copy/verify results are in the common results directory

clc; clear;
fprintf('Running a1_0_collate_results from %s\n', pwd);

results_dir = '/Volumes/AS-Filer/EES/jbrussel/SharedData/HVSR-Joint-Inversion/Results_and_Figures/real_data_correlated_results';
datapack_path = '/Volumes/AS-Filer/EES/jbrussel/SharedData/HVSR-Joint-Inversion/DATA/sta_datapack.mat';

desired_chains = 10;
desired_iter = 24000;

%% Read all stations from datapack
load(datapack_path);
nsta_total = length(sta);
sta_list_all = cell(nsta_total, 1);
for is = 1:nsta_total
    sta_list_all{is} = strrep(sta(is).name, 'Z9.', '');
end

%% Check which stations have results
have_res = false(nsta_total, 1);
n_chains_found = zeros(nsta_total, 1);
best_logL = nan(nsta_total, 1);

for is = 1:nsta_total
    sta_name = sta_list_all{is};
    sta_dir = fullfile(results_dir, sta_name);
    if ~exist(sta_dir, 'dir'), continue; end

    files = dir(fullfile(sta_dir, 'hbi_corr_*.mat'));
    files = files(~contains({files.name}, 'checkpoint'));
    n_chains_found(is) = length(files);

    if n_chains_found(is) >= desired_chains
        have_res(is) = true;
        logLs = [];
        for k = 1:length(files)
            tmp = load(fullfile(files(k).folder, files(k).name));
            logLs(end+1) = tmp.results.best_logL;
        end
        best_logL(is) = max(logLs);
    end
end

%% Report results
fprintf('\n=== RESULTS SUMMARY ===\n');
fprintf('Total stations in datapack: %d\n', nsta_total);
fprintf('Stations with %d+ chains: %d\n', desired_chains, sum(have_res));
fprintf('Stations incomplete: %d\n', sum(n_chains_found > 0 & ~have_res));
fprintf('Stations with no results: %d\n', sum(n_chains_found == 0));

fprintf('\n--- Complete stations ---\n');
for is = 1:nsta_total
    if have_res(is)
        fprintf('  %s: %d chains, best logL=%.1f\n', sta_list_all{is}, n_chains_found(is), best_logL(is));
    end
end

fprintf('\n--- Incomplete stations ---\n');
for is = 1:nsta_total
    if n_chains_found(is) > 0 && ~have_res(is)
        fprintf('  %s: %d/%d chains\n', sta_list_all{is}, n_chains_found(is), desired_chains);
    end
end

fprintf('\n--- Missing stations ---\n');
missing = {};
for is = 1:nsta_total
    if n_chains_found(is) == 0
        fprintf('  %s\n', sta_list_all{is});
        missing{end+1} = sta_list_all{is};
    end
end

%% Save list of stations still needed
fid = fopen(fullfile(results_dir, 'station_list_still_need.txt'), 'w');
for is = 1:length(missing)
    fprintf(fid, '%s\n', missing{is});
end
fclose(fid);
fprintf('\nSaved %d missing stations to station_list_still_need.txt\n', length(missing));