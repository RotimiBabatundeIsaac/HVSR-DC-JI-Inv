%% test_single_station.m
clear; clc; close all;

clear c3_compute_likelihood c2_initialize_hbi_model forward_model

cps_bin = '/Users/birotimi/Downloads/HVSR-Joint-Inversion/AGU_models/bin_v3.30';
hv_bin = '/Users/birotimi/Downloads/HVSR-Joint-Inversion/AGU_models/bin_v3.30_new';
setenv('PATH', [cps_bin, ':', hv_bin, ':/usr/bin:/bin:/usr/sbin:/sbin']);

base_path = '/Users/birotimi/Downloads/HVSR-Joint-Inversion/Codes';
addpath(genpath(base_path));

old_functions_path = '/Users/birotimi/Downloads/HVSR-Joint-Inversion/Codes/Mod5_HBI/functions';
if exist(old_functions_path, 'dir')
    rmpath(old_functions_path);
end

new_functions_path = '/Users/birotimi/Downloads/HVSR-Joint-Inversion/Codes/Mod5_HBI/functions_hbi_new';
if exist(new_functions_path, 'dir')
    addpath(new_functions_path);
else
    error('New functions path not found: %s', new_functions_path);
end

addpath('/Users/birotimi/Downloads/HVSR-Joint-Inversion/AGU_models/functions');
addpath('/Users/birotimi/Downloads/HVSR-Joint-Inversion/Codes/Mod1_forward_models');

station_name = 'D20';

fprintf('MCMC Joint Inversion Test - Station %s\n', station_name);
fprintf('==========================================\n\n');

[~, result] = system('which surf96');
fprintf('surf96 binary: %s', result);

prior = c1_define_hbi_prior();

if ~isfield(prior, 'paths') || ~isfield(prior.paths, 'cps_bin')
    error('prior.paths.cps_bin not defined');
end
if ~isfield(prior.paths, 'hv_dfa_bin')
    error('prior.paths.hv_dfa_bin not defined');
end
if ~exist(prior.paths.cps_bin, 'dir')
    error('cps_bin directory not found: %s', prior.paths.cps_bin);
end
if ~exist(prior.paths.hv_dfa_bin, 'dir')
    error('hv_dfa_bin directory not found: %s', prior.paths.hv_dfa_bin);
end

fprintf('Prior bounds:\n');
for i = 1:prior.n_params
    fprintf('  %s: [%.4f, %.4f] %s\n', prior.param_names{i}, ...
        prior.bounds.lower(i), prior.bounds.upper(i), prior.param_units{i});
end
fprintf('\n');

datapack_path = '/Users/birotimi/Downloads/HVSR-Joint-Inversion/DATA/sta_datapack.mat';
if ~exist(datapack_path, 'file')
    error('Datapack not found: %s', datapack_path);
end
load(datapack_path, 'sta');

idx_sta = find(contains({sta.name}, station_name, 'IgnoreCase', true));
if isempty(idx_sta)
    error('Station %s not found in datapack', station_name);
end
data = sta(idx_sta(1));
fprintf('Station: %s (lat=%.4f, lon=%.4f)\n\n', data.name, data.lat, data.lon);

fprintf('Data summary:\n');
if isfield(data, 'hvsr') && isfield(data.hvsr, 'f') && ~isempty(data.hvsr.f)
    fprintf('  HVSR: %d points, f=[%.2f, %.2f] Hz\n', ...
        numel(data.hvsr.f), min(data.hvsr.f), max(data.hvsr.f));
end
if isfield(data, 'ellip') && isfield(data.ellip, 'T') && ~isempty(data.ellip.T)
    fprintf('  Ellipticity: %d points, T=[%.2f, %.2f] s\n', ...
        numel(data.ellip.T), min(data.ellip.T), max(data.ellip.T));
end
if isfield(data, 'cph') && isfield(data.cph, 'T') && ~isempty(data.cph.T)
    fprintf('  Phase velocity: %d points, T=[%.2f, %.2f] s\n', ...
        numel(data.cph.T), min(data.cph.T), max(data.cph.T));
end
if isfield(data, 'ugr') && isfield(data.ugr, 'T') && ~isempty(data.ugr.T)
    fprintf('  Group velocity: %d points, T=[%.2f, %.2f] s\n', ...
        numel(data.ugr.T), min(data.ugr.T), max(data.ugr.T));
end
fprintf('\n');

fprintf('Quick dispersion verification...\n');
T_test = [10, 15, 20, 25, 30];
h_test = [0.03, 0.17, 0.53, 10.27, 24.0];
vs_test = [0.50, 0.80, 0.99, 3.53, 3.73, 4.4];
model_test = zeros(6,4);
for i = 1:6
    if i <= 5, model_test(i,1) = h_test(i); else, model_test(i,1) = 0; end
    model_test(i,3) = vs_test(i);
    vp = 0.9409 + 2.0947*vs_test(i) - 0.8206*vs_test(i)^2 + 0.2683*vs_test(i)^3 - 0.0251*vs_test(i)^4;
    rho = 1.6612*vp - 0.4721*vp^2 + 0.0671*vp^3 - 0.0043*vp^4 + 0.000106*vp^5;
    model_test(i,2) = max(vp, vs_test(i)*1.5);
    model_test(i,4) = max(rho, 1.5);
end
model_ext = model_test;
model_ext(6,1) = 5;
for i = 1:10
    model_ext(end+1,:) = [10, 7.832, 4.4, 3.0];
end
model_ext(end,1) = 0;
c_R = dispR_surf96(T_test, model_ext, 0);
fprintf('  Phase vel at T=[10-30s]: [%.2f, %.2f] km/s (expected ~3.3-3.8)\n', min(c_R), max(c_R));
if min(c_R) < 2.0
    warning('Dispersion test failed! Check binary paths.');
else
    fprintf('  Dispersion test PASSED\n');
end
fprintf('\n');

fprintf('MCMC Settings:\n');
fprintf('  Iterations: %d\n', prior.mcmc.n_iterations);
fprintf('  Burn-in: %d\n', prior.mcmc.burn_in);
fprintf('  Thinning: %d\n', prior.mcmc.thin);
fprintf('  Chains: %d\n', prior.mcmc.n_chains);
if isfield(prior.mcmc, 'base_seed')
    fprintf('  Base seed: %d\n', prior.mcmc.base_seed);
end

n_samples = floor((prior.mcmc.n_iterations - prior.mcmc.burn_in) / prior.mcmc.thin);
fprintf('  Expected samples per chain: %d\n', n_samples);
fprintf('  Expected total samples: %d\n\n', n_samples * prior.mcmc.n_chains);

fprintf('Starting MCMC...\n');
fprintf('==========================================\n');
tic;

results = c8_run_multi_chain(data, prior);

elapsed = toc;
fprintf('\nMCMC complete in %.1f seconds (%.2f s/chain)\n', elapsed, elapsed/prior.mcmc.n_chains);

fprintf('\n==========================================\n');
fprintf('POSTERIOR SUMMARY\n');
fprintf('==========================================\n\n');

fprintf('Parameter estimates (mean +/- std) [90%% CI]:\n');
for i = 1:prior.n_params
    fprintf('  %s: %.4f +/- %.4f [%.4f, %.4f] %s\n', ...
        prior.param_names{i}, results.theta_mean(i), results.theta_std(i), ...
        results.theta_q05(i), results.theta_q95(i), prior.param_units{i});
end

fprintf('\nNoise scale estimates:\n');
for d = 1:prior.noise.n_noise
    fprintf('  %s: %.4f +/- %.4f\n', prior.noise.names{d}, ...
        results.sigma_mean(d), results.sigma_std(d));
end

fprintf('\nMoho depth: %.2f +/- %.2f km\n', ...
    sum(results.theta_mean(1:5)), ...
    sqrt(sum(results.theta_std(1:5).^2)));

fprintf('\nMAP estimate (max log-likelihood = %.2f):\n', results.map_logL);
fprintf('  h: %.3f %.3f %.3f %.3f %.3f km\n', results.theta_map(1:5));
fprintf('  Vs: %.3f %.3f %.3f %.3f %.3f %.3f km/s\n', results.theta_map(6:11));

fprintf('\nConvergence (Gelman-Rubin Rhat):\n');
for i = 1:prior.n_params
    if isnan(results.Rhat(i))
        status = 'N/A';
    elseif results.Rhat(i) > 1.1
        status = 'NOT CONVERGED';
    else
        status = 'OK';
    end
    fprintf('  %s: %.3f %s\n', prior.param_names{i}, results.Rhat(i), status);
end

n_converged = sum(results.Rhat < 1.1 & ~isnan(results.Rhat));
fprintf('\nConvergence: %d/%d parameters have Rhat < 1.1\n', n_converged, prior.n_params);

fprintf('\nComputing predictions for MAP model (best fit)...\n');
[~, i_map] = max(results.all_logL);
theta_map = results.all_theta(i_map, :)';
sigma_map = results.all_sigma(i_map, :)';
[~, ~, pred_map, ~, ~] = c3_compute_likelihood(theta_map, sigma_map, data, prior);

fprintf('\nVerifying MAP dispersion predictions:\n');
if ~isempty(pred_map.cph) && any(isfinite(pred_map.cph))
    fprintf('  Cph pred: [%.3f, %.3f] km/s, obs: [%.3f, %.3f] km/s\n', ...
        min(pred_map.cph(isfinite(pred_map.cph))), max(pred_map.cph(isfinite(pred_map.cph))), ...
        min(data.cph.obs), max(data.cph.obs));
end
if ~isempty(pred_map.ugr) && any(isfinite(pred_map.ugr))
    fprintf('  Ugr pred: [%.3f, %.3f] km/s, obs: [%.3f, %.3f] km/s\n', ...
        min(pred_map.ugr(isfinite(pred_map.ugr))), max(pred_map.ugr(isfinite(pred_map.ugr))), ...
        min(data.ugr.obs), max(data.ugr.obs));
end

fprintf('\nComputing spaghetti predictions...\n');
n_spag = min(100, results.n_total);
if n_spag < 1
    n_spag = 1;
end
spag_idx = randperm(results.n_total, n_spag);

pred_spag = cell(n_spag, 1);
for ii = 1:n_spag
    th_i = results.all_theta(spag_idx(ii), :)';
    sig_i = results.all_sigma(spag_idx(ii), :)';
    [~, ~, pred_spag{ii}, ~, ~] = c3_compute_likelihood(th_i, sig_i, data, prior);
    if mod(ii, 20) == 0
        fprintf('  Spaghetti prediction %d/%d\n', ii, n_spag);
    end
end
fprintf('Spaghetti predictions complete.\n');

f1 = figure('Position', [100 100 1200 400], 'Name', 'Chain_Diagnostics');

subplot(1,3,1);
hold on;
colors = lines(results.n_chains);
for c = 1:results.n_chains
    if isfield(results.chain_results{c}, 'samples_logL') && ~isempty(results.chain_results{c}.samples_logL)
        plot(results.chain_results{c}.samples_logL, 'Color', colors(c,:), 'LineWidth', 0.5);
    end
end
xlabel('Sample'); ylabel('Log-Likelihood');
title('Log-Likelihood Traces');
legend(arrayfun(@(x) sprintf('Chain %d', x), 1:results.n_chains, 'UniformOutput', false), 'Location', 'best');
grid on;

subplot(1,3,2);
hold on;
for c = 1:results.n_chains
    if isfield(results.chain_results{c}, 'samples_theta') && size(results.chain_results{c}.samples_theta, 2) >= 6
        plot(results.chain_results{c}.samples_theta(:,6), 'Color', colors(c,:), 'LineWidth', 0.5);
    end
end
xlabel('Sample'); ylabel('Vs1 (km/s)');
title('Vs1 Traces');
grid on;

subplot(1,3,3);
bar(results.Rhat);
hold on;
yline(1.1, 'r--', 'LineWidth', 1.5);
xlabel('Parameter Index'); ylabel('Rhat');
title('Gelman-Rubin Diagnostic');
xticks(1:prior.n_params);
xticklabels(prior.param_names);
xtickangle(45);
grid on;

sgtitle(sprintf('%s - Chain Diagnostics', data.name));

f2 = figure('Position', [100 100 1400 900], 'Name', 'Data_Fit_Spaghetti');

subplot(2,2,1);
if isfield(data, 'hvsr') && isfield(data.hvsr, 'f') && ~isempty(data.hvsr.f)
    f_obs = data.hvsr.f(:);
    hold on;
    for ii = 1:n_spag
        if ~isempty(pred_spag{ii}) && isfield(pred_spag{ii}, 'hvsr') && ~isempty(pred_spag{ii}.hvsr)
            y = pred_spag{ii}.hvsr(:);
            n = min(numel(y), numel(f_obs));
            ok = isfinite(y(1:n));
            if any(ok)
                plot(f_obs(ok), y(ok), '-', 'LineWidth', 0.5, 'Color', [0.75 0.75 0.75]);
            end
        end
    end
    if ~isempty(pred_map) && isfield(pred_map, 'hvsr') && ~isempty(pred_map.hvsr)
        y_map = pred_map.hvsr(:);
        n = min(numel(y_map), numel(f_obs));
        ok = isfinite(y_map(1:n));
        if any(ok)
            plot(f_obs(ok), y_map(ok), 'b-', 'LineWidth', 2);
        end
    end
    errorbar(f_obs, data.hvsr.obs(:), data.hvsr.sig(:), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', 1);
    set(gca, 'XScale', 'log');
    xlabel('Frequency (Hz)'); ylabel('HVSR Amplitude');
    title(sprintf('HVSR (%d points)', numel(f_obs)));
    legend('Posterior samples', 'MAP', 'Observed', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'No HVSR data', 'HorizontalAlignment', 'center');
    title('HVSR'); axis off;
end

subplot(2,2,2);
if isfield(data, 'ellip') && isfield(data.ellip, 'T') && ~isempty(data.ellip.T)
    T_obs = data.ellip.T(:);
    hold on;
    for ii = 1:n_spag
        if ~isempty(pred_spag{ii}) && isfield(pred_spag{ii}, 'ellip') && ~isempty(pred_spag{ii}.ellip)
            y = pred_spag{ii}.ellip(:);
            n = min(numel(y), numel(T_obs));
            ok = isfinite(y(1:n));
            if any(ok)
                plot(T_obs(ok), y(ok), '-', 'LineWidth', 0.5, 'Color', [0.75 0.75 0.75]);
            end
        end
    end
    if ~isempty(pred_map) && isfield(pred_map, 'ellip') && ~isempty(pred_map.ellip)
        y_map = pred_map.ellip(:);
        n = min(numel(y_map), numel(T_obs));
        ok = isfinite(y_map(1:n));
        if any(ok)
            plot(T_obs(ok), y_map(ok), 'r-', 'LineWidth', 2);
        end
    end
    errorbar(T_obs, data.ellip.obs(:), data.ellip.sig(:), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', 1);
    set(gca, 'XScale', 'log');
    xlabel('Period (s)'); ylabel('Ellipticity');
    title(sprintf('Ellipticity (%d points)', numel(T_obs)));
    legend('Posterior samples', 'MAP', 'Observed', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'No Ellipticity data', 'HorizontalAlignment', 'center');
    title('Ellipticity'); axis off;
end

subplot(2,2,3);
if isfield(data, 'cph') && isfield(data.cph, 'T') && ~isempty(data.cph.T)
    T_obs = data.cph.T(:);
    hold on;
    for ii = 1:n_spag
        if ~isempty(pred_spag{ii}) && isfield(pred_spag{ii}, 'cph') && ~isempty(pred_spag{ii}.cph)
            y = pred_spag{ii}.cph(:);
            n = min(numel(y), numel(T_obs));
            ok = isfinite(y(1:n));
            if any(ok)
                plot(T_obs(ok), y(ok), '-', 'LineWidth', 0.5, 'Color', [0.75 0.75 0.75]);
            end
        end
    end
    if ~isempty(pred_map) && isfield(pred_map, 'cph') && ~isempty(pred_map.cph)
        y_map = pred_map.cph(:);
        n = min(numel(y_map), numel(T_obs));
        ok = isfinite(y_map(1:n));
        if any(ok)
            plot(T_obs(ok), y_map(ok), 'g-', 'LineWidth', 2);
        end
    end
    errorbar(T_obs, data.cph.obs(:), data.cph.sig(:), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', 1);
    xlabel('Period (s)'); ylabel('Phase Velocity (km/s)');
    title(sprintf('Phase Velocity (%d points)', numel(T_obs)));
    legend('Posterior samples', 'MAP', 'Observed', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'No Phase Velocity data', 'HorizontalAlignment', 'center');
    title('Phase Velocity'); axis off;
end

subplot(2,2,4);
if isfield(data, 'ugr') && isfield(data.ugr, 'T') && ~isempty(data.ugr.T)
    T_obs = data.ugr.T(:);
    hold on;
    for ii = 1:n_spag
        if ~isempty(pred_spag{ii}) && isfield(pred_spag{ii}, 'ugr') && ~isempty(pred_spag{ii}.ugr)
            y = pred_spag{ii}.ugr(:);
            n = min(numel(y), numel(T_obs));
            ok = isfinite(y(1:n));
            if any(ok)
                plot(T_obs(ok), y(ok), '-', 'LineWidth', 0.5, 'Color', [0.75 0.75 0.75]);
            end
        end
    end
    if ~isempty(pred_map) && isfield(pred_map, 'ugr') && ~isempty(pred_map.ugr)
        y_map = pred_map.ugr(:);
        n = min(numel(y_map), numel(T_obs));
        ok = isfinite(y_map(1:n));
        if any(ok)
            plot(T_obs(ok), y_map(ok), 'm-', 'LineWidth', 2);
        end
    end
    errorbar(T_obs, data.ugr.obs(:), data.ugr.sig(:), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', 1);
    xlabel('Period (s)'); ylabel('Group Velocity (km/s)');
    title(sprintf('Group Velocity (%d points)', numel(T_obs)));
    legend('Posterior samples', 'MAP', 'Observed', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'No Group Velocity data', 'HorizontalAlignment', 'center');
    title('Group Velocity'); axis off;
end

sgtitle(sprintf('%s - Predicted vs Observed', data.name));

f3 = figure('Position', [100 100 600 700], 'Name', 'Velocity_Profile_Spaghetti');

h_map = theta_map(1:5);
vs_map = theta_map(6:11);

n_plot = min(200, results.n_total);
if n_plot < 1
    n_plot = 1;
end
plot_idx = randperm(results.n_total, n_plot);

hold on;
for i = 1:n_plot
    h_i = results.all_theta(plot_idx(i), 1:5);
    vs_i = results.all_theta(plot_idx(i), 6:11);
    plot_layer_profile(h_i(:), vs_i(:), 55, '-', 0.5, [0.75 0.75 0.75]);
end

plot_layer_profile(h_map(:), vs_map(:), 55, '-', 2.5, [1 0 0]);

set(gca, 'YDir', 'reverse');
xlabel('Vs (km/s)', 'FontSize', 12);
ylabel('Depth (km)', 'FontSize', 12);
title(sprintf('%s - Posterior Velocity Profiles', data.name), 'FontSize', 14);
xlim([0 5.5]);
ylim([0 55]);
grid on;
legend({'Posterior samples', 'MAP'}, 'Location', 'southwest');

fprintf('\nComputing posterior Vs(z) PDF...\n');

zmax = 50;
dz = 0.25;
z_grid = (0:dz:zmax)';
Nz = numel(z_grid);

vs_edges = 0:0.02:5.5;
vs_vec = 0.5 * (vs_edges(1:end-1) + vs_edges(2:end));
Nv = numel(vs_vec);

Ns = results.n_total;
h_samps = results.all_theta(:, 1:5);
vs_samps = results.all_theta(:, 6:11);

logL = results.all_logL;
logL_shift = logL - max(logL);
logL_shift(logL_shift < -700) = -700;
Lw = exp(logL_shift);
Lw = Lw / sum(Lw);

pdf_zvs = zeros(Nz, Nv);

for ii = 1:Ns
    h_i = h_samps(ii, :)';
    vs_i = vs_samps(ii, :)';
    vs_z = layer_vs_on_grid(z_grid, h_i, vs_i);
    ind_bin = discretize(vs_z, vs_edges);
    wi = Lw(ii);
    if ~isfinite(wi) || wi <= 0
        continue;
    end
    for iz = 1:Nz
        b = ind_bin(iz);
        if ~isnan(b) && b >= 1 && b <= Nv
            pdf_zvs(iz, b) = pdf_zvs(iz, b) + wi;
        end
    end
end

for iz = 1:Nz
    s = sum(pdf_zvs(iz, :));
    if s > 0
        pdf_zvs(iz, :) = pdf_zvs(iz, :) / s;
    end
end

vs_mean_at_z = nan(Nz, 1);
vs_std_at_z = nan(Nz, 1);
for iz = 1:Nz
    p = pdf_zvs(iz, :);
    if sum(p) <= 0
        continue;
    end
    p = p / sum(p);
    mu = sum(p .* vs_vec);
    vs_mean_at_z(iz) = mu;
    vs_std_at_z(iz) = sqrt(sum(p .* (vs_vec - mu).^2));
end

vs_map_profile = layer_vs_on_grid(z_grid, h_map(:), vs_map(:));

vs_2sig_lo = vs_mean_at_z - 2 * vs_std_at_z;
vs_2sig_hi = vs_mean_at_z + 2 * vs_std_at_z;
vs_2sig_lo(vs_2sig_lo < 0) = 0;
vs_2sig_hi(vs_2sig_hi > 5.5) = 5.5;

f4 = figure('Position', [150 80 700 800], 'Name', 'Posterior_Vs_PDF');
hold on; box on;

pdf_plot = pdf_zvs;
pdf_plot(pdf_plot == 0) = NaN;
pdf_log = log10(pdf_plot);

imagesc(vs_vec, z_grid, pdf_log);
set(gca, 'YDir', 'reverse');

colormap(parula(256));
cb = colorbar;
ylabel(cb, 'log_{10}(Probability)', 'FontSize', 11);

valid_pdf = pdf_log(isfinite(pdf_log));
if ~isempty(valid_pdf)
    cmax = max(valid_pdf);
    cmin = cmax - 5;
    caxis([cmin, cmax]);
end

plot(vs_2sig_lo, z_grid, 'w--', 'LineWidth', 1.5);
plot(vs_2sig_hi, z_grid, 'w--', 'LineWidth', 1.5);
plot(vs_map_profile, z_grid, 'r-', 'LineWidth', 2.5);

xlabel('V_S (km/s)', 'FontSize', 12);
ylabel('Depth (km)', 'FontSize', 12);
title(sprintf('%s - Posterior V_S(z) Distribution', data.name), 'FontSize', 14);
xlim([0, 5.5]);
ylim([0, zmax]);
grid on;
legend({'', '2-sigma bounds', '', 'MAP model'}, 'Location', 'southwest', 'TextColor', 'w');

fprintf('\n==========================================\n');
fprintf('DATA FIT STATISTICS (MAP MODEL)\n');
fprintf('==========================================\n');

fields = {'hvsr', 'ellip', 'cph', 'ugr'};
names = {'HVSR', 'Ellipticity', 'Phase Vel', 'Group Vel'};

for k = 1:4
    fld = fields{k};
    if isfield(pred_map, fld) && ~isempty(pred_map.(fld)) && isfield(data, fld)
        pred_k = pred_map.(fld)(:);
        obs_k = data.(fld).obs(:);
        sig_k = data.(fld).sig(:);
        n = min(numel(pred_k), numel(obs_k));
        valid = isfinite(pred_k(1:n)) & isfinite(obs_k(1:n));
        if any(valid)
            rms = sqrt(mean((pred_k(valid) - obs_k(valid)).^2));
            chi2 = sum(((pred_k(valid) - obs_k(valid))./sig_k(valid)).^2) / sum(valid);
            fprintf('%s: RMS=%.4f, Reduced Chi2=%.2f, N=%d\n', names{k}, rms, chi2, sum(valid));
        end
    end
end

fprintf('\n==========================================\n');
fprintf('SAVING RESULTS\n');
fprintf('==========================================\n');

fig_dir = '/Users/birotimi/Downloads/HVSR-Joint-Inversion/Codes/Mod5_HBI/figures_new_test';
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

station_tag = strrep(data.name, '.', '_');
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

results_file = fullfile(fig_dir, sprintf('%s_mcmc_results_%s.mat', station_tag, timestamp));
save(results_file, 'results', 'prior', 'data', 'pred_map', 'pdf_zvs', 'z_grid', 'vs_vec');
fprintf('Results saved: %s\n', results_file);

fig_handles = [f1, f2, f3, f4];
fig_names = {'Chain_Diagnostics', 'Data_Fit_Spaghetti', 'Velocity_Profile_Spaghetti', 'Posterior_Vs_PDF'};

for i = 1:length(fig_handles)
    pdf_path = fullfile(fig_dir, sprintf('%s_%s_%s.pdf', station_tag, fig_names{i}, timestamp));
    try
        exportgraphics(fig_handles(i), pdf_path, 'ContentType', 'vector', 'Resolution', 300);
        fprintf('Saved: %s\n', pdf_path);
    catch ME
        fprintf('Warning: Could not save %s: %s\n', fig_names{i}, ME.message);
    end
end

fprintf('\nAll outputs saved to: %s\n', fig_dir);
fprintf('\n==========================================\n');
fprintf('TEST COMPLETE\n');
fprintf('==========================================\n');

function vs_z = layer_vs_on_grid(z_grid, h, vs)
z_grid = z_grid(:);
h = h(:);
vs = vs(:);
nlay = numel(vs);
z_tops = [0; cumsum(h)];
z_bot = [cumsum(h); inf];
vs_z = nan(size(z_grid));
for k = 1:nlay
    in = (z_grid >= z_tops(k)) & (z_grid < z_bot(k));
    vs_z(in) = vs(k);
end
if any(isnan(vs_z))
    vs_z(isnan(vs_z)) = vs(end);
end
end

function plot_layer_profile(h, vs, zmax, ls, lw, col)
if nargin < 6, col = [0 0 1]; end
if nargin < 5, lw = 1.5; end
if nargin < 4, ls = '-'; end
h = h(:);
vs = vs(:);
nlay = numel(vs);
z = [0; cumsum(h)];
if numel(z) <= nlay
    z(end+1) = zmax;
end
zz = [];
vv = [];
for k = 1:nlay
    z0 = z(k);
    z1 = min(z(k+1), zmax);
    if z0 >= zmax
        break;
    end
    zz = [zz; z0; z1];
    vv = [vv; vs(k); vs(k)];
end
if ~isempty(zz)
    plot(vv, zz, ls, 'LineWidth', lw, 'Color', col);
end
end