% test_e1_run_hbi_mcmc.m
clear; clc; close all;

base_path = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/Codes';
addpath(genpath(base_path));

%addpath('/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/Codes/Mod5_HBI/functions_hbi');
%rmpath('/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/Codes/Mod5_HBI/functions');

% Ensure we are using the updated function set
rmpath('/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/Codes/Mod5_HBI/functions');
addpath('/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/Codes/Mod5_HBI/functions_hbi');
rehash;

fprintf('MCMC Test Run\n');

%% -------------------- PRIOR / MODEL (new c1-c6 compliant) --------------------
prior = c1_define_hbi_prior();

% Ensure fields expected by c2/c4/c5 are present (new c1 should already set these)
if ~isfield(prior, 'idx'); prior.idx = struct(); end
if ~isfield(prior.idx, 'h'); prior.idx.h = (1:5)'; end
if ~isfield(prior.idx, 'vs'); prior.idx.vs = (6:11)'; end

if ~isfield(prior, 'enforce_monotonicity'); prior.enforce_monotonicity = true; end
if ~isfield(prior, 'monotonic_start_idx'); prior.monotonic_start_idx = 6; end

model = c2_initialize_hbi_model(prior);

fprintf('========================================\n');
if prior.use_hierarchical
    fprintf('MODE: HIERARCHICAL BAYESIAN (Gaussian prior)\n');
    mode_str = 'HBI Gaussian';
    mode_tag = 'HBI';
else
    fprintf('MODE: SIMPLE MCMC (Uniform prior)\n');
    mode_str = 'Uniform Prior';
    mode_tag = 'Uniform';
end
fprintf('========================================\n\n');

%% -------------------- MCMC settings (use c1 unless overridden) --------------------
% Use c1 values by default; allow overrides here for test runs
n_iterations = prior.mcmc.n_iterations;
burn_in      = prior.mcmc.burn_in;
thin         = prior.mcmc.thin;

% Optional override (uncomment to force custom test settings)
n_iterations = 1500;
burn_in = 400;
thin = 25;

n_samples = floor((n_iterations - burn_in) / thin);

fprintf('Test MCMC Settings:\n');
fprintf('  Iterations: %d\n', n_iterations);
fprintf('  Burn-in: %d\n', burn_in);
fprintf('  Thinning: %d\n', thin);
fprintf('  Output samples: %d\n\n', n_samples);

%% -------------------- Indices (consistent with new c1/c2) --------------------
idx_h  = prior.idx.h(:)';
idx_vs = prior.idx.vs(:)';

n_h  = numel(idx_h);
n_vs = numel(idx_vs);

%% -------------------- Load data --------------------
datapack_path = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/DATA/sta_datapack.mat';
load(datapack_path, 'sta');

idx_sta = find(contains({sta.name}, 'D17', 'IgnoreCase', true));
if isempty(idx_sta)
    error('No station found matching D17 in sta_datapack.mat');
end
data = sta(idx_sta(1));
fprintf('Station: %s\n\n', data.name);

%% -------------------- Allocate outputs --------------------
samples.theta          = zeros(n_samples, prior.n_params);
samples.sigma_e        = zeros(n_samples, prior.noise.n_noise);
samples.log_likelihood = zeros(n_samples, 1);

trace.log_likelihood = zeros(n_iterations, 1);
trace.accept         = zeros(n_iterations, 1);

theta   = model.theta(:);
sigma_e = model.sigma_e(:);

if prior.use_hierarchical
    mu_theta    = model.mu_theta(:);
    Sigma_theta = model.Sigma_theta;
else
    mu_theta    = [];
    Sigma_theta = [];
end

%% -------------------- Initial likelihood / prior --------------------
fprintf('Computing initial likelihood...\n');
[log_L, ~, pred] = c3_compute_likelihood(theta, sigma_e, data, prior);

if prior.use_hierarchical
    [log_prior_val, ~] = c4_compute_log_prior_theta(theta, mu_theta, Sigma_theta, prior);
else
    [log_prior_val, ~] = c4_compute_log_prior_theta(theta, [], [], prior);
end
log_posterior = log_L + log_prior_val;

fprintf('Initial log-likelihood: %.2f\n', log_L);
fprintf('Initial log-prior: %.2f\n', log_prior_val);
fprintf('Initial log-posterior: %.2f\n\n', log_posterior);

%% -------------------- Run MCMC --------------------
fprintf('Running MCMC test (%d iterations)...\n', n_iterations);
accept_count = 0;
sample_idx = 0;
tic;

for iter = 1:n_iterations
    % Proposal step (must be consistent with bounds + monotonicity in prior)
    [theta_prop, is_valid, ~] = c5_perturb_theta(theta, prior);

    if is_valid
        [log_L_prop, ~, pred_prop] = c3_compute_likelihood(theta_prop, sigma_e, data, prior);

        if prior.use_hierarchical
            [log_prior_prop, ~] = c4_compute_log_prior_theta(theta_prop, mu_theta, Sigma_theta, prior);
        else
            [log_prior_prop, ~] = c4_compute_log_prior_theta(theta_prop, [], [], prior);
        end

        log_posterior_prop = log_L_prop + log_prior_prop;
        log_alpha = log_posterior_prop - log_posterior;

        if log(rand()) < log_alpha
            theta         = theta_prop;
            log_L         = log_L_prop;
            log_prior_val = log_prior_prop;
            log_posterior = log_posterior_prop;
            pred          = pred_prop;
            accept_count  = accept_count + 1;
            trace.accept(iter) = 1;
        end
    end

    % Hierarchical Gibbs updates
    if prior.use_hierarchical
        mu_theta    = c6_gibbs_update_mu_theta(theta, Sigma_theta, prior);
        Sigma_theta = c7_gibbs_update_Sigma_theta(theta, mu_theta, prior);
        sigma_e     = c8_gibbs_update_sigma_e(pred, data, prior);

        [log_prior_val, ~] = c4_compute_log_prior_theta(theta, mu_theta, Sigma_theta, prior);
        log_posterior = log_L + log_prior_val;
    end

    trace.log_likelihood(iter) = log_L;

    % Save thinned samples
    if iter > burn_in && mod(iter - burn_in, thin) == 0
        sample_idx = sample_idx + 1;
        if sample_idx <= n_samples
            samples.theta(sample_idx, :) = theta(:)';
            samples.sigma_e(sample_idx, :) = sigma_e(:)';
            samples.log_likelihood(sample_idx) = log_L;
        end
    end

    if mod(iter, 100) == 0
        fprintf('  Iter %d/%d, logL=%.1f, accept=%.1f%%\n', ...
            iter, n_iterations, log_L, 100*accept_count/iter);
    end
end

elapsed = toc;
fprintf('\nTest complete in %.1f seconds\n', elapsed);
fprintf('Acceptance rate: %.1f%%\n', 100*accept_count/n_iterations);
fprintf('Time per iteration: %.2f s\n\n', elapsed/n_iterations);

%% -------------------- Posterior summaries --------------------
if sample_idx < 5
    warning('Very few posterior samples saved (%d). Consider reducing thinning or increasing iterations.', sample_idx);
end

theta_mean   = mean(samples.theta(1:sample_idx, :), 1)';
theta_std    = std(samples.theta(1:sample_idx, :), 0, 1)';
sigma_e_mean = mean(samples.sigma_e(1:sample_idx, :), 1)';

h_mean  = theta_mean(idx_h);
vs_mean = theta_mean(idx_vs);

fprintf('\nPosterior mean model:\n');
fprintf('Layer thicknesses (km): %.3f %.3f %.3f %.3f %.3f\n', h_mean);
fprintf('Vs values (km/s): %.3f %.3f %.3f %.3f %.3f %.3f\n', vs_mean);
fprintf('Noise std (sigma_e): %.3f %.3f %.3f %.3f\n', sigma_e_mean);

[~, ~, pred_final] = c3_compute_likelihood(theta_mean, sigma_e_mean, data, prior);

%% -------------------- Trace / mean model plots --------------------
figure('Position', [100, 100, 1200, 400]);

subplot(1,3,1);
plot(trace.log_likelihood, 'b-');
xlabel('Iteration');
ylabel('Log-Likelihood');
title('Log-Likelihood Trace');
xline(burn_in, 'r--', 'LineWidth', 1.5);
grid on;

subplot(1,3,2);
plot(samples.theta(1:sample_idx, idx_vs(1)), 'b-');
hold on;
plot(samples.theta(1:sample_idx, idx_vs(min(4, n_vs))), 'r-');
xlabel('Sample');
ylabel('Velocity (km/s)');
title('Velocity Traces');
legend('Vs1', sprintf('Vs%d', min(4,n_vs)));
grid on;

subplot(1,3,3);
depths = [0; cumsum(h_mean)];
max_depth = sum(h_mean) + 10;
hold on;
for i = 1:n_vs
    if i < n_vs
        plot([vs_mean(i), vs_mean(i)], [depths(i), depths(i+1)], 'b-', 'LineWidth', 2);
        plot([vs_mean(i), vs_mean(i+1)], [depths(i+1), depths(i+1)], 'b-', 'LineWidth', 2);
    else
        plot([vs_mean(i), vs_mean(i)], [depths(end), max_depth], 'b-', 'LineWidth', 2);
    end
end
set(gca, 'YDir', 'reverse');
xlabel('Vs (km/s)');
ylabel('Depth (km)');
title('Mean Velocity Profile');
xlim([0 5]);
ylim([0 max_depth]);
grid on;

sgtitle(sprintf('%s - %s (%.1f%% accept, %.1f s)', data.name, mode_str, 100*accept_count/n_iterations, elapsed));

%% -------------------- Predicted vs Observed --------------------
figure('Position', [100 100 1400 900]);

subplot(2,2,1);
if isfield(data, 'hvsr') && isfield(data.hvsr, 'f') && ~isempty(data.hvsr.f)
    f_obs = data.hvsr.f(:);
    obs_hvsr = data.hvsr.obs(:);
    sig_hvsr = data.hvsr.sig(:);

    errorbar(f_obs, obs_hvsr, sig_hvsr, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', 1);
    hold on;

    if isfield(pred_final, 'hvsr') && ~isempty(pred_final.hvsr)
        pred_hvsr = pred_final.hvsr(:);
        valid_pred = isfinite(pred_hvsr);
        if any(valid_pred)
            plot(f_obs(valid_pred), pred_hvsr(valid_pred), 'b-', 'LineWidth', 2);
        end
    end

    set(gca, 'XScale', 'log');
    xlabel('Frequency (Hz)');
    ylabel('HVSR Amplitude');
    title('HVSR');
    legend('Observed', 'Predicted', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'No HVSR data', 'HorizontalAlignment', 'center');
    title('HVSR');
end

subplot(2,2,2);
if isfield(data, 'ellip') && isfield(data.ellip, 'T') && ~isempty(data.ellip.T)
    T_obs = data.ellip.T(:);
    obs_ellip = data.ellip.obs(:);
    sig_ellip = data.ellip.sig(:);

    errorbar(T_obs, obs_ellip, sig_ellip, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', 1);
    hold on;

    if isfield(pred_final, 'ellip') && ~isempty(pred_final.ellip)
        pred_ellip = pred_final.ellip(:);
        if length(pred_ellip) == length(T_obs)
            plot(T_obs, pred_ellip, 'r-', 'LineWidth', 2);
        else
            n_pred = length(pred_ellip);
            T_pred = logspace(log10(min(T_obs)), log10(max(T_obs)), n_pred)';
            plot(T_pred, pred_ellip, 'r-', 'LineWidth', 2);
        end
    end

    set(gca, 'XScale', 'log');
    xlabel('Period (s)');
    ylabel('Ellipticity (H/V)');
    title('Rayleigh Wave Ellipticity');
    legend('Observed', 'Predicted', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'No Ellipticity data', 'HorizontalAlignment', 'center');
    title('Ellipticity');
end

subplot(2,2,3);
if isfield(data, 'cph') && isfield(data.cph, 'T') && ~isempty(data.cph.T)
    T_obs = data.cph.T(:);
    obs_cph = data.cph.obs(:);
    sig_cph = data.cph.sig(:);

    errorbar(T_obs, obs_cph, sig_cph, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', 1);
    hold on;

    if isfield(pred_final, 'cph') && ~isempty(pred_final.cph)
        pred_cph = pred_final.cph(:);
        if length(pred_cph) == length(T_obs)
            plot(T_obs, pred_cph, 'g-', 'LineWidth', 2);
        else
            n_pred = length(pred_cph);
            T_pred = logspace(log10(min(T_obs)), log10(max(T_obs)), n_pred)';
            plot(T_pred, pred_cph, 'g-', 'LineWidth', 2);
        end
    end

    set(gca, 'XScale', 'log');
    xlabel('Period (s)');
    ylabel('Phase Velocity (km/s)');
    title('Phase Velocity Dispersion');
    legend('Observed', 'Predicted', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'No Phase Velocity data', 'HorizontalAlignment', 'center');
    title('Phase Velocity');
end

subplot(2,2,4);
if isfield(data, 'ugr') && isfield(data.ugr, 'T') && ~isempty(data.ugr.T)
    T_obs = data.ugr.T(:);
    obs_ugr = data.ugr.obs(:);
    sig_ugr = data.ugr.sig(:);

    errorbar(T_obs, obs_ugr, sig_ugr, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', 1);
    hold on;

    if isfield(pred_final, 'ugr') && ~isempty(pred_final.ugr)
        pred_ugr = pred_final.ugr(:);
        if length(pred_ugr) == length(T_obs)
            plot(T_obs, pred_ugr, 'm-', 'LineWidth', 2);
        else
            n_pred = length(pred_ugr);
            T_pred = logspace(log10(min(T_obs)), log10(max(T_obs)), n_pred)';
            plot(T_pred, pred_ugr, 'm-', 'LineWidth', 2);
        end
    end

    set(gca, 'XScale', 'log');
    xlabel('Period (s)');
    ylabel('Group Velocity (km/s)');
    title('Group Velocity Dispersion');
    legend('Observed', 'Predicted', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'No Group Velocity data', 'HorizontalAlignment', 'center');
    title('Group Velocity');
end

sgtitle(sprintf('Predicted vs Observed - %s [%s]', data.name, mode_str));

%% -------------------- Posterior velocity profile spaghetti --------------------
figure('Position', [100 100 600 700]);

depths = [0; cumsum(h_mean)];
max_depth = sum(h_mean) + 10;

n_plot = min(200, sample_idx);
if n_plot > 0
    plot_idx = randperm(sample_idx, n_plot);

    hold on;
    for i = 1:n_plot
        h_i  = samples.theta(plot_idx(i), idx_h)';
        vs_i = samples.theta(plot_idx(i), idx_vs)';
        d_i = [0; cumsum(h_i)];
        for j = 1:n_vs
            if j < n_vs
                plot([vs_i(j) vs_i(j)], [d_i(j) d_i(j+1)], 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.5);
                plot([vs_i(j) vs_i(j+1)], [d_i(j+1) d_i(j+1)], 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.5);
            else
                plot([vs_i(j) vs_i(j)], [d_i(end) max_depth], 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.5);
            end
        end
    end
end

for j = 1:n_vs
    if j < n_vs
        plot([vs_mean(j) vs_mean(j)], [depths(j) depths(j+1)], 'b-', 'LineWidth', 2.5);
        plot([vs_mean(j) vs_mean(j+1)], [depths(j+1) depths(j+1)], 'b-', 'LineWidth', 2.5);
    else
        plot([vs_mean(j) vs_mean(j)], [depths(end) max_depth], 'b-', 'LineWidth', 2.5);
    end
end

set(gca, 'YDir', 'reverse');
xlabel('Vs (km/s)');
ylabel('Depth (km)');
title(sprintf('Posterior Velocity Profile - %s [%s]', data.name, mode_str));
grid on;
xlim([0 5]);
ylim([0 max_depth]);

%% -------------------- Fit summary (unchanged) --------------------
fprintf('\nData fit summary:\n');
if isfield(pred_final, 'hvsr') && ~isempty(pred_final.hvsr) && isfield(data,'hvsr') && isfield(data.hvsr,'obs') ...
        && length(pred_final.hvsr) == length(data.hvsr.obs)
    hvsr_rms = sqrt(mean((pred_final.hvsr - data.hvsr.obs).^2));
    hvsr_chi2 = sum(((pred_final.hvsr - data.hvsr.obs)./data.hvsr.sig).^2) / length(data.hvsr.obs);
    fprintf('HVSR: RMS=%.3f, Reduced Chi2=%.2f\n', hvsr_rms, hvsr_chi2);
end
if isfield(pred_final, 'ellip') && ~isempty(pred_final.ellip) && isfield(data,'ellip') && isfield(data.ellip,'obs') ...
        && length(pred_final.ellip) == length(data.ellip.obs)
    ellip_rms = sqrt(mean((pred_final.ellip - data.ellip.obs).^2));
    ellip_chi2 = sum(((pred_final.ellip - data.ellip.obs)./data.ellip.sig).^2) / length(data.ellip.obs);
    fprintf('Ellipticity: RMS=%.3f, Reduced Chi2=%.2f\n', ellip_rms, ellip_chi2);
end
if isfield(pred_final, 'cph') && ~isempty(pred_final.cph) && isfield(data,'cph') && isfield(data.cph,'obs') ...
        && length(pred_final.cph) == length(data.cph.obs)
    cph_rms = sqrt(mean((pred_final.cph - data.cph.obs).^2));
    cph_chi2 = sum(((pred_final.cph - data.cph.obs)./data.cph.sig).^2) / length(data.cph.obs);
    fprintf('Phase Vel: RMS=%.3f km/s, Reduced Chi2=%.2f\n', cph_rms, cph_chi2);
end
if isfield(pred_final, 'ugr') && ~isempty(pred_final.ugr) && isfield(data,'ugr') && isfield(data.ugr,'obs') ...
        && length(pred_final.ugr) == length(data.ugr.obs)
    ugr_rms = sqrt(mean((pred_final.ugr - data.ugr.obs).^2));
    ugr_chi2 = sum(((pred_final.ugr - data.ugr.obs)./data.ugr.sig).^2) / length(data.ugr.obs);
    fprintf('Group Vel: RMS=%.3f km/s, Reduced Chi2=%.2f\n', ugr_rms, ugr_chi2);
end

%% -------------------- Save figures (paths unchanged) --------------------
fig_dir = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/Codes/Mod5_HBI/figures_new_structure';
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

station_name = strrep(data.name, '.', '_');
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

saveas(figure(1), fullfile(fig_dir, sprintf('%s_%s_traces_%s.png', station_name, mode_tag, timestamp)));
saveas(figure(2), fullfile(fig_dir, sprintf('%s_%s_obs_vs_pred_%s.png', station_name, mode_tag, timestamp)));
saveas(figure(3), fullfile(fig_dir, sprintf('%s_%s_vs_profile_%s.png', station_name, mode_tag, timestamp)));

fprintf('\nFigures saved to: %s\n', fig_dir);
fprintf('\nRun complete.\n');
