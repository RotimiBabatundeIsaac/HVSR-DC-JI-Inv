%% test_e1_run_hbi_mcmc.m
clear; clc; close all;

base_path = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/Codes';
addpath(genpath(base_path));

fprintf('MCMC Test Run - Station E13\n');
fprintf('================================\n\n');

%% Load prior and initialize model
prior = c1_define_hbi_prior();
model = c2_initialize_hbi_model(prior);

%% Display mode
fprintf('========================================\n');
if prior.use_hierarchical
    fprintf('MODE: HIERARCHICAL BAYESIAN (Gaussian prior)\n');
    mode_str = 'HBI Gaussian';
    mode_tag = 'HBI';
else
    fprintf('MODE: SIMPLE MCMC (Uniform prior, fixed sigma_e)\n');
    mode_str = 'Uniform Prior';
    mode_tag = 'Uniform';
end
fprintf('========================================\n\n');

%% Override MCMC settings for quick test, 
% remember to set this with the init BIR 2026 
n_iterations = 1000;
burn_in = 40;
thin = 2;
n_samples = floor((n_iterations - burn_in) / thin);

fprintf('Test MCMC Settings:\n');
fprintf('  Iterations: %d\n', n_iterations);
fprintf('  Burn-in: %d\n', burn_in);
fprintf('  Thinning: %d\n', thin);
fprintf('  Output samples: %d\n\n', n_samples);

%% Load station data
datapack_path = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/DATA/sta_datapack.mat';
load(datapack_path, 'sta');
idx_E13 = find(contains({sta.name}, 'E13', 'IgnoreCase', true));
data = sta(idx_E13(1));
fprintf('Station: %s\n\n', data.name);

%% Preallocate
samples.theta = zeros(n_samples, prior.n_params);
samples.sigma_e = zeros(n_samples, prior.noise.n_noise);
samples.log_likelihood = zeros(n_samples, 1);

trace.log_likelihood = zeros(n_iterations, 1);
trace.accept = zeros(n_iterations, 1);

%% Initialize
theta = model.theta;
sigma_e = model.sigma_e;

if prior.use_hierarchical
    mu_theta = model.mu_theta;
    Sigma_theta = model.Sigma_theta;
else
    mu_theta = [];
    Sigma_theta = [];
end

fprintf('Computing initial likelihood...\n');
[log_L, ~, pred] = c3_compute_likelihood(theta, sigma_e, data, prior);
[log_prior, ~] = c4_compute_log_prior_theta(theta, mu_theta, Sigma_theta, prior);
log_posterior = log_L + log_prior;

fprintf('Initial log-likelihood: %.2f\n', log_L);
fprintf('Initial log-prior: %.2f\n', log_prior);
fprintf('Initial log-posterior: %.2f\n\n', log_posterior);

%% MCMC loop
fprintf('Running MCMC test (%d iterations)...\n', n_iterations);
accept_count = 0;
sample_idx = 0;
tic;

for iter = 1:n_iterations
    [theta_prop, is_valid, ~] = c5_perturb_theta(theta, prior);
    
    if is_valid
        [log_L_prop, ~, pred_prop] = c3_compute_likelihood(theta_prop, sigma_e, data, prior);
        [log_prior_prop, ~] = c4_compute_log_prior_theta(theta_prop, mu_theta, Sigma_theta, prior);
        log_posterior_prop = log_L_prop + log_prior_prop;
        
        log_alpha = log_posterior_prop - log_posterior;
        
        if log(rand()) < log_alpha
            theta = theta_prop;
            log_L = log_L_prop;
            log_prior = log_prior_prop;
            log_posterior = log_posterior_prop;
            pred = pred_prop;
            accept_count = accept_count + 1;
            trace.accept(iter) = 1;
        end
    end
    
    if prior.use_hierarchical
        mu_theta = c6_gibbs_update_mu_theta(theta, Sigma_theta, prior);
        Sigma_theta = c7_gibbs_update_Sigma_theta(theta, mu_theta, prior);
        sigma_e = c8_gibbs_update_sigma_e(pred, data, prior);
        [log_prior, ~] = c4_compute_log_prior_theta(theta, mu_theta, Sigma_theta, prior);
        log_posterior = log_L + log_prior;
    end
    
    trace.log_likelihood(iter) = log_L;
    
    if iter > burn_in && mod(iter - burn_in, thin) == 0
        sample_idx = sample_idx + 1;
        samples.theta(sample_idx, :) = theta';
        samples.sigma_e(sample_idx, :) = sigma_e';
        samples.log_likelihood(sample_idx) = log_L;
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

%% Diagnostic plots
figure('Position', [100, 100, 1200, 400]);

subplot(1,3,1);
plot(trace.log_likelihood, 'b-');
xlabel('Iteration');
ylabel('Log-Likelihood');
title('Log-Likelihood Trace');
xline(burn_in, 'r--');
grid on;

subplot(1,3,2);
plot(samples.theta(:,6), 'b-');
hold on;
plot(samples.theta(:,9), 'r-');
xlabel('Sample');
ylabel('Velocity (km/s)');
title('Velocity Traces');
legend('Vs1', 'Vs4');
grid on;

subplot(1,3,3);
h_mean = mean(samples.theta(:,1:5), 1);
vs_mean = mean(samples.theta(:,6:11), 1);
depths = [0, cumsum(h_mean)];
max_depth = sum(h_mean) + 10;
hold on;
for i = 1:6
    if i <= 5
        plot([vs_mean(i), vs_mean(i)], [depths(i), depths(i+1)], 'b-', 'LineWidth', 2);
        if i < 6
            plot([vs_mean(i), vs_mean(i+1)], [depths(i+1), depths(i+1)], 'b-', 'LineWidth', 2);
        end
    else
        plot([vs_mean(i), vs_mean(i)], [depths(6), max_depth], 'b-', 'LineWidth', 2);
    end
end
set(gca, 'YDir', 'reverse');
xlabel('Vs (km/s)');
ylabel('Depth (km)');
title('Mean Velocity Profile');
xlim([0 5]);
grid on;

sgtitle(sprintf('%s (%.1f%% accept, %.1f s)', mode_str, 100*accept_count/n_iterations, elapsed));

%% PREDICTED VS OBSERVED
fprintf('\n========== PREDICTED VS OBSERVED ==========\n');

theta_mean = mean(samples.theta, 1)';
sigma_e_mean = mean(samples.sigma_e, 1)';

fprintf('\nPosterior mean model:\n');
fprintf('Layer thicknesses (km): %.3f %.3f %.3f %.3f %.3f\n', theta_mean(1:5));
fprintf('Vs values (km/s): %.3f %.3f %.3f %.3f %.3f %.3f\n', theta_mean(6:11));
fprintf('Noise std (sigma_e): %.3f %.3f %.3f %.3f\n', sigma_e_mean);

[~, ~, pred_final] = c3_compute_likelihood(theta_mean, sigma_e_mean, data, prior);

figure('Position', [100 100 1400 900], 'Name', 'Predicted vs Observed Data Fit');

subplot(2,2,1);
if isfield(data, 'hvsr') && ~isempty(data.hvsr.f) && ~isempty(pred_final.hvsr)
    errorbar(data.hvsr.f, data.hvsr.obs, data.hvsr.sig, ...
        'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', 1);
    hold on;
    plot(data.hvsr.f, pred_final.hvsr, 'b-', 'LineWidth', 2);
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
if isfield(data, 'ellip') && ~isempty(data.ellip.T) && ~isempty(pred_final.ellip)
    errorbar(data.ellip.T, data.ellip.obs, data.ellip.sig, ...
        'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', 1);
    hold on;
    plot(data.ellip.T, pred_final.ellip, 'r-', 'LineWidth', 2);
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
if isfield(data, 'cph') && ~isempty(data.cph.T) && ~isempty(pred_final.cph)
    errorbar(data.cph.T, data.cph.obs, data.cph.sig, ...
        'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', 1);
    hold on;
    plot(data.cph.T, pred_final.cph, 'g-', 'LineWidth', 2);
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
if isfield(data, 'ugr') && ~isempty(data.ugr.T) && ~isempty(pred_final.ugr)
    errorbar(data.ugr.T, data.ugr.obs, data.ugr.sig, ...
        'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'LineWidth', 1);
    hold on;
    plot(data.ugr.T, pred_final.ugr, 'm-', 'LineWidth', 2);
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

sgtitle(sprintf('Predicted vs Observed [%s]', mode_str), 'FontSize', 14, 'FontWeight', 'bold');

%% Data Fit Statistics
fprintf('\nData fit summary:\n');
if ~isempty(pred_final.hvsr)
    hvsr_rms = sqrt(mean((pred_final.hvsr - data.hvsr.obs).^2));
    hvsr_chi2 = sum(((pred_final.hvsr - data.hvsr.obs)./data.hvsr.sig).^2) / length(data.hvsr.obs);
    fprintf('HVSR: RMS=%.3f, Reduced Chi2=%.2f\n', hvsr_rms, hvsr_chi2);
end
if ~isempty(pred_final.ellip)
    ellip_rms = sqrt(mean((pred_final.ellip - data.ellip.obs).^2));
    ellip_chi2 = sum(((pred_final.ellip - data.ellip.obs)./data.ellip.sig).^2) / length(data.ellip.obs);
    fprintf('Ellipticity: RMS=%.3f, Reduced Chi2=%.2f\n', ellip_rms, ellip_chi2);
end
if ~isempty(pred_final.cph)
    cph_rms = sqrt(mean((pred_final.cph - data.cph.obs).^2));
    cph_chi2 = sum(((pred_final.cph - data.cph.obs)./data.cph.sig).^2) / length(data.cph.obs);
    fprintf('Phase Vel: RMS=%.3f km/s, Reduced Chi2=%.2f\n', cph_rms, cph_chi2);
end
if ~isempty(pred_final.ugr)
    ugr_rms = sqrt(mean((pred_final.ugr - data.ugr.obs).^2));
    ugr_chi2 = sum(((pred_final.ugr - data.ugr.obs)./data.ugr.sig).^2) / length(data.ugr.obs);
    fprintf('Group Vel: RMS=%.3f km/s, Reduced Chi2=%.2f\n', ugr_rms, ugr_chi2);
end

%% Velocity Profile with Uncertainty
figure('Position', [100 100 500 600], 'Name', 'Velocity Profile with Uncertainty');

h_mean = mean(samples.theta(:,1:5), 1);
vs_mean = mean(samples.theta(:,6:11), 1);

n_plot = min(200, size(samples.theta, 1));
plot_idx = randperm(size(samples.theta, 1), n_plot);

hold on;
for i = 1:n_plot
    h_i = samples.theta(plot_idx(i), 1:5);
    vs_i = samples.theta(plot_idx(i), 6:11);
    d_i = [0, cumsum(h_i)];
    for j = 1:6
        if j <= 5
            plot([vs_i(j) vs_i(j)], [d_i(j) d_i(j+1)], 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.5);
            plot([vs_i(j) vs_i(j+1)], [d_i(j+1) d_i(j+1)], 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.5);
        else
            plot([vs_i(j) vs_i(j)], [d_i(6) d_i(6)+10], 'Color', [0.7 0.7 0.7 0.3], 'LineWidth', 0.5);
        end
    end
end

depths_final = [0, cumsum(h_mean)];
for j = 1:6
    if j <= 5
        plot([vs_mean(j) vs_mean(j)], [depths_final(j) depths_final(j+1)], 'b-', 'LineWidth', 2);
        plot([vs_mean(j) vs_mean(j+1)], [depths_final(j+1) depths_final(j+1)], 'b-', 'LineWidth', 2);
    else
        plot([vs_mean(j) vs_mean(j)], [depths_final(6) depths_final(6)+10], 'b-', 'LineWidth', 2);
    end
end

set(gca, 'YDir', 'reverse');
xlabel('Vs (km/s)');
ylabel('Depth (km)');
title(sprintf('Posterior Velocity Profile [%s]', mode_str));
grid on;
xlim([0 5]);

%% Save Figures
fig_dir = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/Codes/Mod5_HBI/figures';
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

station_name = strrep(data.name, '.', '_');
timestamp = datestr(now, 'yyyymmdd_HHMMSS');

figs = findobj('Type', 'figure');
for i = 1:length(figs)
    fig = figs(i);
    fig_name = get(fig, 'Name');
    if isempty(fig_name)
        fig_name = sprintf('diagnostic_%d', fig.Number);
    else
        fig_name = strrep(fig_name, ' ', '_');
        fig_name = strrep(fig_name, '/', '_');
    end
    
    pdf_name = sprintf('%s_%s_%s_%s.pdf', station_name, mode_tag, fig_name, timestamp);
    pdf_path = fullfile(fig_dir, pdf_name);
    
    exportgraphics(fig, pdf_path, 'ContentType', 'vector', 'Resolution', 300);
    fprintf('Saved: %s\n', pdf_name);
end

fprintf('\nAll figures saved to: %s\n', fig_dir);
fprintf('\n============================================\n');
fprintf('Run complete.\n');