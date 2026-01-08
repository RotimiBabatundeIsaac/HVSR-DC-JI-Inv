% e1_run_hbi_mcmc.m
% Main MCMC loop for Hierarchical Bayesian Inversion on station WD20
%
% Gibbs-within-Metropolis algorithm:
%   Step 1: Update theta via Metropolis-Hastings
%   Step 2: Update mu_theta via Gibbs
%   Step 3: Update Sigma_theta via Gibbs
%   Step 4: Update sigma_e via Gibbs

clear; clc; close all;

% Add paths
base_path = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/Codes';
addpath(genpath(base_path));

fprintf('Hierarchical Bayesian Inversion - Station D14\n');
fprintf('==============================================\n\n');

%% Load prior and initialize model
prior = c1_define_hbi_prior();
model = c2_initialize_hbi_model(prior);

%% Load station data
datapack_path = '/Users/birotimi/Library/CloudStorage/OneDrive-SyracuseUniversity/Desktop/PHD/HVSR-project/HVSR-Joint-Inversion/DATA/sta_datapack.mat';
fprintf('Loading station data...\n');
load(datapack_path, 'sta');
idx_w21 = find(contains({sta.name}, 'D14', 'IgnoreCase', true));
data = sta(idx_w01(1));
fprintf('Station: %s (lat=%.4f, lon=%.4f)\n\n', data.name, data.lat, data.lon);

%% MCMC settings
n_iterations = prior.mcmc.n_iterations;
burn_in = prior.mcmc.burn_in;
thin = prior.mcmc.thin;
n_samples = floor((n_iterations - burn_in) / thin);

fprintf('MCMC Settings:\n');
fprintf('  Iterations: %d\n', n_iterations);
fprintf('  Burn-in: %d\n', burn_in);
fprintf('  Thinning: %d\n', thin);
fprintf('  Output samples: %d\n\n', n_samples);

%% Preallocate storage
samples.theta = zeros(n_samples, prior.n_params);
samples.mu_theta = zeros(n_samples, prior.n_params);
samples.Sigma_theta = zeros(n_samples, prior.n_params, prior.n_params);
samples.sigma_e = zeros(n_samples, prior.noise.n_noise);
samples.log_likelihood = zeros(n_samples, 1);
samples.log_posterior = zeros(n_samples, 1);

% Traces for diagnostics (store every iteration)
trace.log_likelihood = zeros(n_iterations, 1);
trace.log_posterior = zeros(n_iterations, 1);
trace.accept = zeros(n_iterations, 1);

%% Initialize current state
theta = model.theta;
mu_theta = model.mu_theta;
Sigma_theta = model.Sigma_theta;
sigma_e = model.sigma_e;

% Compute initial likelihood and prior
fprintf('Computing initial likelihood...\n');
[log_L, ~, pred] = c3_compute_likelihood(theta, sigma_e, data, prior);
[log_prior, ~] = c4_compute_log_prior_theta(theta, mu_theta, Sigma_theta, prior);
log_posterior = log_L + log_prior;

fprintf('Initial log-likelihood: %.2f\n', log_L);
fprintf('Initial log-prior: %.2f\n', log_prior);
fprintf('Initial log-posterior: %.2f\n\n', log_posterior);

%% MCMC loop
fprintf('Starting MCMC...\n');
fprintf('Progress: ');

accept_count = 0;
sample_idx = 0;
tic;

for iter = 1:n_iterations
    
    % Step 1: Update theta via Metropolis-Hastings
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
    
    % Step 2: Update mu_theta via Gibbs
    mu_theta = c6_gibbs_update_mu_theta(theta, Sigma_theta, prior);
    
    % Step 3: Update Sigma_theta via Gibbs
    Sigma_theta = c7_gibbs_update_Sigma_theta(theta, mu_theta, prior);
    
    % Step 4: Update sigma_e via Gibbs
    sigma_e = c8_gibbs_update_sigma_e(pred, data, prior);
    
    % Recompute log_prior with updated hyperparameters
    [log_prior, ~] = c4_compute_log_prior_theta(theta, mu_theta, Sigma_theta, prior);
    log_posterior = log_L + log_prior;
    
    % Store trace
    trace.log_likelihood(iter) = log_L;
    trace.log_posterior(iter) = log_posterior;
    
    % Store samples after burn-in with thinning
    if iter > burn_in && mod(iter - burn_in, thin) == 0
        sample_idx = sample_idx + 1;
        samples.theta(sample_idx, :) = theta';
        samples.mu_theta(sample_idx, :) = mu_theta';
        samples.Sigma_theta(sample_idx, :, :) = Sigma_theta;
        samples.sigma_e(sample_idx, :) = sigma_e';
        samples.log_likelihood(sample_idx) = log_L;
        samples.log_posterior(sample_idx) = log_posterior;
    end
    
    % Progress display
    if mod(iter, floor(n_iterations/10)) == 0
        fprintf('%d%% ', round(100*iter/n_iterations));
    end
end

elapsed_time = toc;
fprintf('\nDone!\n\n');

%% Summary statistics
fprintf('MCMC Summary:\n');
fprintf('  Total time: %.1f minutes\n', elapsed_time/60);
fprintf('  Time per iteration: %.3f s\n', elapsed_time/n_iterations);
fprintf('  Acceptance rate: %.1f%%\n\n', 100*accept_count/n_iterations);

fprintf('Parameter estimates (posterior mean +/- std):\n');
fprintf('  %-6s  %10s  %10s  %10s\n', 'Param', 'Mean', 'Std', 'Initial');
for i = 1:prior.n_params
    fprintf('  %-6s  %10.4f  %10.4f  %10.4f\n', prior.param_names{i}, ...
        mean(samples.theta(:,i)), std(samples.theta(:,i)), model.theta(i));
end

fprintf('\nNoise estimates:\n');
for i = 1:prior.noise.n_noise
    fprintf('  %-12s: %.4f +/- %.4f\n', prior.noise.names{i}, ...
        mean(samples.sigma_e(:,i)), std(samples.sigma_e(:,i)));
end

%% Save results
output_file = 'hbi_results_W01.mat';
save(output_file, 'samples', 'trace', 'prior', 'data', 'model', 'elapsed_time', 'accept_count');
fprintf('\nResults saved to: %s\n', output_file);

%% Quick diagnostic plots
figure('Position', [50, 50, 1400, 900]);

% Plot 1: Log-likelihood trace
subplot(2,3,1);
plot(trace.log_likelihood, 'b-', 'LineWidth', 0.5);
hold on;
xline(burn_in, 'r--', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Log-Likelihood');
title('Log-Likelihood Trace');
grid on;

% Plot 2: Log-posterior trace
subplot(2,3,2);
plot(trace.log_posterior, 'b-', 'LineWidth', 0.5);
hold on;
xline(burn_in, 'r--', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Log-Posterior');
title('Log-Posterior Trace');
grid on;

% Plot 3: Acceptance rate (running)
subplot(2,3,3);
window = 1000;
running_accept = movmean(trace.accept, window) * 100;
plot(running_accept, 'b-', 'LineWidth', 1);
hold on;
yline(23.4, 'r--', 'LineWidth', 1.5);
xline(burn_in, 'k--', 'LineWidth', 1);
xlabel('Iteration');
ylabel('Acceptance Rate (%)');
title(sprintf('Running Acceptance (window=%d)', window));
legend('Acceptance', 'Optimal (23.4%)', 'Burn-in', 'Location', 'best');
grid on;

% Plot 4: Theta traces (selected parameters)
subplot(2,3,4);
plot(samples.theta(:,1), 'LineWidth', 0.5); hold on;
plot(samples.theta(:,6), 'LineWidth', 0.5);
plot(samples.theta(:,9), 'LineWidth', 0.5);
xlabel('Sample');
ylabel('Value');
title('Parameter Traces');
legend('h1', 'Vs1', 'Vs4', 'Location', 'best');
grid on;

% Plot 5: Velocity profile with uncertainty
subplot(2,3,5);
h_mean = mean(samples.theta(:,1:5), 1);
vs_mean = mean(samples.theta(:,6:11), 1);
h_std = std(samples.theta(:,1:5), 0, 1);
vs_std = std(samples.theta(:,6:11), 0, 1);

depths = [0, cumsum(h_mean)];
max_depth = sum(h_mean) + 10;

hold on;
for i = 1:6
    if i <= 5
        d_top = depths(i);
        d_bot = depths(i+1);
    else
        d_top = depths(6);
        d_bot = max_depth;
    end
    fill([vs_mean(i)-vs_std(i), vs_mean(i)+vs_std(i), vs_mean(i)+vs_std(i), vs_mean(i)-vs_std(i)], ...
         [d_top, d_top, d_bot, d_bot], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    plot([vs_mean(i), vs_mean(i)], [d_top, d_bot], 'b-', 'LineWidth', 2);
end
set(gca, 'YDir', 'reverse');
xlabel('Vs (km/s)');
ylabel('Depth (km)');
title('Posterior Velocity Profile');
xlim([0 5]);
grid on;

% Plot 6: Noise parameter histograms
subplot(2,3,6);
for i = 1:4
    subplot(2,3,6);
    hold on;
end
histogram(samples.sigma_e(:,1), 20, 'FaceAlpha', 0.5, 'DisplayName', 'HVSR');
histogram(samples.sigma_e(:,2), 20, 'FaceAlpha', 0.5, 'DisplayName', 'Ellip');
histogram(samples.sigma_e(:,3), 20, 'FaceAlpha', 0.5, 'DisplayName', 'Phase');
histogram(samples.sigma_e(:,4), 20, 'FaceAlpha', 0.5, 'DisplayName', 'Group');
xlabel('sigma\_e');
ylabel('Count');
title('Noise Parameter Distributions');
legend('Location', 'best');
grid on;

sgtitle(sprintf('HBI MCMC Results - %s (%.1f min, %.1f%% accept)', ...
    data.name, elapsed_time/60, 100*accept_count/n_iterations));

saveas(gcf, 'hbi_mcmc_diagnostics_D14.png');
fprintf('Diagnostic plot saved to: hbi_mcmc_diagnostics_D14.png\n');

fprintf('\nHBI MCMC Complete!\n');