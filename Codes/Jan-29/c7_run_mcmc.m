function results = c7_run_mcmc(data, prior)

n_iter  = prior.mcmc.n_iterations;
burn_in = prior.mcmc.burn_in;
thin    = prior.mcmc.thin;
n_params = prior.n_params;
n_noise  = prior.noise.n_noise;

n_saved = floor((n_iter - burn_in) / thin);
if n_saved < 1
    error('c7_run_mcmc: No samples will be saved. Check burn_in and thin settings.');
end

samples_theta = zeros(n_saved, n_params);
samples_sigma = zeros(n_saved, n_noise);
samples_logL  = zeros(n_saved, 1);

model   = c2_initialize_hbi_model(prior);
theta   = model.theta(:);
sigma_e = model.sigma_e(:);

[log_L, ~, pred, Phi, N_d] = c3_compute_likelihood(theta, sigma_e, data, prior);
if ~isfinite(log_L)
    log_L = -1e10;
end

accept_count = zeros(n_params, 1);
total_count  = zeros(n_params, 1);

best = struct();
best.theta = theta;
best.sigma = sigma_e;
best.logL  = log_L;

sample_idx = 0;

fprintf('[c7] Starting MCMC: %d iterations, burn_in=%d, thin=%d\n', n_iter, burn_in, thin);
fprintf('[c7] Initial log_L = %.4f\n', log_L);

for iter = 1:n_iter
    j = randi(n_params);
    [theta_prop, ~, is_feasible, ~] = c5_perturb_theta(theta, prior, j);
    total_count(j) = total_count(j) + 1;
    
    if is_feasible
        [log_L_prop, ~, pred_prop, Phi_prop, N_d_prop] = c3_compute_likelihood(theta_prop, sigma_e, data, prior);
        if isfinite(log_L_prop)
            log_alpha = log_L_prop - log_L;
            if log(rand()) < log_alpha
                theta = theta_prop;
                log_L = log_L_prop;
                pred  = pred_prop;
                Phi   = Phi_prop;
                N_d   = N_d_prop;
                accept_count(j) = accept_count(j) + 1;
            end
        end
    end
    
    sigma_e = c6_gibbs_update_sigma(Phi, N_d, sigma_e, prior);
    [log_L, ~, pred, Phi, N_d] = c3_compute_likelihood(theta, sigma_e, data, prior);
    if ~isfinite(log_L)
        log_L = -1e10;
    end
    
    if log_L > best.logL
        best.theta = theta;
        best.sigma = sigma_e;
        best.logL  = log_L;
    end
    
    if iter > burn_in && mod(iter - burn_in, thin) == 0
        sample_idx = sample_idx + 1;
        samples_theta(sample_idx, :) = theta';
        samples_sigma(sample_idx, :) = sigma_e';
        samples_logL(sample_idx)     = log_L;
    end
    
    if mod(iter, 1000) == 0
        acc_rate = sum(accept_count) / max(sum(total_count), 1);
        sigma_str = sprintf('%.3f ', sigma_e);
        fprintf('[c7] Iter %d/%d, logL=%.2f, acc=%.3f, last_j=%d (%s), sigma=[%s]\n', ...
            iter, n_iter, log_L, acc_rate, j, prior.param_names{j}, strtrim(sigma_str));
    end
end

samples_theta = samples_theta(1:sample_idx, :);
samples_sigma = samples_sigma(1:sample_idx, :);
samples_logL  = samples_logL(1:sample_idx);

results.samples_theta = samples_theta;
results.samples_sigma = samples_sigma;
results.samples_logL  = samples_logL;
results.accept_count    = accept_count;
results.total_count     = total_count;
results.acceptance_rate = accept_count ./ max(total_count, 1);
results.final_theta = theta;
results.final_sigma = sigma_e;
results.final_logL  = log_L;
results.best_theta = best.theta;
results.best_sigma = best.sigma;
results.best_logL  = best.logL;
results.prior     = prior;
results.n_samples = sample_idx;

fprintf('[c7] MCMC complete. Samples saved: %d\n', sample_idx);
fprintf('[c7] Per-parameter acceptance rates:\n');
for i = 1:n_params
    fprintf('  %s: %.3f (%d/%d)\n', prior.param_names{i}, results.acceptance_rate(i), ...
        accept_count(i), total_count(i));
end
fprintf('[c7] Best (max logL) state seen: logL=%.3f\n', results.best_logL);

end