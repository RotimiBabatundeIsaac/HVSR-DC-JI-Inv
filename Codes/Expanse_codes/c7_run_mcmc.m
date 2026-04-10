function results = c7_run_mcmc(data, prior)

n_iter  = prior.mcmc.n_iterations;
burn_in = prior.mcmc.burn_in;
thin    = prior.mcmc.thin;
n_params = prior.n_params;
n_noise  = prior.noise.n_noise;

use_corr = isfield(prior.noise, 'use_correlated') && prior.noise.use_correlated;
if use_corr
    n_corrlen = n_noise;
    n_total_params = n_params + n_noise + n_corrlen;
else
    n_corrlen = 0;
    n_total_params = n_params + n_noise;
end

n_saved = floor((n_iter - burn_in) / thin);
if n_saved < 1
    error('c7_run_mcmc: No samples will be saved. Check burn_in and thin settings.');
end

samples_theta = zeros(n_saved, n_params);
samples_sigma = zeros(n_saved, n_noise);
samples_logL  = zeros(n_saved, 1);
if use_corr
    samples_corrlen = zeros(n_saved, n_noise);
end

n_hvsr  = 0; n_ellip = 0; n_cph = 0; n_ugr = 0;
if isfield(data, 'hvsr')  && isfield(data.hvsr, 'f'),  n_hvsr  = length(data.hvsr.f);  end
if isfield(data, 'ellip') && isfield(data.ellip, 'T'), n_ellip = length(data.ellip.T); end
if isfield(data, 'cph')   && isfield(data.cph, 'T'),   n_cph   = length(data.cph.T);   end
if isfield(data, 'ugr')   && isfield(data.ugr, 'T'),   n_ugr   = length(data.ugr.T);   end

samples_pred_hvsr  = NaN(n_saved, n_hvsr);
samples_pred_ellip = NaN(n_saved, n_ellip);
samples_pred_cph   = NaN(n_saved, n_cph);
samples_pred_ugr   = NaN(n_saved, n_ugr);

model   = c2_initialize_hbi_model(prior);
theta   = model.theta(:);
sigma_e = model.sigma_e(:);

if use_corr
    corrlen = prior.noise.corrlen_init(:);
    corrlen_step = prior.noise.corrlen_step(:);
    corrlen_lb = prior.noise.corrlen_lb(:);
    corrlen_ub = prior.noise.corrlen_ub(:);
else
    corrlen = [];
end

sigma_step = [0.05; 0.05; 0.05; 0.05];
sigma_lb = 0.01;
sigma_ub = 50.0;

[log_L, ~, pred, Phi, N_d, chol_cache] = c3_compute_likelihood(theta, sigma_e, data, prior, corrlen);
log_prior_sigma = compute_log_prior_sigma(sigma_e, prior);
if ~isfinite(log_L)
    log_L = -1e10;
end

accept_count_theta = zeros(n_params, 1);
total_count_theta  = zeros(n_params, 1);
accept_count_sigma = zeros(n_noise, 1);
total_count_sigma  = zeros(n_noise, 1);
if use_corr
    accept_count_corrlen = zeros(n_noise, 1);
    total_count_corrlen  = zeros(n_noise, 1);
end

best = struct();
best.theta = theta;
best.sigma = sigma_e;
best.logL  = log_L;
if use_corr
    best.corrlen = corrlen;
end

sample_idx = 0;

fprintf('[c7] Starting MCMC: %d iterations, burn_in=%d, thin=%d\n', n_iter, burn_in, thin);
if use_corr
    fprintf('[c7] Correlated noise: ON (Cholesky caching enabled)\n');
else
    fprintf('[c7] Correlated noise: OFF (diagonal covariance)\n');
end
fprintf('[c7] Delayed rejection: ON\n');
fprintf('[c7] Initial log_L = %.4f\n', log_L);

for iter = 1:n_iter
    j = randi(n_total_params);

    if j <= n_params
        [theta_prop, ~, is_feasible, ~] = c5_perturb_theta(theta, prior, j);
        total_count_theta(j) = total_count_theta(j) + 1;

        accepted = false;
        if is_feasible
            [log_L_prop, ~, pred_prop, Phi_prop, N_d_prop, ~] = c3_compute_likelihood(theta_prop, sigma_e, data, prior, corrlen, chol_cache);
            if isfinite(log_L_prop)
                log_alpha = log_L_prop - log_L;
                if log(rand()) < log_alpha
                    theta = theta_prop;
                    log_L = log_L_prop;
                    pred  = pred_prop;
                    Phi   = Phi_prop;
                    N_d   = N_d_prop;
                    accept_count_theta(j) = accept_count_theta(j) + 1;
                    accepted = true;
                end
            end
        end

        if ~accepted
            j2 = randi(n_params);
            if is_feasible
                base_theta = theta_prop;
            else
                base_theta = theta;
            end
            [theta_prop2, ~, is_feasible2, ~] = c5_perturb_theta(base_theta, prior, j2);
            if is_feasible2
                [log_L_prop2, ~, pred_prop2, Phi_prop2, N_d_prop2, ~] = c3_compute_likelihood(theta_prop2, sigma_e, data, prior, corrlen, chol_cache);
                if isfinite(log_L_prop2)
                    log_alpha2 = log_L_prop2 - log_L;
                    if log(rand()) < log_alpha2
                        theta = theta_prop2;
                        log_L = log_L_prop2;
                        pred  = pred_prop2;
                        Phi   = Phi_prop2;
                        N_d   = N_d_prop2;
                        accept_count_theta(j) = accept_count_theta(j) + 1;
                    end
                end
            end
        end

    elseif j <= n_params + n_noise
        d = j - n_params;
        total_count_sigma(d) = total_count_sigma(d) + 1;

        log_sigma_prop = log(sigma_e(d)) + sigma_step(d) * randn;
        sigma_prop = exp(log_sigma_prop);

        if sigma_prop >= sigma_lb && sigma_prop <= sigma_ub
            sigma_e_prop = sigma_e;
            sigma_e_prop(d) = sigma_prop;

            [log_L_prop, ~, pred_prop, Phi_prop, N_d_prop, chol_cache_prop] = c3_compute_likelihood(theta, sigma_e_prop, data, prior, corrlen);
            log_prior_sigma_prop = compute_log_prior_sigma(sigma_e_prop, prior);

            if isfinite(log_L_prop)
                log_alpha = (log_L_prop - log_L) ...
                          + (log_prior_sigma_prop - log_prior_sigma) ...
                          + (log_sigma_prop - log(sigma_e(d)));

                if log(rand()) < log_alpha
                    sigma_e = sigma_e_prop;
                    log_L = log_L_prop;
                    pred  = pred_prop;
                    Phi   = Phi_prop;
                    N_d   = N_d_prop;
                    log_prior_sigma = log_prior_sigma_prop;
                    chol_cache = chol_cache_prop;
                    accept_count_sigma(d) = accept_count_sigma(d) + 1;
                end
            end
        end

    else
        d = j - n_params - n_noise;
        total_count_corrlen(d) = total_count_corrlen(d) + 1;

        log_L_prop_val = log(corrlen(d)) + corrlen_step(d) * randn;
        corrlen_prop_val = exp(log_L_prop_val);

        if corrlen_prop_val >= corrlen_lb(d) && corrlen_prop_val <= corrlen_ub(d)
            corrlen_prop = corrlen;
            corrlen_prop(d) = corrlen_prop_val;

            [log_L_prop, ~, pred_prop, Phi_prop, N_d_prop, chol_cache_prop] = c3_compute_likelihood(theta, sigma_e, data, prior, corrlen_prop);

            if isfinite(log_L_prop)
                log_alpha = (log_L_prop - log_L) ...
                          + (log_L_prop_val - log(corrlen(d)));

                if log(rand()) < log_alpha
                    corrlen = corrlen_prop;
                    log_L = log_L_prop;
                    pred  = pred_prop;
                    Phi   = Phi_prop;
                    N_d   = N_d_prop;
                    chol_cache = chol_cache_prop;
                    accept_count_corrlen(d) = accept_count_corrlen(d) + 1;
                end
            end
        end
    end

    if log_L > best.logL
        best.theta = theta;
        best.sigma = sigma_e;
        best.logL  = log_L;
        if use_corr
            best.corrlen = corrlen;
        end
    end

    if iter > burn_in && mod(iter - burn_in, thin) == 0
        sample_idx = sample_idx + 1;
        samples_theta(sample_idx, :) = theta';
        samples_sigma(sample_idx, :) = sigma_e';
        samples_logL(sample_idx)     = log_L;
        if use_corr
            samples_corrlen(sample_idx, :) = corrlen';
        end

        if ~isempty(pred.hvsr) && n_hvsr > 0
            pv = pred.hvsr(:);
            nc = min(length(pv), n_hvsr);
            samples_pred_hvsr(sample_idx, 1:nc) = pv(1:nc)';
        end
        if ~isempty(pred.ellip) && n_ellip > 0
            pv = pred.ellip(:);
            nc = min(length(pv), n_ellip);
            samples_pred_ellip(sample_idx, 1:nc) = pv(1:nc)';
        end
        if ~isempty(pred.cph) && n_cph > 0
            pv = pred.cph(:);
            nc = min(length(pv), n_cph);
            samples_pred_cph(sample_idx, 1:nc) = pv(1:nc)';
        end
        if ~isempty(pred.ugr) && n_ugr > 0
            pv = pred.ugr(:);
            nc = min(length(pv), n_ugr);
            samples_pred_ugr(sample_idx, 1:nc) = pv(1:nc)';
        end
    end

    if mod(iter, 1000) == 0
        acc_theta = sum(accept_count_theta) / max(sum(total_count_theta), 1);
        acc_sigma = sum(accept_count_sigma) / max(sum(total_count_sigma), 1);
        sigma_str = sprintf('%.3f ', sigma_e);
        if use_corr
            acc_corrlen = sum(accept_count_corrlen) / max(sum(total_count_corrlen), 1);
            corrlen_str = sprintf('%.3f ', corrlen);
            fprintf('[c7] Iter %d/%d, logL=%.2f, acc_th=%.3f, acc_sig=%.3f, acc_L=%.3f, sigma=[%s], L=[%s]\n', ...
                iter, n_iter, log_L, acc_theta, acc_sigma, acc_corrlen, strtrim(sigma_str), strtrim(corrlen_str));
        else
            fprintf('[c7] Iter %d/%d, logL=%.2f, acc_theta=%.3f, acc_sigma=%.3f, sigma=[%s]\n', ...
                iter, n_iter, log_L, acc_theta, acc_sigma, strtrim(sigma_str));
        end
    end
end

samples_theta = samples_theta(1:sample_idx, :);
samples_sigma = samples_sigma(1:sample_idx, :);
samples_logL  = samples_logL(1:sample_idx);
samples_pred_hvsr  = samples_pred_hvsr(1:sample_idx, :);
samples_pred_ellip = samples_pred_ellip(1:sample_idx, :);
samples_pred_cph   = samples_pred_cph(1:sample_idx, :);
samples_pred_ugr   = samples_pred_ugr(1:sample_idx, :);

results.samples_theta = samples_theta;
results.samples_sigma = samples_sigma;
results.samples_logL  = samples_logL;
results.samples_pred_hvsr  = samples_pred_hvsr;
results.samples_pred_ellip = samples_pred_ellip;
results.samples_pred_cph   = samples_pred_cph;
results.samples_pred_ugr   = samples_pred_ugr;
results.accept_count    = accept_count_theta;
results.total_count     = total_count_theta;
results.acceptance_rate = accept_count_theta ./ max(total_count_theta, 1);
results.accept_count_sigma = accept_count_sigma;
results.total_count_sigma  = total_count_sigma;
results.acceptance_rate_sigma = accept_count_sigma ./ max(total_count_sigma, 1);
results.final_theta = theta;
results.final_sigma = sigma_e;
results.final_logL  = log_L;
results.best_theta = best.theta;
results.best_sigma = best.sigma;
results.best_logL  = best.logL;
results.prior     = prior;
results.n_samples = sample_idx;

if use_corr
    samples_corrlen = samples_corrlen(1:sample_idx, :);
    results.samples_corrlen = samples_corrlen;
    results.accept_count_corrlen = accept_count_corrlen;
    results.total_count_corrlen  = total_count_corrlen;
    results.acceptance_rate_corrlen = accept_count_corrlen ./ max(total_count_corrlen, 1);
    results.final_corrlen = corrlen;
    results.best_corrlen  = best.corrlen;
end

fprintf('[c7] MCMC complete. Samples saved: %d\n', sample_idx);
fprintf('[c7] Per-parameter acceptance rates (theta):\n');
for i = 1:n_params
    fprintf('  %s: %.3f (%d/%d)\n', prior.param_names{i}, results.acceptance_rate(i), ...
        accept_count_theta(i), total_count_theta(i));
end
fprintf('[c7] Per-noise acceptance rates (sigma):\n');
for d = 1:n_noise
    fprintf('  %s: %.3f (%d/%d)\n', prior.noise.names{d}, results.acceptance_rate_sigma(d), ...
        accept_count_sigma(d), total_count_sigma(d));
end
if use_corr
    fprintf('[c7] Per-noise acceptance rates (corrlen):\n');
    for d = 1:n_noise
        fprintf('  %s: %.3f (%d/%d)\n', prior.noise.corrlen_names{d}, results.acceptance_rate_corrlen(d), ...
            accept_count_corrlen(d), total_count_corrlen(d));
    end
end
fprintf('[c7] Best (max logL) state seen: logL=%.3f\n', results.best_logL);

end


function log_p = compute_log_prior_sigma(sigma_e, prior)

log_p = 0;
for d = 1:prior.noise.n_noise
    alpha = prior.noise.alpha_0(d);
    beta  = prior.noise.beta_0(d);
    s = sigma_e(d);
    if s <= 0
        log_p = -Inf;
        return;
    end
    tau = 1 / s^2;
    log_p = log_p + (alpha - 1) * log(tau) - beta * tau - 2 * log(s);
end

end
