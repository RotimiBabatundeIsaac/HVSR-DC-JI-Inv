function results = c8_run_multi_chain(data, prior)

n_chains = prior.mcmc.n_chains;
n_params = prior.n_params;
n_noise  = prior.noise.n_noise;

if n_chains < 2
    fprintf('[c8] Running single chain\n');
    results = c7_run_mcmc(data, prior);
    results.Rhat = NaN(n_params, 1);
    return;
end

fprintf('[c8] Running %d chains\n', n_chains);

chain_results = cell(n_chains, 1);

if isfield(prior, 'mcmc') && isfield(prior.mcmc, 'base_seed') && ~isempty(prior.mcmc.base_seed)
    base_seed = prior.mcmc.base_seed;
else
    base_seed = sum(100*clock);
end

streamType = 'Threefry';

use_parallel = ~isempty(ver('parallel')) && ~isempty(gcp('nocreate'));

if use_parallel
    parfor c = 1:n_chains
        s = RandStream(streamType, 'Seed', base_seed);
        s.Substream = c;
        RandStream.setGlobalStream(s);

        fprintf('[c8] Starting chain %d (parallel)\n', c);
        chain_results{c} = c7_run_mcmc(data, prior);
    end
else
    for c = 1:n_chains
        s = RandStream(streamType, 'Seed', base_seed);
        s.Substream = c;
        RandStream.setGlobalStream(s);

        fprintf('[c8] Starting chain %d (serial)\n', c);
        chain_results{c} = c7_run_mcmc(data, prior);
    end
end

n_samples_chain = zeros(n_chains, 1);
for c = 1:n_chains
    n_samples_chain(c) = chain_results{c}.n_samples;
end

if any(n_samples_chain == 0)
    error('c8_run_multi_chain: One or more chains produced no samples');
end

min_samples = min(n_samples_chain);
n_total = min_samples * n_chains;

all_theta = zeros(n_total, n_params);
all_sigma = zeros(n_total, n_noise);
all_logL  = zeros(n_total, 1);

for c = 1:n_chains
    idx_start = (c-1) * min_samples + 1;
    idx_end   = c * min_samples;

    all_theta(idx_start:idx_end, :) = chain_results{c}.samples_theta(1:min_samples, :);
    all_sigma(idx_start:idx_end, :) = chain_results{c}.samples_sigma(1:min_samples, :);
    all_logL(idx_start:idx_end)     = chain_results{c}.samples_logL(1:min_samples);
end

results.all_theta = all_theta;
results.all_sigma = all_sigma;
results.all_logL  = all_logL;

results.chain_results = chain_results;
results.prior = prior;

results.n_chains = n_chains;
results.n_samples_per_chain = min_samples;
results.n_total = n_total;

results.theta_mean   = mean(all_theta, 1)';
results.theta_std    = std(all_theta, 0, 1)';
results.theta_median = median(all_theta, 1)';
results.theta_q05    = prctile(all_theta, 5, 1)';
results.theta_q95    = prctile(all_theta, 95, 1)';

results.sigma_mean = mean(all_sigma, 1)';
results.sigma_std  = std(all_sigma, 0, 1)';

[results.map_logL, map_idx] = max(all_logL);
results.theta_map = all_theta(map_idx, :)';
results.sigma_map = all_sigma(map_idx, :)';

results.Rhat = compute_gelman_rubin(chain_results, n_params, min_samples);

fprintf('\n[c8] Multi-chain results:\n');
fprintf('  Base seed: %g\n', base_seed);
fprintf('  Chains: %d, Samples per chain: %d, Total: %d\n', n_chains, min_samples, n_total);
fprintf('  MAP log-likelihood: %.2f\n', results.map_logL);

fprintf('\n[c8] Parameter estimates (mean +/- std):\n');
for i = 1:n_params
    fprintf('  %s: %.4f +/- %.4f %s (Rhat=%.3f)\n', ...
        prior.param_names{i}, results.theta_mean(i), results.theta_std(i), ...
        prior.param_units{i}, results.Rhat(i));
end

fprintf('\n[c8] Noise scale estimates (mean +/- std):\n');
for d = 1:n_noise
    fprintf('  %s: %.4f +/- %.4f\n', prior.noise.names{d}, results.sigma_mean(d), results.sigma_std(d));
end

n_converged = sum(results.Rhat < 1.1);
fprintf('\n[c8] Convergence: %d/%d parameters have Rhat < 1.1\n', n_converged, n_params);

end


function Rhat = compute_gelman_rubin(chain_results, n_params, n_samples)

n_chains = length(chain_results);
Rhat = zeros(n_params, 1);

if n_chains < 2 || n_samples < 2
    Rhat(:) = NaN;
    return;
end

for p = 1:n_params
    chain_means = zeros(n_chains, 1);
    chain_vars  = zeros(n_chains, 1);

    for c = 1:n_chains
        samples = chain_results{c}.samples_theta(1:n_samples, p);
        chain_means(c) = mean(samples);
        chain_vars(c)  = var(samples, 0);
    end

    grand_mean = mean(chain_means);
    B = n_samples * sum((chain_means - grand_mean).^2) / (n_chains - 1);
    W = mean(chain_vars);

    if W > 1e-10
        var_plus = ((n_samples - 1) / n_samples) * W + (1 / n_samples) * B;
        Rhat(p) = sqrt(var_plus / W);
    else
        Rhat(p) = NaN;
    end
end

end