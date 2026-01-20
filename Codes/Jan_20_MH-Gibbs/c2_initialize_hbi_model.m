function model = c2_initialize_hbi_model(prior)
% Initialize theta, mu_theta, Sigma_theta, sigma_e with valid starting values
%
% Input:
%   prior - structure from c1_define_hbi_prior
%
% Output:
%   model - structure containing initialized model state

model.theta = prior.mu_theta_init;
model.theta = enforce_monotonicity(model.theta, prior);

if prior.use_hierarchical
    model.mu_theta = model.theta;
    model.Sigma_theta = prior.Psi_0;
else
    model.mu_theta = [];
    model.Sigma_theta = [];
end

model.sigma_e = prior.noise.init;
model.n_params = prior.n_params;
model.n_noise = prior.noise.n_noise;

model.log_likelihood = -Inf;
model.log_prior_theta = -Inf;
model.log_posterior = -Inf;

model.accept_count = 0;
model.reject_count = 0;

model.is_valid = check_constraints(model.theta, prior);

if ~model.is_valid
    warning('Initial theta violates constraints. Check prior bounds.');
end

end

function theta = enforce_monotonicity(theta, prior)
% Enforce monotonicity for all velocities (Vs1-Vs20)

if ~prior.enforce_monotonicity
    return;
end

mono_start = prior.monotonic_start_idx;
vs_mono = theta(mono_start:end);

for i = 2:length(vs_mono)
    if vs_mono(i) < vs_mono(i-1)
        vs_mono(i) = vs_mono(i-1) + 0.01;
    end
end

theta(mono_start:end) = vs_mono;

end

function is_valid = check_constraints(theta, prior)
% Check all constraints on theta

is_valid = true;

% Hard bounds
if any(theta < prior.hard_min) || any(theta > prior.hard_max)
    is_valid = false;
    return;
end

% Monotonicity for all velocities (Vs1-Vs20)
if prior.enforce_monotonicity
    mono_start = prior.monotonic_start_idx;
    vs_mono = theta(mono_start:end);
    if any(diff(vs_mono) < 0)
        is_valid = false;
        return;
    end
end

% Positive thicknesses
h = theta(prior.idx.h);
if any(h <= 0)
    is_valid = false;
    return;
end

end