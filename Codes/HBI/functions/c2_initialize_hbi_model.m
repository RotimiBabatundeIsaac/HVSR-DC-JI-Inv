function model = c2_initialize_hbi_model(prior)
% Initialize theta, mu_theta, Sigma_theta, sigma_e with valid starting values
%
% Input:
%   prior - structure from c1_define_hbi_prior
%
% Output:
%   model - structure containing initialized model state

model.theta = prior.mu_theta_init;

model.theta = enforce_monotonicity(model.theta);

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


function theta = enforce_monotonicity(theta)
vs_idx = 6:11;
vs = theta(vs_idx);
for i = 2:length(vs)
    if vs(i) < vs(i-1)
        vs(i) = vs(i-1) + 0.01;
    end
end
theta(vs_idx) = vs;
end


function is_valid = check_constraints(theta, prior)
is_valid = true;

if any(theta < prior.hard_min) || any(theta > prior.hard_max)
    is_valid = false;
    return;
end

vs = theta(6:11);
if any(diff(vs) < 0)
    is_valid = false;
    return;
end

h = theta(1:5);
if any(h <= 0)
    is_valid = false;
    return;
end

end