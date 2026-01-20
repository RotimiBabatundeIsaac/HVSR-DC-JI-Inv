function mu_theta_new = c6_gibbs_update_mu_theta(theta, Sigma_theta, prior)
% c6_gibbs_update_mu_theta
% Numerically stable Gibbs update for mu_theta:
%   theta | mu_theta, Sigma_theta ~ N(mu_theta, Sigma_theta)
%   mu_theta ~ N(mu0, Sigma0)
%
% Robustness:
%   - avoids inv()
%   - jittered Cholesky for PD enforcement
%   - strict dimension checks
%   - clamps to bounds (pragmatic safeguard)

% ---- Guard: only meaningful in hierarchical mode ----
if isfield(prior, 'use_hierarchical') && ~prior.use_hierarchical
    error('c6_gibbs_update_mu_theta called with prior.use_hierarchical = false.');
end

% Force column vector
theta = theta(:);
n_params = length(theta);

% Bounds (required)
[lb, ub] = get_bounds(prior);
lb = lb(:); ub = ub(:);

if length(lb) ~= n_params || length(ub) ~= n_params
    error('Bounds size mismatch: lb/ub length must match theta length (%d).', n_params);
end
if any(~isfinite(lb)) || any(~isfinite(ub)) || any(lb >= ub)
    error('Invalid bounds: check prior.hard_min/hard_max (finite and lb < ub).');
end

% Sigma_theta required
if isempty(Sigma_theta)
    error('Sigma_theta is empty. Ensure hierarchical initialization sets model.Sigma_theta.');
end
if ~isequal(size(Sigma_theta), [n_params, n_params])
    error('Sigma_theta must be %dx%d.', n_params, n_params);
end

% -------------------- Prior for mu_theta: N(mu0, Sigma0) --------------------
if isfield(prior, 'hyper') && isfield(prior.hyper, 'mu0') && isfield(prior.hyper, 'Sigma0')
    mu0 = prior.hyper.mu0(:);
    Sigma0 = prior.hyper.Sigma0;
elseif isfield(prior, 'mu_theta_init')
    mu0 = prior.mu_theta_init(:);
    if isfield(prior, 'prior_std')
        ps = prior.prior_std(:);
        if length(ps) ~= n_params
            error('prior.prior_std length must match theta length (%d).', n_params);
        end
        Sigma0 = diag(ps.^2);
    else
        Sigma0 = diag(((ub - lb) / 4).^2);
    end
else
    mu0 = lb + 0.5 * (ub - lb);
    Sigma0 = diag(((ub - lb) / 4).^2);
end

if numel(mu0) ~= n_params
    error('mu0 length (%d) does not match n_params (%d).', numel(mu0), n_params);
end
if ~isequal(size(Sigma0), [n_params, n_params])
    error('Sigma0 must be %dx%d.', n_params, n_params);
end

% Symmetrize covariances (protect against numerical asymmetry)
Sigma0     = 0.5 * (Sigma0 + Sigma0');
Sigma_theta = 0.5 * (Sigma_theta + Sigma_theta');

% -------------------- Stable posterior update using precision --------------------
% Posterior precision:
%   P = Sigma_post^{-1} = Sigma0^{-1} + Sigma_theta^{-1}
% Posterior mean:
%   mu_post = P^{-1} * (Sigma0^{-1}*mu0 + Sigma_theta^{-1}*theta)

R0 = safe_chol(Sigma0, 'Sigma0');            % Sigma0 = R0*R0'
Rt = safe_chol(Sigma_theta, 'Sigma_theta');  % Sigma_theta = Rt*Rt'

% Compute Sigma0^{-1}*mu0 and Sigma_theta^{-1}*theta via triangular solves
v0 = R0' \ (R0 \ mu0);
vt = Rt' \ (Rt \ theta);

% Build precisions explicitly (n_params small; stable via solves)
I = eye(n_params);
Sigma0_inv = R0' \ (R0 \ I);
SigmaT_inv = Rt' \ (Rt \ I);

P = Sigma0_inv + SigmaT_inv;
P = 0.5 * (P + P');  % symmetrize

RP = safe_chol(P, 'Posterior precision P'); % P = RP*RP'

b = v0 + vt;

% Solve P * mu_post = b
mu_post = RP' \ (RP \ b);

% Sample from N(mu_post, Sigma_post) where Sigma_post = P^{-1}
z = randn(n_params, 1);
y = RP \ z;                     % Cov(y) = P^{-1}

mu_theta_new = mu_post + y;

% Clamp to bounds (pragmatic safeguard)
mu_theta_new = max(mu_theta_new, lb);
mu_theta_new = min(mu_theta_new, ub);

% Final sanity check
if any(~isfinite(mu_theta_new))
    error('mu_theta_new contains non-finite values. Check Sigma_theta / Sigma0 conditioning.');
end

end

function R = safe_chol(A, name)
% Robust Cholesky with diagonal jitter
A = 0.5 * (A + A');
n = size(A, 1);

jitter = 0;
max_tries = 10;

for k = 1:max_tries
    [R, p] = chol(A + jitter * eye(n), 'lower');
    if p == 0
        return;
    end
    if k == 1
        scale = trace(A) / max(1, n);
        if ~isfinite(scale) || scale <= 0
            scale = 1;
        end
        jitter = 1e-10 * scale;
    else
        jitter = jitter * 10;
    end
end

error('Cholesky failed for %s after jitter up to %g.', name, jitter);
end

function [lb, ub] = get_bounds(prior)
if isfield(prior, 'bounds')
    lb = prior.bounds.lower;
    ub = prior.bounds.upper;
elseif isfield(prior, 'hard_min')
    lb = prior.hard_min;
    ub = prior.hard_max;
else
    error('No bounds found in prior structure');
end
end
