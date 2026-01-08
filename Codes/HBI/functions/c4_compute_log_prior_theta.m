function [log_prior, is_valid] = c4_compute_log_prior_theta(theta, mu_theta, Sigma_theta, prior)
% Compute log prior probability of theta
%
% Inputs:
%   theta       - model parameters [11 x 1]
%   mu_theta    - mean of parameter distribution [11 x 1] (ignored if use_hierarchical=false)
%   Sigma_theta - covariance matrix [11 x 11] (ignored if use_hierarchical=false)
%   prior       - structure from c1_define_hbi_prior
%
% Outputs:
%   log_prior - log prior probability
%   is_valid  - boolean, true if theta satisfies all constraints

is_valid = check_constraints(theta, prior);

if ~is_valid
    log_prior = -Inf;
    return;
end

if prior.use_hierarchical
    % Gaussian prior: N(theta | mu_theta, Sigma_theta)
    n = length(theta);
    diff = theta - mu_theta;
    
    [L, p] = chol(Sigma_theta, 'lower');
    if p > 0
        log_prior = -Inf;
        is_valid = false;
        return;
    end
    
    log_det = 2 * sum(log(diag(L)));
    v = L \ diff;
    mahal_sq = sum(v.^2);
    log_prior = -0.5 * n * log(2 * pi) - 0.5 * log_det - 0.5 * mahal_sq;
else
    % Uniform prior: flat within bounds
    log_prior = 0;
end

end

function is_valid = check_constraints(theta, prior)
% Check all constraints on theta
%
% Constraints:
%   1. Within hard bounds
%   2. Velocity monotonicity: Vs1 <= Vs2 <= Vs3 <= Vs4 <= Vs5 <= Vs6
%   3. Positivity of thicknesses

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